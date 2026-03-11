use std::f64::consts::PI;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use crate::atoms::Configuration;
use crate::cells::CellList;
use crate::constraints::{Constraints, PrecomputedConstraints};
use crate::log_println;
use crate::potential::PotentialSet;
use crate::rdf::compute_histograms;
use crate::xray::form_factor;

/// Estimate per-point sigma from data using windowed second finite differences.
///
/// The second difference `d2[i] = y[i+1] - 2*y[i] + y[i-1]` removes smooth
/// signal, leaving noise.  For white noise with std σ, `Var(d2) = 6σ²`.
/// A sliding window of half-width `half_w` gives a local noise estimate.
/// A minimum floor is applied so sigma never drops to zero.
pub fn estimate_sigma(y: &[f64], half_w: usize) -> Vec<f64> {
    let n = y.len();
    if n < 3 {
        let fallback = y.iter().map(|v| v.abs()).sum::<f64>() / n.max(1) as f64 * 0.01;
        return vec![fallback.max(1e-6); n];
    }

    // Second finite differences
    let mut d2 = vec![0.0; n];
    for i in 1..n - 1 {
        d2[i] = y[i + 1] - 2.0 * y[i] + y[i - 1];
    }
    // Pad edges
    d2[0] = d2[1];
    d2[n - 1] = d2[n - 2];

    // Windowed RMS of d2, scaled by 1/sqrt(6)
    let inv_sqrt6 = 1.0 / 6.0_f64.sqrt();
    let mut sigma = vec![0.0; n];
    for i in 0..n {
        let lo = i.saturating_sub(half_w);
        let hi = if i + half_w < n { i + half_w } else { n - 1 };
        let count = (hi - lo + 1) as f64;
        let rms = (d2[lo..=hi].iter().map(|v| v * v).sum::<f64>() / count).sqrt();
        sigma[i] = rms * inv_sqrt6;
    }

    // Floor: at least 1% of the global RMS of d2
    let global_rms = (d2.iter().map(|v| v * v).sum::<f64>() / n as f64).sqrt() * inv_sqrt6;
    let floor = (global_rms * 0.01).max(1e-6);
    for s in &mut sigma {
        if *s < floor {
            *s = floor;
        }
    }

    sigma
}

/// Serialisable RMC state (for checkpointing).
#[derive(Debug, Clone)]
pub struct RmcState {
    pub move_count: u64,
    pub accepted: u64,
    pub chi2: f64,
    pub max_step: f64,
    pub seed: u64,
}

/// Experimental dataset to fit against.
pub struct ExperimentalData {
    pub q: Vec<f64>,
    pub sq: Vec<f64>,
    pub sigma: Vec<f64>,
    pub weight: f64,
    pub kind: DataKind,
    pub fit_min: f64,
    pub fit_max: f64,
}

pub enum DataKind {
    Xray,
    Neutron,
}

/// Experimental g(r) dataset to fit against (real-space target).
pub struct ExperimentalGrData {
    pub r: Vec<f64>,
    pub gr: Vec<f64>,
    pub sigma: Vec<f64>,
    pub weight: f64,
    pub fit_min: f64,
    pub fit_max: f64,
    /// Q_max used when deriving the experimental g(r) from S(Q).
    pub qmax: f64,
    /// Apply Lorch modification W(Q) = sin(πQ/Qmax)/(πQ/Qmax) in inverse FT.
    pub lorch: bool,
}

/// RMC refinement parameters.
pub struct RmcParams {
    pub max_moves: u64,
    pub max_step: f64,
    pub checkpoint_every: u64,
    pub seed: u64,
    pub rdf_cutoff: f64,
    pub rdf_nbins: usize,
    pub q_grid: Vec<f64>,
    pub lorch: bool,
    pub print_every: u64,
    pub target_acceptance: f64,
    pub adjust_step_every: u64,
    pub anneal_start: f64,
    pub anneal_end: f64,
    pub anneal_steps: u64,
    pub convergence_threshold: f64,
    pub convergence_window: u64,
}

impl Default for RmcParams {
    fn default() -> Self {
        let nq = 500;
        let qmax = 20.0;
        let dq = qmax / nq as f64;
        RmcParams {
            max_moves: 1_000_000,
            max_step: 0.1,
            checkpoint_every: 50_000,
            seed: 42,
            rdf_cutoff: 10.0,
            rdf_nbins: 500,
            q_grid: (0..nq).map(|i| 0.3 + i as f64 * dq).collect(),
            lorch: true,
            print_every: 1000,
            target_acceptance: 0.3,
            adjust_step_every: 5000,
            anneal_start: 1.0,
            anneal_end: 1.0,
            anneal_steps: 0, // 0 = use max_moves
            convergence_threshold: 0.0,
            convergence_window: 50_000,
        }
    }
}

pub type StatusCallback = Box<dyn Fn(u64, u64, f64, f64, f64)>;

/// Single-threaded atom histogram using cell list for O(N_neighbors) scaling.
/// Writes into pre-allocated buffer `hist` (must be zeroed by caller).
#[inline]
fn atom_histogram_st(
    config: &Configuration,
    atom_idx: usize,
    pos: &[f64; 3],
    nbins: usize,
    cutoff2: f64,
    inv_dr: f64,
    n_species: usize,
    hist: &mut [f64],
    cell_list: &CellList,
    pos_cell: usize,
) {
    let ti = config.atoms[atom_idx].type_id;
    let box_lengths = &config.box_lengths;

    let neighbor_cells = cell_list.neighbor_cells(pos_cell);
    for &nc in &neighbor_cells {
        for j in cell_list.atoms_in_cell(nc) {
            if j == atom_idx {
                continue;
            }

            let pj = &config.atoms[j].position;
            let mut r2 = 0.0f64;
            for d in 0..3 {
                let mut delta = pj[d] - pos[d];
                let l = box_lengths[d];
                delta -= l * (delta / l).round();
                r2 += delta * delta;
            }

            if r2 < cutoff2 {
                let r = r2.sqrt();
                let bin = (r * inv_dr) as usize;
                if bin < nbins {
                    let tj = config.atoms[j].type_id;
                    let (a, b) = if ti <= tj { (ti, tj) } else { (tj, ti) };
                    let pair_idx = a * n_species - a * (a + 1) / 2 + b;
                    hist[pair_idx * nbins + bin] += 1.0;
                }
            }
        }
    }
}

/// Check min distance constraints using cell list (O(N_neighbors) instead of O(N)).
#[inline]
fn check_min_distances_cell(
    config: &Configuration,
    atom_idx: usize,
    new_pos: &[f64; 3],
    pc: &PrecomputedConstraints,
    cell_list: &CellList,
    pos_cell: usize,
) -> bool {
    let ti = config.atoms[atom_idx].type_id;
    let box_lengths = &config.box_lengths;

    let neighbor_cells = cell_list.neighbor_cells(pos_cell);
    for &nc in &neighbor_cells {
        for j in cell_list.atoms_in_cell(nc) {
            if j == atom_idx {
                continue;
            }

            let min_d2 = pc.min_dist_sq_for(ti, config.atoms[j].type_id);
            if min_d2 <= 0.0 {
                continue;
            }

            let pj = &config.atoms[j].position;
            let mut r2 = 0.0f64;
            for d in 0..3 {
                let mut delta = pj[d] - new_pos[d];
                let l = box_lengths[d];
                delta -= l * (delta / l).round();
                r2 += delta * delta;
            }

            if r2 < min_d2 {
                return false;
            }
        }
    }
    true
}

/// Check coordination constraints using cell list.
#[inline]
fn check_coordination_cell(
    config: &Configuration,
    atom_idx: usize,
    new_pos: &[f64; 3],
    pc: &PrecomputedConstraints,
    cell_list: &CellList,
) -> bool {
    let moved_type = config.atoms[atom_idx].type_id;
    let box_lengths = &config.box_lengths;

    for cc in &pc.coordination {
        if moved_type != cc.type_a && moved_type != cc.type_b {
            continue;
        }

        // We need to check atoms of type_a near the moved atom.
        // Use cell list centered on new_pos to find affected atoms.
        let new_cell = cell_list.cell_of[atom_idx]; // atom already has old cell; use new_pos cell
        let new_pos_cell = {
            let cx =
                ((new_pos[0] / box_lengths[0]).fract() * cell_list.nc[0] as f64).floor() as usize;
            let cy =
                ((new_pos[1] / box_lengths[1]).fract() * cell_list.nc[1] as f64).floor() as usize;
            let cz =
                ((new_pos[2] / box_lengths[2]).fract() * cell_list.nc[2] as f64).floor() as usize;
            let cx = cx.min(cell_list.nc[0] - 1);
            let cy = cy.min(cell_list.nc[1] - 1);
            let cz = cz.min(cell_list.nc[2] - 1);
            cz * cell_list.nc[0] * cell_list.nc[1] + cy * cell_list.nc[0] + cx
        };

        // Collect atoms of type sp_a that could be affected (within 2*cutoff of moved atom)
        // We need to check all type_a atoms, but only those near old or new position
        let check_cutoff2 = cc.cutoff2 * 4.0; // 2*cutoff squared

        // Check all type_a atoms near old or new position
        let old_cell = new_cell; // cell_list still has old assignment
        let old_neighbors = cell_list.neighbor_cells(old_cell);
        let new_neighbors = cell_list.neighbor_cells(new_pos_cell);

        // Merge neighbor cells (may have duplicates, that's fine)
        let mut atoms_to_check: Vec<usize> = Vec::new();
        // Always check the moved atom itself if it's type_a
        if moved_type == cc.type_a {
            atoms_to_check.push(atom_idx);
        }

        for &nc in old_neighbors.iter().chain(new_neighbors.iter()) {
            for j in cell_list.atoms_in_cell(nc) {
                if j == atom_idx {
                    continue;
                }
                if config.atoms[j].type_id != cc.type_a {
                    continue;
                }

                let pj = &config.atoms[j].position;
                let old_pos = &config.atoms[atom_idx].position;
                let old_r2 = min_image_r2_inline(old_pos, pj, box_lengths);
                let new_r2 = min_image_r2_inline(new_pos, pj, box_lengths);
                if old_r2 < check_cutoff2 || new_r2 < check_cutoff2 {
                    atoms_to_check.push(j);
                }
            }
        }

        // Deduplicate
        atoms_to_check.sort_unstable();
        atoms_to_check.dedup();

        // For each atom of type_a to check, count type_b neighbors
        for &i in &atoms_to_check {
            let pos_i = if i == atom_idx {
                *new_pos
            } else {
                config.atoms[i].position
            };

            // Find cell for pos_i
            let ci = {
                let cx =
                    ((pos_i[0] / box_lengths[0]).fract() * cell_list.nc[0] as f64).floor() as usize;
                let cy =
                    ((pos_i[1] / box_lengths[1]).fract() * cell_list.nc[1] as f64).floor() as usize;
                let cz =
                    ((pos_i[2] / box_lengths[2]).fract() * cell_list.nc[2] as f64).floor() as usize;
                let cx = cx.min(cell_list.nc[0] - 1);
                let cy = cy.min(cell_list.nc[1] - 1);
                let cz = cz.min(cell_list.nc[2] - 1);
                cz * cell_list.nc[0] * cell_list.nc[1] + cy * cell_list.nc[0] + cx
            };

            let mut count = 0usize;
            let nbr_cells = cell_list.neighbor_cells(ci);
            for &nc in &nbr_cells {
                for j in cell_list.atoms_in_cell(nc) {
                    if j == i {
                        continue;
                    }
                    if config.atoms[j].type_id != cc.type_b {
                        continue;
                    }

                    let pos_j = if j == atom_idx {
                        *new_pos
                    } else {
                        config.atoms[j].position
                    };
                    let r2 = min_image_r2_inline(&pos_i, &pos_j, box_lengths);
                    if r2 < cc.cutoff2 {
                        count += 1;
                    }
                }
            }

            if count < cc.min || count > cc.max {
                return false;
            }
        }
    }
    true
}

#[inline]
fn min_image_r2_inline(a: &[f64; 3], b: &[f64; 3], box_lengths: &[f64; 3]) -> f64 {
    let mut r2 = 0.0f64;
    for d in 0..3 {
        let mut delta = b[d] - a[d];
        let l = box_lengths[d];
        delta -= l * (delta / l).round();
        r2 += delta * delta;
    }
    r2
}

/// Run the RMC refinement loop with incremental S(Q) updates.
pub fn run_rmc(
    config: &mut Configuration,
    experiments: &[ExperimentalData],
    gr_data: &[ExperimentalGrData],
    constraints: &Constraints,
    params: &RmcParams,
    potential: Option<&PotentialSet>,
    checkpoint_fn: Option<Box<dyn Fn(&RmcState, &Configuration)>>,
    resume_state: Option<RmcState>,
) -> RmcState {
    // On resume, re-seed with seed + move_count (deterministic but different trajectory)
    let rng_seed = if let Some(ref rs) = resume_state {
        rs.seed.wrapping_add(rs.move_count)
    } else {
        params.seed
    };
    let mut rng = StdRng::seed_from_u64(rng_seed);
    let n_atoms = config.atoms.len();
    let nbins = params.rdf_nbins;
    let nq = params.q_grid.len();
    let dr = params.rdf_cutoff / nbins as f64;
    let inv_dr = 1.0 / dr;
    let rho0 = config.number_density();
    let n_types = config.species.len();
    let n_pairs = config.num_type_pairs();
    let volume = config.volume();

    // Precompute tables
    log_println!("Precomputing lookup tables...");
    let r_grid: Vec<f64> = (0..nbins).map(|i| (i as f64 + 0.5) * dr).collect();
    let r_max = r_grid[nbins - 1];

    // sin(Q_k * r_i), row-major [nbins][nq]
    let mut sin_table = vec![0.0f64; nbins * nq];
    for i in 0..nbins {
        let row = i * nq;
        for k in 0..nq {
            sin_table[row + k] = (params.q_grid[k] * r_grid[i]).sin();
        }
    }

    // Lorch window
    let lorch_w: Vec<f64> = (0..nbins)
        .map(|i| {
            let r = r_grid[i];
            if params.lorch && r > 0.0 {
                let arg = PI * r / r_max;
                arg.sin() / arg
            } else {
                1.0
            }
        })
        .collect();

    // r * W(r) with endpoint weight (midpoint rule, all weight 1)
    let rw: Vec<f64> = (0..nbins).map(|i| r_grid[i] * lorch_w[i]).collect();

    // Normalisation denominators per pair per bin
    // For unlike (a!=b): norm = N_a * rho_b * 4πr²dr (hist counted once per pair in i<j loop)
    //   g = hist / norm
    // For like (a==b): norm = N_a * rho_a * 4πr²dr / 2
    //   g = 2 * hist / (N_a * rho_a * 4πr²dr) = hist / (norm/2)... let's just store 1/norm.
    let mut inv_norm: Vec<f64> = vec![0.0; n_pairs * nbins]; // 1/norm so g = inv_norm * hist * like_factor
    let mut like_factor: Vec<f64> = vec![1.0; n_pairs]; // 2 for like pairs, 1 for unlike
    for a in 0..n_types {
        for b in a..n_types {
            let pair_idx = config.pair_index(a, b);
            let n_a = config.count_type(a) as f64;
            let rho_b = config.count_type(b) as f64 / volume;
            if a == b {
                like_factor[pair_idx] = 2.0;
            }
            let base = pair_idx * nbins;
            for i in 0..nbins {
                let r = r_grid[i];
                let shell = n_a * rho_b * 4.0 * PI * r * r * dr;
                if shell > 0.0 {
                    inv_norm[base + i] = 1.0 / shell;
                }
            }
        }
    }

    // X-ray weights: w_ab(Q_k)
    let conc: Vec<f64> = (0..n_types).map(|t| config.concentration(t)).collect();
    let form_factors: Vec<Vec<f64>> = config
        .species
        .iter()
        .map(|s| form_factor(s, &params.q_grid))
        .collect();

    let mut xray_w = vec![0.0f64; n_pairs * nq];
    for k in 0..nq {
        let f_avg: f64 = (0..n_types).map(|a| conc[a] * form_factors[a][k]).sum();
        let f_avg_sq = f_avg * f_avg;
        if f_avg_sq < 1e-30 {
            continue;
        }
        for a in 0..n_types {
            for b in a..n_types {
                let pair_idx = config.pair_index(a, b);
                let dab = if a == b { 1.0 } else { 2.0 };
                xray_w[pair_idx * nq + k] =
                    dab * conc[a] * conc[b] * form_factors[a][k] * form_factors[b][k] / f_avg_sq;
            }
        }
    }

    // Precompute experimental interpolation indices
    let mut exp_interp: Vec<Vec<(usize, f64, usize)>> = Vec::new(); // (lo, t, exp_idx)
    for exp in experiments.iter() {
        let mut interp = Vec::new();
        for (i, &qe) in exp.q.iter().enumerate() {
            if qe < params.q_grid[0] || qe > params.q_grid[nq - 1] {
                continue;
            }
            if qe < exp.fit_min || qe > exp.fit_max {
                continue;
            }
            let mut lo = 0;
            let mut hi = nq - 1;
            while hi - lo > 1 {
                let mid = (lo + hi) / 2;
                if params.q_grid[mid] <= qe {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let t = (qe - params.q_grid[lo]) / (params.q_grid[hi] - params.q_grid[lo]);
            interp.push((lo, t, i));
        }
        exp_interp.push(interp);
    }

    // --- Compute initial state ---
    log_println!("Computing initial RDF histograms...");
    let hist_map = compute_histograms(config, nbins, params.rdf_cutoff);
    let mut flat_hist = vec![0.0f64; n_pairs * nbins];
    for (&p, bins) in &hist_map {
        let start = p * nbins;
        flat_hist[start..start + nbins].copy_from_slice(bins);
    }

    // Compute initial partial S(Q) from histograms
    // S_ab(Q_k) = 1 + (4πρ₀*dr/Q_k) * Σ_i rw_i * [g_ab(r_i) - 1] * sin(Q_k*r_i)
    let prefactor_sq = 4.0 * PI * rho0 * dr;
    let inv_q: Vec<f64> = params
        .q_grid
        .iter()
        .map(|&q| if q > 0.05 { 1.0 / q } else { 0.0 })
        .collect();

    let mut partial_sq = vec![0.0f64; n_pairs * nq]; // stores S_ab(Q) - 1 (the integral part)

    for p in 0..n_pairs {
        let hist_base = p * nbins;
        let norm_base = p * nbins;
        let sq_base = p * nq;
        let lf = like_factor[p];

        for i in 0..nbins {
            let g = lf * flat_hist[hist_base + i] * inv_norm[norm_base + i];
            let contrib = rw[i] * (g - 1.0); // r * W * (g - 1)
            if contrib.abs() < 1e-30 {
                continue;
            }
            let sin_row = i * nq;
            for k in 0..nq {
                partial_sq[sq_base + k] += contrib * sin_table[sin_row + k];
            }
        }

        // Apply prefactor and add 1
        for k in 0..nq {
            partial_sq[sq_base + k] = 1.0 + prefactor_sq * partial_sq[sq_base + k] * inv_q[k];
        }
    }

    // Compute initial total X-ray S(Q) and chi2
    let mut total_sq = vec![0.0f64; nq];
    for k in 0..nq {
        let mut s = 0.0;
        for p in 0..n_pairs {
            s += xray_w[p * nq + k] * partial_sq[p * nq + k];
        }
        total_sq[k] = s;
    }

    let mut sq_chi2_current = 0.0;
    for (ei, exp) in experiments.iter().enumerate() {
        let mut chi2 = 0.0;
        for &(lo, t, i) in &exp_interp[ei] {
            let sq_calc = total_sq[lo] + t * (total_sq[lo + 1] - total_sq[lo]);
            let diff = sq_calc - exp.sq[i];
            chi2 += (diff * diff) / (exp.sigma[i] * exp.sigma[i]);
        }
        sq_chi2_current += exp.weight * chi2;
    }

    // --- g(r) precomputation via inverse Fourier transform of total S_X(Q) ---
    let has_gr = !gr_data.is_empty();
    let dq = if nq > 1 {
        params.q_grid[1] - params.q_grid[0]
    } else {
        1.0
    };

    // For each g(r) dataset: build FT matrix M[n_r × nq], compute initial g(r)
    let mut gr_ft_matrices: Vec<Vec<f64>> = Vec::new();
    let mut gr_cached: Vec<Vec<f64>> = Vec::new();
    let mut gr_scratch: Vec<Vec<f64>> = Vec::new();

    for gd in gr_data.iter() {
        let n_r = gd.r.len();
        // M[i][k] = dq * Q_k * [W(Q_k)] * sin(Q_k * r_i) / (2π²ρ₀ * r_i)
        // Only include Q_k <= qmax; optionally apply Lorch W(Q) = sin(πQ/Qmax)/(πQ/Qmax)
        let qmax_gr = gd.qmax;
        let use_lorch = gd.lorch;
        let mut matrix = vec![0.0f64; n_r * nq];
        for i in 0..n_r {
            let ri = gd.r[i];
            if ri < 1e-10 {
                continue;
            }
            let row = i * nq;
            let prefactor = dq / (2.0 * PI * PI * rho0 * ri);
            for k in 0..nq {
                let qk = params.q_grid[k];
                if qk > qmax_gr {
                    break;
                }
                let window = if use_lorch {
                    let arg = PI * qk / qmax_gr;
                    if arg > 1e-10 {
                        arg.sin() / arg
                    } else {
                        1.0
                    }
                } else {
                    1.0
                };
                matrix[row + k] = prefactor * qk * window * (qk * ri).sin();
            }
        }
        gr_ft_matrices.push(matrix);
        gr_cached.push(vec![0.0; n_r]);
        gr_scratch.push(vec![0.0; n_r]);
    }

    // Precompute g(r) fitting indices (only points within fit_min..fit_max)
    let mut gr_fit_indices: Vec<Vec<usize>> = Vec::new();
    for gd in gr_data.iter() {
        let indices: Vec<usize> = (0..gd.r.len())
            .filter(|&i| gd.r[i] >= gd.fit_min && gd.r[i] <= gd.fit_max)
            .collect();
        gr_fit_indices.push(indices);
    }

    if has_gr {
        for (di, gd) in gr_data.iter().enumerate() {
            log_println!(
                "g(r) dataset {}: {} total points, {} in fit range [{:.2}, {:.2}] A",
                di,
                gd.r.len(),
                gr_fit_indices[di].len(),
                gd.fit_min,
                gd.fit_max
            );
        }
    }

    // Compute initial model g(r) from initial total_sq
    let mut gr_chi2_current = 0.0;
    for (di, gd) in gr_data.iter().enumerate() {
        let n_r = gd.r.len();
        let matrix = &gr_ft_matrices[di];
        for i in 0..n_r {
            let row = i * nq;
            let mut val = 1.0;
            for k in 0..nq {
                val += matrix[row + k] * (total_sq[k] - 1.0);
            }
            gr_cached[di][i] = val;
        }
        let mut chi2 = 0.0;
        for &i in &gr_fit_indices[di] {
            let diff = gr_cached[di][i] - gd.gr[i];
            chi2 += (diff * diff) / (gd.sigma[i] * gd.sigma[i]);
        }
        gr_chi2_current += gd.weight * chi2;
    }

    let mut current_chi2 = sq_chi2_current + gr_chi2_current;

    if has_gr {
        log_println!(
            "Initial chi2 = {:.6} (sq: {:.6}, gr: {:.6})",
            current_chi2,
            sq_chi2_current,
            gr_chi2_current
        );
    } else {
        log_println!("Initial chi2 = {:.6}", current_chi2);
    }

    let mut state = if let Some(rs) = resume_state {
        log_println!(
            "Resuming from move {} (accepted {}, max_step {:.4})",
            rs.move_count,
            rs.accepted,
            rs.max_step
        );
        RmcState {
            move_count: rs.move_count,
            accepted: rs.accepted,
            chi2: current_chi2, // always recomputed from checkpoint configuration
            max_step: rs.max_step,
            seed: rs.seed,
        }
    } else {
        RmcState {
            move_count: 0,
            accepted: 0,
            chi2: current_chi2,
            max_step: params.max_step,
            seed: params.seed,
        }
    };

    let mut recent_accepted = 0u64;
    let mut recent_total = 0u64;

    // Scratch buffers for incremental updates
    let mut delta_partial_sq = vec![0.0f64; n_pairs * nq];
    let mut new_total_sq = vec![0.0f64; nq];
    let mut old_hist_buf = vec![0.0f64; n_pairs * nbins];
    let mut new_hist_buf = vec![0.0f64; n_pairs * nbins];
    let mut delta_sq_buf = vec![0.0f64; nq];

    // Precompute constraint lookup table (eliminates string allocs in hot loop)
    let pc = PrecomputedConstraints::from_constraints(constraints, config);
    let cutoff2 = params.rdf_cutoff * params.rdf_cutoff;

    // Build cell list for O(1)-per-neighbor spatial lookups
    let positions: Vec<[f64; 3]> = config.atoms.iter().map(|a| a.position).collect();
    let mut cell_list = CellList::new(&positions, &config.box_lengths, params.rdf_cutoff);
    log_println!(
        "Cell list: {}x{}x{} = {} cells (cell size: {:.2} A)",
        cell_list.nc[0],
        cell_list.nc[1],
        cell_list.nc[2],
        cell_list.n_cells,
        cell_list.cell_size[0]
    );

    // Potential energy initialization
    let energy_weight = potential.map_or(0.0, |p| p.weight);
    let mut current_energy = if let Some(pot) = potential {
        let e = pot.total_energy(config, &cell_list);
        log_println!(
            "Initial potential energy = {:.6} eV (weight = {:.6})",
            e,
            energy_weight
        );
        e
    } else {
        0.0
    };

    // Annealing setup
    let annealing = (params.anneal_start - params.anneal_end).abs() > 1e-10;
    let anneal_n = if params.anneal_steps > 0 {
        params.anneal_steps
    } else {
        params.max_moves
    };
    // Precompute log ratio to replace powf with exp (faster per-move)
    let anneal_log_ratio = if annealing {
        (params.anneal_end / params.anneal_start).ln()
    } else {
        0.0
    };
    if annealing {
        log_println!(
            "Simulated annealing: T = {:.2} -> {:.2} over {} moves ({:.0}% of run)",
            params.anneal_start,
            params.anneal_end,
            anneal_n,
            100.0 * anneal_n as f64 / params.max_moves as f64
        );
    }

    // Best-structure tracking: save atom positions at lowest chi2
    let mut best_chi2 = current_chi2;
    let mut best_positions: Vec<[f64; 3]> = config.atoms.iter().map(|a| a.position).collect();
    let mut best_move = 0u64;

    // Convergence detection (offset by resume start)
    let conv_active = params.convergence_threshold > 0.0;
    let mut conv_chi2 = current_chi2;
    let mut conv_next_check = state.move_count + params.convergence_window;

    // Calibration: accumulate |delta_chi2| and |delta_E| to suggest weight
    let calibration_moves: u64 = if potential.is_some() { 1000 } else { 0 };
    let mut calib_sum_dchi2 = 0.0f64;
    let mut calib_sum_de = 0.0f64;
    let mut calib_count = 0u64;

    // === Main RMC loop ===
    let start_move = state.move_count;
    for move_num in start_move..params.max_moves {
        let atom_idx = rng.gen_range(0..n_atoms);

        let dx: f64 = rng.gen_range(-state.max_step..state.max_step);
        let dy: f64 = rng.gen_range(-state.max_step..state.max_step);
        let dz: f64 = rng.gen_range(-state.max_step..state.max_step);

        let old_pos = config.atoms[atom_idx].position;
        let mut new_pos = [old_pos[0] + dx, old_pos[1] + dy, old_pos[2] + dz];
        for d in 0..3 {
            let l = config.box_lengths[d];
            new_pos[d] -= l * (new_pos[d] / l).floor();
        }

        // Determine cell for new position (for cell-list lookups)
        let new_pos_cell = {
            let cx = ((new_pos[0] / config.box_lengths[0]).fract() * cell_list.nc[0] as f64).floor()
                as usize;
            let cy = ((new_pos[1] / config.box_lengths[1]).fract() * cell_list.nc[1] as f64).floor()
                as usize;
            let cz = ((new_pos[2] / config.box_lengths[2]).fract() * cell_list.nc[2] as f64).floor()
                as usize;
            let cx = cx.min(cell_list.nc[0] - 1);
            let cy = cy.min(cell_list.nc[1] - 1);
            let cz = cz.min(cell_list.nc[2] - 1);
            cz * cell_list.nc[0] * cell_list.nc[1] + cy * cell_list.nc[0] + cx
        };

        // Check constraints using cell list (O(N_neighbors) instead of O(N))
        if !check_min_distances_cell(config, atom_idx, &new_pos, &pc, &cell_list, new_pos_cell) {
            recent_total += 1;
            state.move_count += 1;
            continue;
        }
        if !pc.coordination.is_empty()
            && !check_coordination_cell(config, atom_idx, &new_pos, &pc, &cell_list)
        {
            recent_total += 1;
            state.move_count += 1;
            continue;
        }

        // --- Incremental histogram + S(Q) update ---
        // Compute old histogram contributions using cell list
        let old_pos_cell = cell_list.cell_of[atom_idx];
        old_hist_buf.fill(0.0);
        atom_histogram_st(
            config,
            atom_idx,
            &old_pos,
            nbins,
            cutoff2,
            inv_dr,
            n_types,
            &mut old_hist_buf,
            &cell_list,
            old_pos_cell,
        );

        // Move atom temporarily
        config.atoms[atom_idx].position = new_pos;

        // Compute new histogram contributions
        new_hist_buf.fill(0.0);
        atom_histogram_st(
            config,
            atom_idx,
            &new_pos,
            nbins,
            cutoff2,
            inv_dr,
            n_types,
            &mut new_hist_buf,
            &cell_list,
            new_pos_cell,
        );

        // Compute ΔS_ab(Q) from histogram delta
        // ΔS_ab(Q_k) = (4πρ₀*dr/Q_k) * Σ_i rw_i * Δg_ab(r_i) * sin(Q_k*r_i)
        delta_partial_sq.fill(0.0);

        for p in 0..n_pairs {
            let hist_base = p * nbins;
            let norm_base = p * nbins;
            let sq_base = p * nq;
            let lf = like_factor[p];

            for i in 0..nbins {
                let dh = new_hist_buf[hist_base + i] - old_hist_buf[hist_base + i];
                if dh == 0.0 {
                    continue;
                }
                let dg = lf * dh * inv_norm[norm_base + i];
                let contrib = rw[i] * dg;
                let sin_row = i * nq;
                for k in 0..nq {
                    delta_partial_sq[sq_base + k] += contrib * sin_table[sin_row + k];
                }
            }

            // Apply prefactor
            for k in 0..nq {
                delta_partial_sq[sq_base + k] *= prefactor_sq * inv_q[k];
            }
        }

        // Compute new total S_X(Q) = Σ w_ab * (S_ab + ΔS_ab)
        for k in 0..nq {
            let mut s = 0.0;
            for p in 0..n_pairs {
                let idx = p * nq + k;
                s += xray_w[idx] * (partial_sq[idx] + delta_partial_sq[idx]);
            }
            new_total_sq[k] = s;
        }

        // Compute new S(Q) chi2
        let mut new_sq_chi2 = 0.0;
        for (ei, exp) in experiments.iter().enumerate() {
            let mut chi2 = 0.0;
            for &(lo, t, i) in &exp_interp[ei] {
                let sq_calc = new_total_sq[lo] + t * (new_total_sq[lo + 1] - new_total_sq[lo]);
                let diff = sq_calc - exp.sq[i];
                chi2 += (diff * diff) / (exp.sigma[i] * exp.sigma[i]);
            }
            new_sq_chi2 += exp.weight * chi2;
        }

        // Compute new g(r) incrementally and its chi2
        // Precompute ΔS(Q) vector once (enables auto-vectorization of matvec)
        let mut new_gr_chi2 = 0.0;
        if has_gr {
            for k in 0..nq {
                delta_sq_buf[k] = new_total_sq[k] - total_sq[k];
            }
            for (di, gd) in gr_data.iter().enumerate() {
                let n_r = gd.r.len();
                let matrix = &gr_ft_matrices[di];
                for i in 0..n_r {
                    let row = i * nq;
                    let mut delta = 0.0;
                    for k in 0..nq {
                        delta += matrix[row + k] * delta_sq_buf[k];
                    }
                    gr_scratch[di][i] = gr_cached[di][i] + delta;
                }
                let mut chi2 = 0.0;
                for &i in &gr_fit_indices[di] {
                    let diff = gr_scratch[di][i] - gd.gr[i];
                    chi2 += (diff * diff) / (gd.sigma[i] * gd.sigma[i]);
                }
                new_gr_chi2 += gd.weight * chi2;
            }
        }

        let new_chi2 = new_sq_chi2 + new_gr_chi2;

        // Compute potential energy delta if potential is active
        let delta_energy = if let Some(pot) = potential {
            let old_e = pot.energy_of_atom(config, atom_idx, &old_pos, &cell_list, old_pos_cell);
            let new_e = pot.energy_of_atom(config, atom_idx, &new_pos, &cell_list, new_pos_cell);
            new_e - old_e
        } else {
            0.0
        };

        // Calibration accumulation
        if calib_count < calibration_moves {
            calib_sum_dchi2 += (new_chi2 - current_chi2).abs();
            calib_sum_de += delta_energy.abs();
            calib_count += 1;
            if calib_count == calibration_moves {
                let avg_dchi2 = calib_sum_dchi2 / calib_count as f64;
                let avg_de = calib_sum_de / calib_count as f64;
                let suggested = if avg_de > 1e-15 {
                    avg_dchi2 / avg_de
                } else {
                    0.0
                };
                log_println!(
                    "Calibration ({} moves): avg |delta_chi2| = {:.6}, avg |delta_E| = {:.6} eV",
                    calib_count,
                    avg_dchi2,
                    avg_de
                );
                log_println!(
                    "  Current weight = {:.6}, suggested weight for equal balance = {:.6}",
                    energy_weight,
                    suggested
                );
                log_println!(
                    "  Ratio current/suggested = {:.4}",
                    energy_weight / suggested.max(1e-30)
                );
            }
        }

        // Metropolis acceptance (with optional annealing temperature)
        // Exponential schedule: T(n) = T_start * (T_end/T_start)^(n/anneal_n)
        // After anneal_n moves, T stays at anneal_end
        let temperature = if annealing {
            let frac = (move_num as f64 / anneal_n as f64).min(1.0);
            params.anneal_start * (frac * anneal_log_ratio).exp()
        } else {
            1.0
        };
        let delta_chi2 = new_chi2 - current_chi2;
        let delta_cost = delta_chi2 + energy_weight * delta_energy;
        let accept = if delta_cost < 0.0 {
            true
        } else {
            let prob = (-delta_cost / (2.0 * temperature)).exp();
            rng.gen::<f64>() < prob
        };

        if accept {
            current_chi2 = new_chi2;
            sq_chi2_current = new_sq_chi2;
            gr_chi2_current = new_gr_chi2;
            current_energy += delta_energy;
            // Commit: update histograms, partial S(Q), total S(Q)
            for i in 0..n_pairs * nbins {
                flat_hist[i] += new_hist_buf[i] - old_hist_buf[i];
            }
            for i in 0..n_pairs * nq {
                partial_sq[i] += delta_partial_sq[i];
            }
            total_sq.copy_from_slice(&new_total_sq);
            // Commit g(r) updates
            for (di, _) in gr_data.iter().enumerate() {
                std::mem::swap(&mut gr_cached[di], &mut gr_scratch[di]);
            }
            // Update cell list for accepted move
            cell_list.move_atom(atom_idx, &new_pos);
            state.accepted += 1;
            recent_accepted += 1;

            // Track best structure
            if new_chi2 < best_chi2 {
                best_chi2 = new_chi2;
                best_move = move_num + 1;
                for (i, atom) in config.atoms.iter().enumerate() {
                    best_positions[i] = atom.position;
                }
            }
        } else {
            // Revert atom position
            config.atoms[atom_idx].position = old_pos;
        }

        state.move_count = move_num + 1;
        state.chi2 = current_chi2;
        recent_total += 1;

        // Print status
        if state.move_count % params.print_every == 0 {
            let ratio = if recent_total > 0 {
                recent_accepted as f64 / recent_total as f64
            } else {
                0.0
            };
            let overall_ratio = state.accepted as f64 / state.move_count as f64;
            let chi2_str = if potential.is_some() && has_gr {
                let cost = current_chi2 + energy_weight * current_energy;
                format!(
                    "cost = {:.4} (chi2: {:.4} [sq: {:.4}, gr: {:.4}], w*E: {:.4} [E: {:.2}])",
                    cost,
                    current_chi2,
                    sq_chi2_current,
                    gr_chi2_current,
                    energy_weight * current_energy,
                    current_energy
                )
            } else if potential.is_some() {
                let cost = current_chi2 + energy_weight * current_energy;
                format!(
                    "cost = {:.4} (chi2: {:.4}, w*E: {:.4} [E: {:.2}])",
                    cost,
                    current_chi2,
                    energy_weight * current_energy,
                    current_energy
                )
            } else if has_gr {
                format!(
                    "chi2 = {:.6} (sq: {:.6}, gr: {:.6})",
                    current_chi2, sq_chi2_current, gr_chi2_current
                )
            } else {
                format!("chi2 = {:.6}", current_chi2)
            };
            let temp_str = if annealing {
                format!(", T = {:.3}", temperature)
            } else {
                String::new()
            };
            log_println!(
                "Move {}/{}: {}, accept = {:.3} (recent {:.3}), step = {:.4} A{}",
                state.move_count,
                params.max_moves,
                chi2_str,
                overall_ratio,
                ratio,
                state.max_step,
                temp_str
            );
        }

        // Adaptive step size
        if state.move_count % params.adjust_step_every == 0 && recent_total > 0 {
            if annealing && move_num < anneal_n {
                // During annealing: scale step with T, but never above initial max_step
                state.max_step = params.max_step * temperature.clamp(params.anneal_end, 1.0);
            } else {
                // After annealing (or no annealing): adapt based on acceptance rate
                let ratio = recent_accepted as f64 / recent_total as f64;
                if ratio > params.target_acceptance + 0.05 {
                    state.max_step *= 1.05;
                } else if ratio < params.target_acceptance - 0.05 {
                    state.max_step *= 0.95;
                }
            }
            state.max_step = state.max_step.clamp(0.001, 2.0);
            recent_accepted = 0;
            recent_total = 0;
        }

        // Checkpoint
        if state.move_count % params.checkpoint_every == 0 {
            if let Some(ref f) = checkpoint_fn {
                f(&state, config);
            }
        }

        // Convergence check (skip while annealing is still active)
        if conv_active
            && state.move_count >= conv_next_check
            && temperature <= params.anneal_end * 1.01
        {
            let improvement = conv_chi2 - current_chi2;
            if improvement < params.convergence_threshold {
                log_println!(
                    "\nConverged at move {}: chi2 improved by {:.6} over last {} moves (threshold: {:.6})",
                    state.move_count, improvement, params.convergence_window, params.convergence_threshold
                );
                break;
            }
            conv_chi2 = current_chi2;
            conv_next_check += params.convergence_window;
        } else if conv_active
            && state.move_count >= conv_next_check
            && temperature > params.anneal_end * 1.01
        {
            // Reset baseline during annealing so first post-anneal check has a fresh reference
            conv_chi2 = current_chi2;
            conv_next_check += params.convergence_window;
        }
    }

    // Restore best structure
    if best_chi2 < state.chi2 {
        log_println!(
            "\nRestoring best structure from move {} (chi2 = {:.6})",
            best_move,
            best_chi2
        );
        for (i, atom) in config.atoms.iter_mut().enumerate() {
            atom.position = best_positions[i];
        }
        state.chi2 = best_chi2;
    }

    log_println!("\nRMC refinement complete.");
    log_println!(
        "Final chi2 = {:.6} (best at move {}), accepted {}/{} ({:.1}%)",
        state.chi2,
        best_move,
        state.accepted,
        state.move_count,
        100.0 * state.accepted as f64 / state.move_count.max(1) as f64
    );

    state
}
