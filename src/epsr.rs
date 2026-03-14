//! Empirical Potential Structure Refinement (EPSR).
//!
//! Implements the Soper (1996) algorithm: an outer loop that iteratively
//! refines a perturbation potential (EP) so the MC simulation naturally
//! reproduces the experimental structure factor data.

use std::f64::consts::PI;
use std::path::Path;

use crate::atoms::Configuration;
use crate::neutron::scattering_length as neutron_scattering_length;
use crate::potential::{PairPotential, PotentialSet};
use crate::xray::form_factor;

/// Cumulative empirical potential state for all pairs.
pub struct EpsrState {
    /// EP tables per pair, flat [n_pairs][n_bins].
    pub ep_tables: Vec<Vec<f64>>,
    /// Grid spacing in Å.
    pub dr: f64,
    /// Number of grid bins.
    pub n_bins: usize,
    /// EP cutoff in Å.
    pub cutoff: f64,
    /// Pair labels (e.g. "Ca-O").
    pub pair_labels: Vec<String>,
    /// (type_a, type_b) for each pair.
    pub pair_types: Vec<(usize, usize)>,
    /// Number of species.
    pub n_types: usize,
    /// Number of unique pairs.
    pub n_pairs: usize,
}

impl EpsrState {
    /// Initialize with zero EP for all pairs.
    pub fn new(species: &[String], cutoff: f64, dr: f64) -> Self {
        let n_types = species.len();
        let n_bins = (cutoff / dr).ceil() as usize + 1;
        let mut pair_labels = Vec::new();
        let mut pair_types = Vec::new();
        let mut ep_tables = Vec::new();

        for a in 0..n_types {
            for b in a..n_types {
                pair_labels.push(format!("{}-{}", species[a], species[b]));
                pair_types.push((a, b));
                ep_tables.push(vec![0.0; n_bins]);
            }
        }

        let n_pairs = pair_labels.len();
        EpsrState {
            ep_tables,
            dr,
            n_bins,
            cutoff,
            pair_labels,
            pair_types,
            n_types,
            n_pairs,
        }
    }

    /// Load previous EP tables from `epsr_ep_{pair}.dat` files in the given directory.
    ///
    /// Each file is a two-column (r, V) table. The data is interpolated onto the
    /// current EP grid and added to the existing (zero-initialized) tables.
    /// Missing files are silently skipped (that pair keeps zero EP).
    pub fn load_potentials(&mut self, dir: &Path) -> Result<usize, Box<dyn std::error::Error>> {
        let mut loaded = 0usize;
        for p in 0..self.n_pairs {
            let filename = format!("epsr_ep_{}.dat", self.pair_labels[p]);
            let path = dir.join(&filename);
            if !path.exists() {
                continue;
            }

            // Read two-column file, skipping comment lines
            let content = std::fs::read_to_string(&path)?;
            let mut r_data = Vec::new();
            let mut v_data = Vec::new();
            for line in content.lines() {
                let line = line.trim();
                if line.is_empty() || line.starts_with('#') {
                    continue;
                }
                let mut parts = line.split_whitespace();
                if let (Some(r_str), Some(v_str)) = (parts.next(), parts.next()) {
                    if let (Ok(r), Ok(v)) = (r_str.parse::<f64>(), v_str.parse::<f64>()) {
                        r_data.push(r);
                        v_data.push(v);
                    }
                }
            }

            if r_data.len() < 2 {
                continue;
            }

            // Interpolate onto our grid
            for i in 0..self.n_bins {
                let r = i as f64 * self.dr;
                self.ep_tables[p][i] = interp_linear(&r_data, &v_data, r);
            }
            loaded += 1;
        }
        Ok(loaded)
    }

    /// Compute X-ray weights w_ab(Q) for given species/concentrations/form factors.
    ///
    /// Returns flat array [n_pairs * nq] with w_ab(Q_k) = δ_ab c_a c_b f_a f_b / <f>².
    pub fn compute_xray_weights(config: &Configuration, q_grid: &[f64]) -> Vec<f64> {
        let n_types = config.species.len();
        let nq = q_grid.len();
        let n_pairs = n_types * (n_types + 1) / 2;
        let conc: Vec<f64> = (0..n_types).map(|t| config.concentration(t)).collect();
        let ff: Vec<Vec<f64>> = config
            .species
            .iter()
            .map(|s| form_factor(s, q_grid))
            .collect();

        let mut weights = vec![0.0f64; n_pairs * nq];
        for k in 0..nq {
            let f_avg: f64 = (0..n_types).map(|a| conc[a] * ff[a][k]).sum();
            let f_avg_sq = f_avg * f_avg;
            if f_avg_sq < 1e-30 {
                continue;
            }
            for a in 0..n_types {
                for b in a..n_types {
                    let pair_idx = pair_index(a, b, n_types);
                    let dab = if a == b { 1.0 } else { 2.0 };
                    weights[pair_idx * nq + k] =
                        dab * conc[a] * conc[b] * ff[a][k] * ff[b][k] / f_avg_sq;
                }
            }
        }
        weights
    }

    /// Compute Faber-Ziman neutron weights w_ab(Q) (Q-independent).
    ///
    /// w_ab = (2-δ_ab) * c_a * c_b * b_a * b_b / <b>²
    pub fn compute_neutron_weights(config: &Configuration, nq: usize) -> Vec<f64> {
        let n_types = config.species.len();
        let n_pairs = n_types * (n_types + 1) / 2;
        let conc: Vec<f64> = (0..n_types).map(|t| config.concentration(t)).collect();
        let b: Vec<f64> = config
            .species
            .iter()
            .map(|s| neutron_scattering_length(s))
            .collect();
        let b_avg: f64 = (0..n_types).map(|a| conc[a] * b[a]).sum();
        let b_avg_sq = b_avg * b_avg;

        let mut weights = vec![0.0f64; n_pairs * nq];
        if b_avg_sq < 1e-30 {
            return weights;
        }
        for a in 0..n_types {
            for b_idx in a..n_types {
                let pair_idx = pair_index(a, b_idx, n_types);
                let dab = if a == b_idx { 1.0 } else { 2.0 };
                let w = dab * conc[a] * conc[b_idx] * b[a] * b[b_idx] / b_avg_sq;
                for k in 0..nq {
                    weights[pair_idx * nq + k] = w;
                }
            }
        }
        weights
    }

    /// Decompose total residual ΔS(Q) into partial residuals ΔS_ab(Q).
    ///
    /// Uses proportional decomposition:
    ///   ΔS_ab(Q) = w_ab(Q) * ΔS(Q) / Σ_cd w_cd(Q)²
    ///
    /// Returns flat array [n_pairs * nq].
    pub fn compute_residual_partials(
        _partial_sq_sim: &[f64],
        total_sq_sim: &[f64],
        total_sq_exp: &[f64],
        xray_weights: &[f64],
        n_pairs: usize,
        nq: usize,
    ) -> Vec<f64> {
        let mut delta_partials = vec![0.0f64; n_pairs * nq];

        for k in 0..nq {
            // Sign: sim - exp (matches EPSR25 convention where the fitted
            // residual is added to the EP to drive the simulation toward experiment)
            let delta_s = total_sq_sim[k] - total_sq_exp[k];

            // Sum of squared weights at this Q
            let mut w2_sum = 0.0f64;
            for p in 0..n_pairs {
                let w = xray_weights[p * nq + k];
                w2_sum += w * w;
            }
            if w2_sum < 1e-30 {
                continue;
            }

            for p in 0..n_pairs {
                let w = xray_weights[p * nq + k];
                // Proportional decomposition
                delta_partials[p * nq + k] = w * delta_s / w2_sum;
            }
        }

        // Also include the per-partial residual: sim partial vs what it "should" be
        // Actually, the simpler EPSR approach uses ΔS_ab = S_ab_target - S_ab_sim
        // where S_ab_target is obtained from the total decomposition.
        // The proportional decomposition above gives the correction per partial.
        // We use: ΔS_ab(Q) ≈ w_ab(Q) * [S_exp(Q) - S_sim(Q)] / Σ w_cd²
        // This is the standard Soper approach.

        delta_partials
    }

    /// Discrete sine transform: ΔS_ab(Q) → Δg_ab(r) on the EP grid.
    ///
    /// Δg_ab(r) = (2/π) * dr_q * Σ_k Q_k * ΔS_ab(Q_k) * sin(Q_k * r) / r
    pub fn fourier_to_real_space(
        delta_sq: &[f64],
        q_grid: &[f64],
        r_grid: &[f64],
        rho0: f64,
    ) -> Vec<f64> {
        let nq = q_grid.len();
        let nr = r_grid.len();
        let dq = if nq > 1 { q_grid[1] - q_grid[0] } else { 1.0 };
        let prefactor = dq / (2.0 * PI * PI * rho0);

        let mut delta_gr = vec![0.0f64; nr];
        for i in 0..nr {
            let r = r_grid[i];
            if r < 1e-10 {
                continue;
            }
            let mut sum = 0.0;
            for k in 0..nq {
                let q = q_grid[k];
                // ΔS_ab here is (S_ab - 1) form delta, so the transform gives Δg
                sum += q * delta_sq[k] * (q * r).sin();
            }
            delta_gr[i] = prefactor * sum / r;
        }
        delta_gr
    }

    /// Apply the EP update: EP_ab(r) += feedback * kT * Δg_ab(r),
    /// Gaussian smooth, zero below min_r.
    ///
    /// Returns `(max_delta, max_ep)`: the max |ΔEP| across all pairs and the
    /// max |EP| after accumulation, for convergence checking.
    pub fn update(
        &mut self,
        delta_partials_sq: &[f64],
        q_grid: &[f64],
        rho0: f64,
        feedback: f64,
        kt: f64,
        smooth_sigma: f64,
        min_r: f64,
    ) -> (f64, f64) {
        let nq = q_grid.len();
        let r_grid: Vec<f64> = (0..self.n_bins).map(|i| i as f64 * self.dr).collect();
        let mut max_delta = 0.0f64;

        for p in 0..self.n_pairs {
            let sq_slice = &delta_partials_sq[p * nq..(p + 1) * nq];

            // Sine transform to real space
            let delta_gr = Self::fourier_to_real_space(sq_slice, q_grid, &r_grid, rho0);

            // Compute raw delta EP
            let mut delta_ep = vec![0.0f64; self.n_bins];
            for i in 0..self.n_bins {
                delta_ep[i] = feedback * kt * delta_gr[i];
            }

            // Zero below min_r
            for i in 0..self.n_bins {
                if r_grid[i] < min_r {
                    delta_ep[i] = 0.0;
                }
            }

            // Gaussian smooth
            if smooth_sigma > 0.0 {
                delta_ep = gaussian_smooth(&delta_ep, smooth_sigma, self.dr);
            }

            // Zero below min_r again after smoothing
            for i in 0..self.n_bins {
                if r_grid[i] < min_r {
                    delta_ep[i] = 0.0;
                }
            }

            // Track max absolute delta
            for &d in &delta_ep {
                if d.abs() > max_delta {
                    max_delta = d.abs();
                }
            }

            // Accumulate into EP
            for i in 0..self.n_bins {
                self.ep_tables[p][i] += delta_ep[i];
            }
        }

        // Compute max |EP| across all pairs after accumulation
        let max_ep = self
            .ep_tables
            .iter()
            .flat_map(|t| t.iter())
            .fold(0.0f64, |m, &v| m.max(v.abs()));

        (max_delta, max_ep)
    }

    /// Build a combined PotentialSet = reference + EP tables.
    ///
    /// Clones the reference potential set and adds the current EP tables.
    /// The EP is interpolated onto the reference potential's grid (dr=0.001)
    /// to ensure correct element-wise addition.
    /// If no reference potential exists, creates a new set from EP alone.
    pub fn build_combined_potential(
        &self,
        reference: Option<&PotentialSet>,
        species: &[String],
        rdf_cutoff: f64,
    ) -> PotentialSet {
        let n_types = species.len();

        let mut combined = if let Some(refp) = reference {
            refp.clone()
        } else {
            PotentialSet {
                potentials: Vec::new(),
                weight: 0.001,
                cutoff: self.cutoff.min(rdf_cutoff),
                potential_index: vec![usize::MAX; n_types * n_types],
                n_types,
            }
        };

        // Determine target grid spacing: use reference potential's dr if available,
        // otherwise use the EP grid spacing.
        let target_dr = if let Some(refp) = reference {
            if let Some(pot) = refp.potentials.first() {
                pot.dr
            } else {
                self.dr
            }
        } else {
            self.dr
        };

        // Add each EP table, interpolated onto the target grid
        for p in 0..self.n_pairs {
            let (type_a, type_b) = self.pair_types[p];
            let target_cutoff = combined.cutoff.min(self.cutoff);
            let target_nbins = (target_cutoff / target_dr).ceil() as usize + 1;

            // Interpolate EP from its native grid onto the target grid
            let mut interp_table = vec![0.0f64; target_nbins];
            for i in 0..target_nbins {
                let r = i as f64 * target_dr;
                interp_table[i] = self.evaluate_ep(p, r);
            }

            let ep_pot = PairPotential::from_vec(
                format!("EP_{}", self.pair_labels[p]),
                type_a,
                type_b,
                interp_table,
                target_cutoff,
                target_dr,
            );
            combined.add_potential(ep_pot);
        }

        combined
    }

    /// Evaluate EP for pair `p` at distance `r` by linear interpolation on the EP grid.
    fn evaluate_ep(&self, p: usize, r: f64) -> f64 {
        if r >= self.cutoff {
            return 0.0;
        }
        let table = &self.ep_tables[p];
        let x = r / self.dr;
        let i = x as usize;
        if i + 1 >= self.n_bins {
            return 0.0;
        }
        let t = x - i as f64;
        table[i] + t * (table[i + 1] - table[i])
    }

    /// Write EP tables to files: `epsr_ep_{pair}.dat`.
    pub fn write_potentials(&self, output_dir: &Path) -> std::io::Result<()> {
        use std::io::Write;

        for p in 0..self.n_pairs {
            let filename = format!("epsr_ep_{}.dat", self.pair_labels[p]);
            let path = output_dir.join(filename);
            let mut file = std::fs::File::create(&path)?;
            writeln!(
                file,
                "# r(A)  EP(eV)  -- EPSR empirical potential for {}",
                self.pair_labels[p]
            )?;
            for i in 0..self.n_bins {
                let r = i as f64 * self.dr;
                writeln!(file, "{:.6} {:.8}", r, self.ep_tables[p][i])?;
            }
        }
        Ok(())
    }
}

/// 1D Gaussian convolution with kernel width sigma (in Å) on a uniform grid with spacing dr.
pub fn gaussian_smooth(data: &[f64], sigma: f64, dr: f64) -> Vec<f64> {
    let n = data.len();
    if sigma < dr * 0.1 || n == 0 {
        return data.to_vec();
    }

    // Kernel half-width: 4 sigma
    let hw = (4.0 * sigma / dr).ceil() as usize;
    let kernel_len = 2 * hw + 1;
    let mut kernel = vec![0.0f64; kernel_len];
    let inv_2s2 = 1.0 / (2.0 * sigma * sigma);
    let mut ksum = 0.0;
    for j in 0..kernel_len {
        let x = (j as f64 - hw as f64) * dr;
        kernel[j] = (-x * x * inv_2s2).exp();
        ksum += kernel[j];
    }
    // Normalise
    for k in &mut kernel {
        *k /= ksum;
    }

    let mut out = vec![0.0f64; n];
    for i in 0..n {
        let mut val = 0.0;
        for j in 0..kernel_len {
            let idx = i as isize + j as isize - hw as isize;
            if idx >= 0 && (idx as usize) < n {
                val += data[idx as usize] * kernel[j];
            }
        }
        out[i] = val;
    }
    out
}

/// Interpolate experimental S(Q) onto the simulation Q grid.
///
/// Uses linear interpolation. Points outside the experimental range are set to 1.0.
pub fn interpolate_exp_to_grid(q_exp: &[f64], sq_exp: &[f64], q_grid: &[f64]) -> Vec<f64> {
    q_grid
        .iter()
        .map(|&q| {
            if q <= q_exp[0] {
                return sq_exp[0];
            }
            if q >= q_exp[q_exp.len() - 1] {
                return sq_exp[sq_exp.len() - 1];
            }
            // Binary search
            let mut lo = 0;
            let mut hi = q_exp.len() - 1;
            while hi - lo > 1 {
                let mid = (lo + hi) / 2;
                if q_exp[mid] <= q {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let t = (q - q_exp[lo]) / (q_exp[hi] - q_exp[lo]);
            sq_exp[lo] + t * (sq_exp[hi] - sq_exp[lo])
        })
        .collect()
}

/// Map (type_a, type_b) with a <= b to a flat pair index.
#[inline]
fn pair_index(a: usize, b: usize, n_types: usize) -> usize {
    let (a, b) = if a <= b { (a, b) } else { (b, a) };
    a * n_types - a * (a + 1) / 2 + b
}

/// Linear interpolation on non-uniform data. Returns 0.0 outside the data range.
fn interp_linear(x: &[f64], y: &[f64], xi: f64) -> f64 {
    if x.is_empty() {
        return 0.0;
    }
    if xi <= x[0] {
        return y[0];
    }
    if xi >= x[x.len() - 1] {
        return y[y.len() - 1];
    }
    let mut lo = 0;
    let mut hi = x.len() - 1;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if x[mid] <= xi {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let t = (xi - x[lo]) / (x[hi] - x[lo]);
    y[lo] + t * (y[hi] - y[lo])
}
