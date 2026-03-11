use std::collections::HashMap;

use crate::atoms::Configuration;

/// Partial pair distribution functions g_ab(r) computed from atomic positions.
pub struct PartialRdfs {
    /// r grid (bin centres).
    pub r: Vec<f64>,
    /// g_ab(r) keyed by pair index (from Configuration::pair_index).
    pub partials: HashMap<usize, Vec<f64>>,
    /// Number of bins.
    pub nbins: usize,
    /// Bin width.
    pub dr: f64,
    /// Cutoff radius.
    pub cutoff: f64,
}

/// Compute all partial g_ab(r) from a Configuration.
///
/// Uses histogram binning with minimum-image convention.
/// Normalisation: g_ab(r) = histogram / (4*pi*r^2*dr * N_a * rho_b)
/// where rho_b = N_b / V.
pub fn compute_partial_rdfs(config: &Configuration, nbins: usize, cutoff: f64) -> PartialRdfs {
    let dr = cutoff / nbins as f64;
    let histograms = compute_histograms(config, nbins, cutoff);
    let partials = normalise_histograms(&histograms, config, nbins, dr);
    let r_grid: Vec<f64> = (0..nbins).map(|i| (i as f64 + 0.5) * dr).collect();

    PartialRdfs {
        r: r_grid,
        partials,
        nbins,
        dr,
        cutoff,
    }
}

/// Compute histogram contributions from a single atom (for incremental RMC updates).
/// Returns histogram increments as a flat Vec indexed by pair_index * nbins + bin.
/// This avoids HashMap overhead in the hot path.
pub fn atom_histogram_flat(
    config: &Configuration,
    atom_idx: usize,
    pos: &[f64; 3],
    nbins: usize,
    cutoff: f64,
    dr: f64,
    n_pairs: usize,
) -> Vec<f64> {
    let ti = config.atoms[atom_idx].type_id;
    let cutoff2 = cutoff * cutoff;
    let inv_dr = 1.0 / dr;
    let box_lengths = &config.box_lengths;

    let n = config.atoms.len();
    let mut hist = vec![0.0f64; n_pairs * nbins];

    for j in 0..n {
        if j == atom_idx {
            continue;
        }
        let tj = config.atoms[j].type_id;
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
                let pair_idx = if ti <= tj {
                    let n_sp = config.species.len();
                    ti * n_sp - ti * (ti + 1) / 2 + tj
                } else {
                    let n_sp = config.species.len();
                    tj * n_sp - tj * (tj + 1) / 2 + ti
                };
                hist[pair_idx * nbins + bin] += 1.0;
            }
        }
    }

    hist
}

/// Rebuild full histograms from scratch (used for initial setup).
/// Returns raw histograms (not normalised) keyed by pair_index.
pub fn compute_histograms(
    config: &Configuration,
    nbins: usize,
    cutoff: f64,
) -> HashMap<usize, Vec<f64>> {
    let dr = cutoff / nbins as f64;
    let inv_dr = 1.0 / dr;
    let n_pairs = config.num_type_pairs();
    let n = config.atoms.len();
    let cutoff2 = cutoff * cutoff;
    let box_lengths = config.box_lengths;

    let mut flat_hist = vec![0.0f64; n_pairs * nbins];
    for i in 0..n {
        let ti = config.atoms[i].type_id;
        let pi = &config.atoms[i].position;
        for j in (i + 1)..n {
            let tj = config.atoms[j].type_id;
            let pj = &config.atoms[j].position;

            let mut r2 = 0.0f64;
            for d in 0..3 {
                let mut delta = pj[d] - pi[d];
                let l = box_lengths[d];
                delta -= l * (delta / l).round();
                r2 += delta * delta;
            }

            if r2 < cutoff2 {
                let r = r2.sqrt();
                let bin = (r * inv_dr) as usize;
                if bin < nbins {
                    let pair_idx = config.pair_index(ti, tj);
                    flat_hist[pair_idx * nbins + bin] += 1.0;
                }
            }
        }
    }

    // Convert flat histogram to HashMap
    let mut histograms: HashMap<usize, Vec<f64>> = HashMap::new();
    for p in 0..n_pairs {
        let start = p * nbins;
        histograms.insert(p, flat_hist[start..start + nbins].to_vec());
    }

    histograms
}

/// Normalise raw histograms to g_ab(r).
pub fn normalise_histograms(
    histograms: &HashMap<usize, Vec<f64>>,
    config: &Configuration,
    nbins: usize,
    dr: f64,
) -> HashMap<usize, Vec<f64>> {
    let volume = config.volume();
    let n_types = config.species.len();
    let r_grid: Vec<f64> = (0..nbins).map(|i| (i as f64 + 0.5) * dr).collect();

    // Collect pair metadata
    let mut pair_info: Vec<(usize, bool, f64, f64)> = Vec::new();
    for a in 0..n_types {
        for b in a..n_types {
            let pair_idx = config.pair_index(a, b);
            let n_a = config.count_type(a) as f64;
            let n_b = config.count_type(b) as f64;
            let rho_b = n_b / volume;
            pair_info.push((pair_idx, a == b, n_a, rho_b));
        }
    }

    let results: Vec<(usize, Vec<f64>)> = pair_info
        .iter()
        .map(|&(pair_idx, is_like, n_a, rho_b)| {
            let hist = &histograms[&pair_idx];
            let gr: Vec<f64> = (0..nbins)
                .map(|bin| {
                    let r = r_grid[bin];
                    let shell_vol = 4.0 * std::f64::consts::PI * r * r * dr;
                    let norm = n_a * rho_b * shell_vol;
                    if norm > 0.0 {
                        if is_like {
                            2.0 * hist[bin] / norm
                        } else {
                            hist[bin] / norm
                        }
                    } else {
                        0.0
                    }
                })
                .collect();
            (pair_idx, gr)
        })
        .collect();

    results.into_iter().collect()
}
