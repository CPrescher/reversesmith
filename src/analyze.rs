use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use crate::atoms::Configuration;
use crate::{log_println, log_print};

/// A pair of species with a distance cutoff for analysis.
#[derive(Debug, Clone)]
pub struct AnalysisPair {
    pub species_a: String,
    pub species_b: String,
    pub cutoff: f64,
}

/// Coordination number distribution for one pair.
#[derive(Debug)]
pub struct CnDistribution {
    pub pair_label: String,
    pub cutoff: f64,
    pub histogram: Vec<usize>,
    pub mean: f64,
    pub std_dev: f64,
    pub min: usize,
    pub max: usize,
}

/// Bond angle distribution for one triplet X-M-Y.
#[derive(Debug)]
pub struct AngleDistribution {
    pub triplet_label: String,
    pub histogram: Vec<f64>,
    pub bin_centres: Vec<f64>,
    pub n_angles: usize,
    pub peak_angle: f64,
    pub mean_angle: f64,
}

/// Parse pair cutoff map into AnalysisPair list.
pub fn build_analysis_pairs(cutoffs: &HashMap<String, f64>) -> Vec<AnalysisPair> {
    cutoffs
        .iter()
        .map(|(pair, &cutoff)| {
            let parts: Vec<&str> = pair.split('-').collect();
            AnalysisPair {
                species_a: parts[0].to_string(),
                species_b: parts.get(1).unwrap_or(&parts[0]).to_string(),
                cutoff,
            }
        })
        .collect()
}

/// Compute coordination number distributions for each pair.
pub fn compute_coordination_numbers(
    config: &Configuration,
    pairs: &[AnalysisPair],
) -> Vec<CnDistribution> {
    let mut results = Vec::new();

    for pair in pairs {
        let type_a = config.species.iter().position(|s| s == &pair.species_a);
        let type_b = config.species.iter().position(|s| s == &pair.species_b);

        let (type_a, type_b) = match (type_a, type_b) {
            (Some(a), Some(b)) => (a, b),
            _ => continue,
        };

        let cutoff2 = pair.cutoff * pair.cutoff;
        let mut cn_values: Vec<usize> = Vec::new();

        for (i, atom_i) in config.atoms.iter().enumerate() {
            if atom_i.type_id != type_a {
                continue;
            }
            let mut count = 0usize;
            for (j, atom_j) in config.atoms.iter().enumerate() {
                if i == j {
                    continue;
                }
                if atom_j.type_id != type_b {
                    continue;
                }
                let dr = config.distance_vec(i, j);
                let r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                if r2 < cutoff2 {
                    count += 1;
                }
            }
            cn_values.push(count);
        }

        if cn_values.is_empty() {
            continue;
        }

        let min_cn = *cn_values.iter().min().unwrap();
        let max_cn = *cn_values.iter().max().unwrap();
        let n = cn_values.len() as f64;
        let mean = cn_values.iter().sum::<usize>() as f64 / n;
        let variance = cn_values.iter().map(|&c| {
            let d = c as f64 - mean;
            d * d
        }).sum::<f64>() / n;
        let std_dev = variance.sqrt();

        let mut histogram = vec![0usize; max_cn + 1];
        for &cn in &cn_values {
            histogram[cn] += 1;
        }

        results.push(CnDistribution {
            pair_label: format!("{}-{}", pair.species_a, pair.species_b),
            cutoff: pair.cutoff,
            histogram,
            mean,
            std_dev,
            min: min_cn,
            max: max_cn,
        });
    }

    results
}

/// Build neighbor list: for each atom, find all bonded neighbors from all pairs.
/// Returns Vec<Vec<usize>> indexed by atom index.
fn build_neighbor_list(config: &Configuration, pairs: &[AnalysisPair]) -> Vec<Vec<usize>> {
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); config.atoms.len()];

    // Map species pairs to cutoffs (both directions for asymmetric pairs)
    let mut pair_cutoff2: HashMap<(usize, usize), f64> = HashMap::new();
    for pair in pairs {
        let ta = config.species.iter().position(|s| s == &pair.species_a);
        let tb = config.species.iter().position(|s| s == &pair.species_b);
        if let (Some(a), Some(b)) = (ta, tb) {
            let c2 = pair.cutoff * pair.cutoff;
            pair_cutoff2.insert((a, b), c2);
            pair_cutoff2.insert((b, a), c2);
        }
    }

    for i in 0..config.atoms.len() {
        let ti = config.atoms[i].type_id;
        for j in (i + 1)..config.atoms.len() {
            let tj = config.atoms[j].type_id;

            // Check if this pair type has a cutoff
            if let Some(&c2) = pair_cutoff2.get(&(ti, tj)) {
                let dr = config.distance_vec(i, j);
                let r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                if r2 < c2 {
                    neighbors[i].push(j);
                    neighbors[j].push(i);
                }
            }
        }
    }

    neighbors
}

/// Compute bond angle distributions for all triplets X-M-Y.
pub fn compute_bond_angles(
    config: &Configuration,
    pairs: &[AnalysisPair],
    nbins: usize,
    triplet_filter: Option<&[String]>,
) -> Vec<AngleDistribution> {
    let neighbors = build_neighbor_list(config, pairs);

    let bin_width = std::f64::consts::PI / nbins as f64;
    let bin_centres: Vec<f64> = (0..nbins)
        .map(|i| ((i as f64 + 0.5) * bin_width).to_degrees())
        .collect();

    // Accumulate angles by triplet label
    let mut triplet_histograms: HashMap<String, Vec<f64>> = HashMap::new();
    let mut triplet_counts: HashMap<String, usize> = HashMap::new();
    let mut triplet_angle_sums: HashMap<String, f64> = HashMap::new();

    for m in 0..config.atoms.len() {
        let nbrs = &neighbors[m];
        if nbrs.len() < 2 {
            continue;
        }

        let species_m = &config.species[config.atoms[m].type_id];

        for ii in 0..nbrs.len() {
            for jj in (ii + 1)..nbrs.len() {
                let x = nbrs[ii];
                let y = nbrs[jj];

                let species_x = &config.species[config.atoms[x].type_id];
                let species_y = &config.species[config.atoms[y].type_id];

                // Sorted end-species for canonical label
                let (end_a, end_b) = if species_x <= species_y {
                    (species_x.as_str(), species_y.as_str())
                } else {
                    (species_y.as_str(), species_x.as_str())
                };
                let label = format!("{}-{}-{}", end_a, species_m, end_b);

                // Apply triplet filter if given
                if let Some(filter) = triplet_filter {
                    if !filter.iter().any(|f| f == &label) {
                        continue;
                    }
                }

                // Compute angle X-M-Y
                let dr_mx = config.distance_vec(m, x);
                let dr_my = config.distance_vec(m, y);

                let dot = dr_mx[0] * dr_my[0] + dr_mx[1] * dr_my[1] + dr_mx[2] * dr_my[2];
                let r_mx = (dr_mx[0] * dr_mx[0] + dr_mx[1] * dr_mx[1] + dr_mx[2] * dr_mx[2]).sqrt();
                let r_my = (dr_my[0] * dr_my[0] + dr_my[1] * dr_my[1] + dr_my[2] * dr_my[2]).sqrt();

                if r_mx < 1e-10 || r_my < 1e-10 {
                    continue;
                }

                let cos_theta = (dot / (r_mx * r_my)).clamp(-1.0, 1.0);
                let theta = cos_theta.acos(); // radians

                let bin = (theta / bin_width).min((nbins - 1) as f64) as usize;

                let hist = triplet_histograms
                    .entry(label.clone())
                    .or_insert_with(|| vec![0.0; nbins]);
                hist[bin] += 1.0;

                *triplet_counts.entry(label.clone()).or_insert(0) += 1;
                *triplet_angle_sums.entry(label.clone()).or_insert(0.0) += theta.to_degrees();
            }
        }
    }

    // Convert to probability density and build results
    let mut results: Vec<AngleDistribution> = Vec::new();
    let mut labels: Vec<String> = triplet_histograms.keys().cloned().collect();
    labels.sort();

    for label in labels {
        let hist = triplet_histograms.get(&label).unwrap();
        let n_angles = triplet_counts[&label];
        let angle_sum = triplet_angle_sums[&label];

        if n_angles == 0 {
            continue;
        }

        // Normalize: probability density P(theta) such that integral P(theta) dtheta = 1
        // Each bin has width bin_width_deg degrees
        let bin_width_deg = 180.0 / nbins as f64;
        let total: f64 = hist.iter().sum();
        let normalized: Vec<f64> = if total > 0.0 {
            hist.iter().map(|&h| h / (total * bin_width_deg)).collect()
        } else {
            hist.clone()
        };

        // Find peak
        let peak_bin = hist
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| i)
            .unwrap_or(0);

        results.push(AngleDistribution {
            triplet_label: label,
            histogram: normalized,
            bin_centres: bin_centres.clone(),
            n_angles,
            peak_angle: bin_centres[peak_bin],
            mean_angle: angle_sum / n_angles as f64,
        });
    }

    results
}

/// Print formatted analysis summary to stdout.
pub fn print_analysis_summary(cn_results: &[CnDistribution], angle_results: &[AngleDistribution]) {
    log_println!("\n=== Coordination Number Analysis ===\n");
    log_println!(
        "{:<10} {:>8} {:>8} {:>8} {:>8} {:>10}",
        "Pair", "Cutoff", "Mean", "StdDev", "Min", "Max"
    );
    log_println!("{}", "-".repeat(58));
    for cn in cn_results {
        log_println!(
            "{:<10} {:>8.2} {:>8.3} {:>8.3} {:>8} {:>10}",
            cn.pair_label, cn.cutoff, cn.mean, cn.std_dev, cn.min, cn.max
        );

        // Print histogram inline
        let n_atoms: usize = cn.histogram.iter().sum();
        log_print!("           CN distribution:");
        for (cn_val, &count) in cn.histogram.iter().enumerate() {
            if count > 0 {
                log_print!(
                    " {}:{} ({:.1}%)",
                    cn_val,
                    count,
                    100.0 * count as f64 / n_atoms as f64
                );
            }
        }
        log_println!();
    }

    if !angle_results.is_empty() {
        log_println!("\n=== Bond Angle Analysis ===\n");
        log_println!(
            "{:<14} {:>10} {:>12} {:>12}",
            "Triplet", "N_angles", "Peak (deg)", "Mean (deg)"
        );
        log_println!("{}", "-".repeat(52));
        for ad in angle_results {
            log_println!(
                "{:<14} {:>10} {:>12.1} {:>12.1}",
                ad.triplet_label, ad.n_angles, ad.peak_angle, ad.mean_angle
            );
        }
    }
}

/// Write CN histograms to file.
pub fn write_cn_histograms(
    path: &Path,
    results: &[CnDistribution],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(path)?;

    // Header
    write!(file, "# CN")?;
    for cn in results {
        write!(file, "  count_{0}  frac_{0}", cn.pair_label)?;
    }
    writeln!(file)?;

    // Find max CN across all distributions
    let max_cn = results.iter().map(|r| r.histogram.len()).max().unwrap_or(0);

    for cn_val in 0..max_cn {
        write!(file, "{}", cn_val)?;
        for r in results {
            let count = r.histogram.get(cn_val).copied().unwrap_or(0);
            let n_atoms: usize = r.histogram.iter().sum();
            let frac = if n_atoms > 0 {
                count as f64 / n_atoms as f64
            } else {
                0.0
            };
            write!(file, "  {}  {:.6}", count, frac)?;
        }
        writeln!(file)?;
    }

    Ok(())
}

/// Write angle histograms to multi-column file.
pub fn write_angle_histograms(
    path: &Path,
    results: &[AngleDistribution],
) -> Result<(), Box<dyn std::error::Error>> {
    if results.is_empty() {
        return Ok(());
    }

    let mut file = File::create(path)?;

    // Header
    write!(file, "# angle_deg")?;
    for ad in results {
        write!(file, "  P_{}", ad.triplet_label)?;
    }
    writeln!(file)?;

    let nbins = results[0].bin_centres.len();
    for bin in 0..nbins {
        write!(file, "{:.2}", results[0].bin_centres[bin])?;
        for ad in results {
            write!(file, "  {:.6}", ad.histogram[bin])?;
        }
        writeln!(file)?;
    }

    Ok(())
}
