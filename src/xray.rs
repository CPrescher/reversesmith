use std::collections::HashMap;

use rayon::prelude::*;

use crate::atoms::Configuration;

/// Cromer-Mann coefficients for X-ray atomic form factors.
/// f(s) = sum_i a_i * exp(-b_i * s^2) + c, where s = Q/(4*pi).
///
/// Coefficients from International Tables for Crystallography, Vol C (2004).
struct CromerMann {
    a: [f64; 4],
    b: [f64; 4],
    c: f64,
}

fn cromer_mann_params(element: &str) -> CromerMann {
    match element {
        "Ca" => CromerMann {
            a: [8.6266, 7.3873, 1.5899, 1.0211],
            b: [10.4421, 0.6599, 85.7484, 178.437],
            c: 1.3751,
        },
        "Si" => CromerMann {
            a: [6.2915, 3.0353, 1.9891, 1.541],
            b: [2.4386, 32.3337, 0.6785, 81.6937],
            c: 1.1407,
        },
        "O" => CromerMann {
            a: [3.0485, 2.2868, 1.5463, 0.867],
            b: [13.2771, 5.7011, 0.3239, 32.9089],
            c: 0.2508,
        },
        "Na" => CromerMann {
            a: [4.7626, 3.1736, 1.2674, 1.1128],
            b: [3.285, 8.8422, 0.3136, 129.424],
            c: 0.676,
        },
        "Al" => CromerMann {
            a: [6.4202, 1.9002, 1.5936, 1.9646],
            b: [3.0387, 0.7426, 31.5472, 85.0886],
            c: 1.1151,
        },
        "Mg" => CromerMann {
            a: [5.4204, 2.1735, 1.2269, 2.3073],
            b: [2.8275, 79.2611, 0.3808, 7.1937],
            c: 0.8584,
        },
        _ => panic!("No Cromer-Mann parameters for element '{}'", element),
    }
}

/// Compute X-ray atomic form factor f(Q) for an element.
pub fn form_factor(element: &str, q: &[f64]) -> Vec<f64> {
    let cm = cromer_mann_params(element);
    q.iter()
        .map(|&qi| {
            let s = qi / (4.0 * std::f64::consts::PI);
            let s2 = s * s;
            let mut f = cm.c;
            for i in 0..4 {
                f += cm.a[i] * (-cm.b[i] * s2).exp();
            }
            f
        })
        .collect()
}

/// Compute total X-ray structure factor S_X(Q) from partial S_ab(Q).
///
/// S_X(Q) = sum_{a<=b} (2-delta_ab) * c_a*c_b*f_a*f_b / <f>^2 * S_ab(Q)
///
/// Uses Faber-Ziman convention. Parallel over Q points.
pub fn compute_xray_sq(
    config: &Configuration,
    partial_sq: &HashMap<usize, Vec<f64>>,
    q: &[f64],
) -> Vec<f64> {
    let n_types = config.species.len();

    // Precompute form factors for each species
    let form_factors: Vec<Vec<f64>> = config
        .species
        .iter()
        .map(|s| form_factor(s, q))
        .collect();

    // Precompute concentrations
    let conc: Vec<f64> = (0..n_types).map(|t| config.concentration(t)).collect();

    // Collect pair info for inner loop
    let mut pairs: Vec<(usize, usize, usize, f64)> = Vec::new();
    for a in 0..n_types {
        for b in a..n_types {
            let pair_idx = config.pair_index(a, b);
            let delta_ab = if a == b { 1.0 } else { 2.0 };
            pairs.push((a, b, pair_idx, delta_ab));
        }
    }

    (0..q.len())
        .into_par_iter()
        .map(|iq| {
            let f_avg: f64 = (0..n_types)
                .map(|a| conc[a] * form_factors[a][iq])
                .sum();
            let f_avg_sq = f_avg * f_avg;

            if f_avg_sq < 1e-30 {
                return 1.0;
            }

            let mut total = 0.0;
            for &(a, b, pair_idx, delta_ab) in &pairs {
                let weight =
                    delta_ab * conc[a] * conc[b] * form_factors[a][iq] * form_factors[b][iq]
                        / f_avg_sq;
                if let Some(sq) = partial_sq.get(&pair_idx) {
                    total += weight * sq[iq];
                }
            }
            total
        })
        .collect()
}

/// Neutron coherent scattering lengths (in fm) for common elements.
/// From NIST neutron scattering length tables.
pub fn neutron_scattering_length(element: &str) -> f64 {
    match element {
        "Ca" => 4.70,
        "Si" => 4.1491,
        "O" => 5.803,
        "Na" => 3.63,
        "Al" => 3.449,
        "Mg" => 5.375,
        _ => panic!("No neutron scattering length for element '{}'", element),
    }
}

/// Compute total neutron structure factor from partial S_ab(Q).
/// Uses Q-independent coherent scattering lengths (Faber-Ziman). Parallel over Q points.
pub fn compute_neutron_sq(
    config: &Configuration,
    partial_sq: &HashMap<usize, Vec<f64>>,
    q: &[f64],
) -> Vec<f64> {
    let n_types = config.species.len();
    let conc: Vec<f64> = (0..n_types).map(|t| config.concentration(t)).collect();
    let b: Vec<f64> = config
        .species
        .iter()
        .map(|s| neutron_scattering_length(s))
        .collect();

    let b_avg: f64 = (0..n_types).map(|a| conc[a] * b[a]).sum();
    let b_avg_sq = b_avg * b_avg;

    // Collect pair info
    let mut pairs: Vec<(usize, usize, usize, f64)> = Vec::new();
    for a in 0..n_types {
        for b_idx in a..n_types {
            let pair_idx = config.pair_index(a, b_idx);
            let delta_ab = if a == b_idx { 1.0 } else { 2.0 };
            let weight = delta_ab * conc[a] * conc[b_idx] * b[a] * b[b_idx] / b_avg_sq;
            pairs.push((pair_idx, 0, 0, weight));
        }
    }

    (0..q.len())
        .into_par_iter()
        .map(|iq| {
            let mut total = 0.0;
            for &(pair_idx, _, _, weight) in &pairs {
                if let Some(sq) = partial_sq.get(&pair_idx) {
                    total += weight * sq[iq];
                }
            }
            total
        })
        .collect()
}
