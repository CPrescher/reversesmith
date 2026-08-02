use std::collections::HashMap;

use crate::atoms::Configuration;

/// Neutron coherent scattering lengths (in fm).
/// From NIST Center for Neutron Research, neutron scattering length tables.
pub fn scattering_length(element: &str) -> f64 {
    match element {
        "H" => -3.7390,
        "He" => 3.26,
        "Li" => -1.90,
        "Be" => 7.79,
        "B" => 5.30, // 11B; natural B has large absorption
        "C" => 6.6460,
        "N" => 9.36,
        "O" => 5.803,
        "F" => 5.654,
        "Ne" => 4.566,
        "Na" => 3.63,
        "Mg" => 5.375,
        "Al" => 3.449,
        "Si" => 4.1491,
        "P" => 5.13,
        "S" => 2.847,
        "Cl" => 9.577,
        "Ar" => 1.909,
        "K" => 3.67,
        "Ca" => 4.70,
        "Sc" => 12.29,
        "Ti" => -3.438,
        "V" => -0.3824,
        "Cr" => 3.635,
        "Mn" => -3.73,
        "Fe" => 9.45,
        "Co" => 2.49,
        "Ni" => 10.3,
        "Cu" => 7.718,
        "Zn" => 5.680,
        "Ga" => 7.288,
        "Ge" => 8.185,
        "As" => 6.58,
        "Se" => 7.970,
        "Br" => 6.795,
        "Kr" => 7.81,
        "Rb" => 7.09,
        "Sr" => 7.02,
        "Y" => 7.75,
        "Zr" => 7.16,
        "Nb" => 7.054,
        "Mo" => 6.715,
        "Ru" => 7.03,
        "Rh" => 5.88,
        "Pd" => 5.91,
        "Ag" => 5.922,
        "Cd" => 4.87, // natural Cd; large absorption from 113Cd
        "In" => 4.065,
        "Sn" => 6.225,
        "Sb" => 5.57,
        "Te" => 5.80,
        "I" => 5.28,
        "Xe" => 4.92,
        "Cs" => 5.42,
        "Ba" => 5.07,
        "La" => 8.24,
        "Ce" => 4.84,
        "Pr" => 4.58,
        "Nd" => 7.69,
        "Sm" => 0.80, // large absorption
        "Eu" => 7.22,
        "Gd" => 6.50, // large absorption from 157Gd
        "Tb" => 7.38,
        "Dy" => 16.9,
        "Ho" => 8.01,
        "Er" => 7.79,
        "Tm" => 7.07,
        "Yb" => 12.43,
        "Lu" => 7.21,
        "Hf" => 7.77,
        "Ta" => 6.91,
        "W" => 4.86,
        "Re" => 9.2,
        "Os" => 10.7,
        "Ir" => 10.6,
        "Pt" => 9.60,
        "Au" => 7.63,
        "Hg" => 12.692,
        "Tl" => 8.776,
        "Pb" => 9.405,
        "Bi" => 8.532,
        "Th" => 10.31,
        "U" => 8.417,
        _ => panic!(
            "No neutron scattering length for element '{}'. See NIST tables.",
            element
        ),
    }
}

/// Compute total neutron structure factor from partial S_ab(Q).
/// Uses Q-independent coherent scattering lengths (Faber-Ziman).
pub fn compute_sq(
    config: &Configuration,
    partial_sq: &HashMap<usize, Vec<f64>>,
    q: &[f64],
) -> Vec<f64> {
    let b: Vec<f64> = config
        .species
        .iter()
        .map(|s| scattering_length(s))
        .collect();
    compute_sq_with_lengths(config, partial_sq, q, &b)
}

/// Compute a total neutron structure factor with caller-supplied coherent
/// scattering lengths. This supports isotope-enriched and null-scattering
/// contrasts without pretending that all samples have natural isotopic abundance.
pub fn compute_sq_with_lengths(
    config: &Configuration,
    partial_sq: &HashMap<usize, Vec<f64>>,
    q: &[f64],
    scattering_lengths: &[f64],
) -> Vec<f64> {
    let n_types = config.species.len();
    assert_eq!(
        scattering_lengths.len(),
        n_types,
        "one neutron scattering length is required per species"
    );
    let conc: Vec<f64> = (0..n_types).map(|t| config.concentration(t)).collect();
    let b = scattering_lengths;

    let b_avg: f64 = (0..n_types).map(|a| conc[a] * b[a]).sum();
    let b_avg_sq = b_avg * b_avg;
    assert!(
        b_avg_sq > 1e-30,
        "the concentration-weighted mean neutron scattering length is zero"
    );

    let mut pairs: Vec<(usize, f64)> = Vec::new();
    for a in 0..n_types {
        for b_idx in a..n_types {
            let pair_idx = config.pair_index(a, b_idx);
            let delta_ab = if a == b_idx { 1.0 } else { 2.0 };
            let weight = delta_ab * conc[a] * conc[b_idx] * b[a] * b[b_idx] / b_avg_sq;
            pairs.push((pair_idx, weight));
        }
    }

    (0..q.len())
        .map(|iq| {
            let mut total = 0.0;
            for &(pair_idx, weight) in &pairs {
                if let Some(sq) = partial_sq.get(&pair_idx) {
                    total += weight * sq[iq];
                }
            }
            total
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atoms::{Atom, Configuration};

    #[test]
    fn isotope_contrast_uses_supplied_scattering_lengths() {
        let config = Configuration {
            atoms: vec![
                Atom {
                    position: [0.0; 3],
                    species: "Ge".to_string(),
                    type_id: 0,
                },
                Atom {
                    position: [1.0; 3],
                    species: "O".to_string(),
                    type_id: 1,
                },
            ],
            box_lengths: [10.0; 3],
            species: vec!["Ge".to_string(), "O".to_string()],
            composition: HashMap::from([("Ge".to_string(), 1), ("O".to_string(), 1)]),
        };
        let partials = HashMap::from([
            (config.pair_index(0, 0), vec![2.0]),
            (config.pair_index(0, 1), vec![3.0]),
            (config.pair_index(1, 1), vec![4.0]),
        ]);

        let result = compute_sq_with_lengths(&config, &partials, &[1.0], &[10.0, 5.0]);
        // c_Ge = c_O = 0.5 and <b> = 7.5:
        // weights are 25/56.25, 25/56.25, and 6.25/56.25.
        let expected = (25.0 * 2.0 + 25.0 * 3.0 + 6.25 * 4.0) / 56.25;
        assert!((result[0] - expected).abs() < 1e-12);
    }
}
