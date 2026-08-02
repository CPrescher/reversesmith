use std::collections::HashMap;
use std::f64::consts::PI;
use std::path::Path;

use rsmith::atoms::{Atom, Configuration};
use rsmith::{io, neutron, rdf, sq, xray};

fn fixture_path() -> &'static Path {
    Path::new("tests/data/scattering/small_sio2.data")
}

fn type_map() -> HashMap<u32, String> {
    HashMap::from([(1, "Si".to_string()), (2, "O".to_string())])
}

#[test]
fn committed_lammps_scattering_fixture_loads() {
    let config = io::read_lammps_data(fixture_path(), &type_map()).unwrap();

    assert_eq!(config.atoms.len(), 6);
    assert_eq!(config.species, ["Si", "O"]);
    assert_eq!(config.composition["Si"], 2);
    assert_eq!(config.composition["O"], 4);
    assert_eq!(config.box_lengths, [10.0, 10.0, 10.0]);
}

#[test]
fn partial_rdf_histograms_and_normalisation_match_the_fixture_geometry() {
    let config = io::read_lammps_data(fixture_path(), &type_map()).unwrap();
    let nbins = 16;
    let cutoff = 4.0;
    let dr = cutoff / nbins as f64;
    let histograms = rdf::compute_histograms(&config, nbins, cutoff);

    let si_o = config.pair_index(0, 1);
    let o_o = config.pair_index(1, 1);
    let si_o_bin = (1.5 / dr) as usize;
    let o_o_distance = (2.0_f64 * 1.5 * 1.5).sqrt();
    let o_o_bin = (o_o_distance / dr) as usize;

    assert_eq!(histograms[&si_o].iter().sum::<f64>(), 4.0);
    assert_eq!(histograms[&si_o][si_o_bin], 4.0);
    assert_eq!(histograms[&o_o].iter().sum::<f64>(), 2.0);
    assert_eq!(histograms[&o_o][o_o_bin], 2.0);

    let partials = rdf::normalise_histograms(&histograms, &config, nbins, dr);
    let volume = config.volume();

    let r_si_o = (si_o_bin as f64 + 0.5) * dr;
    let shell_si_o = 4.0 * PI * r_si_o * r_si_o * dr;
    let expected_si_o = 4.0 / (2.0 * (4.0 / volume) * shell_si_o);
    assert!((partials[&si_o][si_o_bin] - expected_si_o).abs() < 1.0e-12);

    let r_o_o = (o_o_bin as f64 + 0.5) * dr;
    let shell_o_o = 4.0 * PI * r_o_o * r_o_o * dr;
    let expected_o_o = 2.0 * 2.0 / (4.0 * (4.0 / volume) * shell_o_o);
    assert!((partials[&o_o][o_o_bin] - expected_o_o).abs() < 1.0e-12);
}

#[test]
fn direct_partial_sq_matches_the_documented_quadrature() {
    let r = vec![0.5, 1.0, 1.5, 2.0];
    let gr = vec![0.2, 1.4, 0.8, 1.1];
    let q = vec![0.4, 1.3, 2.7];
    let rho0 = 0.067;
    let dr = 0.5;
    let rmax = 2.0;

    for lorch in [false, true] {
        let actual = sq::compute_partial_sq(&r, &gr, rho0, &q, lorch);
        for (k, &qk) in q.iter().enumerate() {
            let sum: f64 = r
                .iter()
                .zip(&gr)
                .map(|(&ri, &gi)| {
                    let window = if lorch {
                        let arg = PI * ri / rmax;
                        arg.sin() / arg
                    } else {
                        1.0
                    };
                    ri * (gi - 1.0) * window * (qk * ri).sin()
                })
                .sum();
            let expected = 1.0 + 4.0 * PI * rho0 * dr * sum / qk;
            assert!((actual[k] - expected).abs() < 1.0e-13);
        }
    }
}

fn two_species_configuration() -> Configuration {
    Configuration {
        atoms: vec![
            Atom {
                position: [1.0, 1.0, 1.0],
                species: "Si".to_string(),
                type_id: 0,
            },
            Atom {
                position: [3.0, 3.0, 3.0],
                species: "O".to_string(),
                type_id: 1,
            },
        ],
        box_lengths: [8.0, 8.0, 8.0],
        species: vec!["Si".to_string(), "O".to_string()],
        composition: HashMap::from([("Si".to_string(), 1), ("O".to_string(), 1)]),
    }
}

#[test]
fn xray_and_neutron_weights_preserve_a_common_partial_curve() {
    let config = two_species_configuration();
    let q = vec![0.5, 1.0, 2.0, 4.0];
    let common = vec![0.8, 1.2, 0.95, 1.05];
    let partials: HashMap<usize, Vec<f64>> = (0..config.num_type_pairs())
        .map(|pair| (pair, common.clone()))
        .collect();

    let xray_total = xray::compute_xray_sq(&config, &partials, &q);
    let neutron_total = neutron::compute_sq(&config, &partials, &q);

    for ((&expected, &xray), &neutron) in common.iter().zip(&xray_total).zip(&neutron_total) {
        assert!((xray - expected).abs() < 1.0e-12);
        assert!((neutron - expected).abs() < 1.0e-12);
    }
}

#[test]
fn minimum_image_distance_crosses_all_three_boundaries() {
    let config = Configuration {
        atoms: vec![
            Atom {
                position: [0.5, 0.5, 0.5],
                species: "Si".to_string(),
                type_id: 0,
            },
            Atom {
                position: [9.5, 9.5, 9.5],
                species: "Si".to_string(),
                type_id: 0,
            },
        ],
        box_lengths: [10.0, 10.0, 10.0],
        species: vec!["Si".to_string()],
        composition: HashMap::from([("Si".to_string(), 2)]),
    };

    assert!((config.distance(0, 1) - 3.0_f64.sqrt()).abs() < 1.0e-12);
}
