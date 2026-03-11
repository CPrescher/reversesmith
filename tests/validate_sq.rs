use std::collections::HashMap;
use std::path::Path;

use reversesmith::atoms::Configuration;
use reversesmith::io;
use reversesmith::rdf;
use reversesmith::sq;
use reversesmith::xray;

#[test]
fn test_load_lammps_data() {
    let mut type_map = HashMap::new();
    type_map.insert(1, "Ca".to_string());
    type_map.insert(2, "Si".to_string());
    type_map.insert(3, "O".to_string());

    let path = Path::new("../run01/CaSiO3_glass.data");
    if !path.exists() {
        eprintln!("Skipping test: {:?} not found", path);
        return;
    }

    let config = io::read_lammps_data(path, &type_map).unwrap();
    assert_eq!(config.atoms.len(), 1000);
    assert_eq!(config.species.len(), 3);
    assert_eq!(config.composition["Ca"], 200);
    assert_eq!(config.composition["Si"], 200);
    assert_eq!(config.composition["O"], 600);

    // Box should be ~23.86 A
    for d in 0..3 {
        assert!((config.box_lengths[d] - 23.859).abs() < 0.01);
    }
}

#[test]
fn test_rdf_peaks() {
    let mut type_map = HashMap::new();
    type_map.insert(1, "Ca".to_string());
    type_map.insert(2, "Si".to_string());
    type_map.insert(3, "O".to_string());

    let path = Path::new("../run01/CaSiO3_glass.data");
    if !path.exists() {
        eprintln!("Skipping test: {:?} not found", path);
        return;
    }

    let config = io::read_lammps_data(path, &type_map).unwrap();
    let rdfs = rdf::compute_partial_rdfs(&config, 550, 11.0);

    // Si-O first peak should be around 1.6 A
    let si_o_idx = config.pair_index(1, 2); // Si=1, O=2
    let gr = &rdfs.partials[&si_o_idx];
    let max_bin = gr.iter().enumerate().max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap()).unwrap().0;
    let peak_r = rdfs.r[max_bin];
    assert!((peak_r - 1.6).abs() < 0.2, "Si-O peak at {} A, expected ~1.6 A", peak_r);

    // O-O first peak should be around 2.6 A
    let o_o_idx = config.pair_index(2, 2); // O=2
    let gr = &rdfs.partials[&o_o_idx];
    let max_bin = gr.iter().enumerate().max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap()).unwrap().0;
    let peak_r = rdfs.r[max_bin];
    assert!((peak_r - 2.6).abs() < 0.2, "O-O peak at {} A, expected ~2.6 A", peak_r);
}

#[test]
fn test_sq_normalisation() {
    // S(Q) should approach 1.0 at high Q
    let mut type_map = HashMap::new();
    type_map.insert(1, "Ca".to_string());
    type_map.insert(2, "Si".to_string());
    type_map.insert(3, "O".to_string());

    let path = Path::new("../run01/CaSiO3_glass.data");
    if !path.exists() {
        eprintln!("Skipping test: {:?} not found", path);
        return;
    }

    let config = io::read_lammps_data(path, &type_map).unwrap();
    let rho0 = config.number_density();
    let rdfs = rdf::compute_partial_rdfs(&config, 550, 11.0);

    let nq = 500;
    let q: Vec<f64> = (0..nq).map(|i| 0.3 + i as f64 * (20.0 - 0.3) / nq as f64).collect();

    let partial_sq = sq::compute_all_partial_sq(&rdfs.r, &rdfs.partials, rho0, &q, true);
    let sx = xray::compute_xray_sq(&config, &partial_sq, &q);

    // At high Q (>15 1/A), S(Q) should oscillate around 1.0
    let high_q_vals: Vec<f64> = sx.iter().skip(400).cloned().collect();
    let avg: f64 = high_q_vals.iter().sum::<f64>() / high_q_vals.len() as f64;
    assert!(
        (avg - 1.0).abs() < 0.1,
        "High-Q average S(Q) = {}, expected ~1.0",
        avg
    );
}

#[test]
fn test_pbc_distance() {
    let config = Configuration {
        atoms: vec![
            reversesmith::atoms::Atom {
                position: [0.5, 0.5, 0.5],
                species: "A".to_string(),
                type_id: 0,
            },
            reversesmith::atoms::Atom {
                position: [9.5, 9.5, 9.5],
                species: "A".to_string(),
                type_id: 0,
            },
        ],
        box_lengths: [10.0, 10.0, 10.0],
        species: vec!["A".to_string()],
        composition: {
            let mut m = HashMap::new();
            m.insert("A".to_string(), 2);
            m
        },
    };

    // Distance should use minimum image: 0.5 to 9.5 across PBC = 1.0 in each dim
    let d = config.distance(0, 1);
    let expected = (3.0_f64).sqrt(); // sqrt(1^2 + 1^2 + 1^2)
    assert!(
        (d - expected).abs() < 1e-10,
        "PBC distance = {}, expected {}",
        d,
        expected
    );
}
