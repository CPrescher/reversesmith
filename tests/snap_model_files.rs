use std::path::{Path, PathBuf};

use rsmith::ml_potential::snap_native::{SnapModelFiles, SnapNeighbor};

fn test_data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/snap")
        .join(name)
}

fn diamond_supercell_neighbors() -> Vec<SnapNeighbor> {
    let lattice_constant = 5.43;
    let box_length = 2.0 * lattice_constant;
    let basis = [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.25, 0.75, 0.75],
        [0.75, 0.25, 0.75],
        [0.75, 0.75, 0.25],
    ];
    let mut neighbors = Vec::with_capacity(63);
    for z_cell in 0..2 {
        for y_cell in 0..2 {
            for x_cell in 0..2 {
                for fractional in basis {
                    let position = [
                        (x_cell as f64 + fractional[0]) * lattice_constant,
                        (y_cell as f64 + fractional[1]) * lattice_constant,
                        (z_cell as f64 + fractional[2]) * lattice_constant,
                    ];
                    if position == [0.0, 0.0, 0.0] {
                        continue;
                    }
                    neighbors.push(SnapNeighbor {
                        displacement: std::array::from_fn(|axis| {
                            position[axis] - box_length * (position[axis] / box_length).round()
                        }),
                        type_index: 0,
                    });
                }
            }
        }
    }
    neighbors
}

#[test]
fn loads_synthetic_multi_element_model_and_mapping() {
    let elements = vec!["O".to_string(), "Si".to_string(), "O".to_string()];
    let model = SnapModelFiles::load(
        &test_data("linear_two_element.snapcoeff"),
        &test_data("linear_two_element.snapparam"),
        &elements,
    )
    .unwrap();

    assert_eq!(model.type_to_element, vec![1, 0, 1]);
    assert_eq!(model.coefficients.elements[0].name, "Si");
    assert_eq!(model.coefficients.ncoeff, 6);
    assert_eq!(model.parameters.twojmax, 2);
    assert_eq!(model.descriptor_count(), 5);
}

#[test]
fn evaluator_reports_mutated_model_invariants_instead_of_truncating_or_panicking() {
    let mut model = SnapModelFiles::load(
        &test_data("linear_two_element.snapcoeff"),
        &test_data("linear_two_element.snapparam"),
        &["Si".to_string()],
    )
    .unwrap();

    model.type_to_element[0] = usize::MAX;
    let error = model.atomic_energy(0, &[]).unwrap_err();
    assert!(error.contains("maps to invalid element index"), "{error}");

    model.type_to_element[0] = 0;
    model.coefficients.elements[0].coefficients.pop();
    let error = model.atomic_energy(0, &[]).unwrap_err();
    assert!(error.contains("coefficients"), "{error}");
    assert!(error.contains("are required"), "{error}");
}

#[test]
fn evaluator_rejects_numerically_zero_neighbor_distances() {
    let model = SnapModelFiles::load(
        &test_data("linear_two_element.snapcoeff"),
        &test_data("linear_two_element.snapparam"),
        &["Si".to_string()],
    )
    .unwrap();
    let error = model
        .atomic_descriptors(
            0,
            &[SnapNeighbor {
                displacement: [1.0e-200, 0.0, 0.0],
                type_index: 0,
            }],
        )
        .unwrap_err();
    assert!(error.contains("distance must be at least"), "{error}");
}

#[test]
fn loads_lammps_example_models_when_available() {
    let Some(directory) = std::env::var_os("RSMITH_LAMMPS_POTENTIALS") else {
        eprintln!("skipping LAMMPS SNAP parser examples: RSMITH_LAMMPS_POTENTIALS is not set");
        return;
    };
    let directory = PathBuf::from(directory);

    let silicon = SnapModelFiles::load(
        &directory.join("Si_Zuo_JPCA2020.snapcoeff"),
        &directory.join("Si_Zuo_JPCA2020.snapparam"),
        &["Si".to_string()],
    )
    .unwrap();
    assert_eq!(silicon.coefficients.elements[0].name, "Si");
    assert_eq!(silicon.coefficients.ncoeff, 56);
    // LAMMPS 22 Jul 2025 Update 4, pair_style snap, 2x2x2 diamond-Si.
    let total_energy = 64.0
        * silicon
            .atomic_energy(0, &diamond_supercell_neighbors())
            .unwrap();
    assert!((total_energy - -346.549_395_358_721_5).abs() < 1.0e-9);

    let inp = SnapModelFiles::load(
        &directory.join("InP_JCPA2020.snapcoeff"),
        &directory.join("InP_JCPA2020.snapparam"),
        &["In".to_string(), "P".to_string()],
    )
    .unwrap();
    assert!(inp.parameters.chemflag);
    assert_eq!(inp.type_to_element, vec![0, 1]);
}
