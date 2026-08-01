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

fn zincblende_supercell_neighbors(central_index: usize) -> Vec<SnapNeighbor> {
    let lattice_constant = 5.87;
    let box_length = 3.0 * lattice_constant;
    let basis = [
        ([0.0, 0.0, 0.0], 0),
        ([0.0, 0.5, 0.5], 0),
        ([0.5, 0.0, 0.5], 0),
        ([0.5, 0.5, 0.0], 0),
        ([0.25, 0.25, 0.25], 1),
        ([0.25, 0.75, 0.75], 1),
        ([0.75, 0.25, 0.75], 1),
        ([0.75, 0.75, 0.25], 1),
    ];
    let mut atoms = Vec::with_capacity(216);
    for z_cell in 0..3 {
        for y_cell in 0..3 {
            for x_cell in 0..3 {
                for (fractional, type_index) in basis {
                    atoms.push((
                        [
                            (x_cell as f64 + fractional[0]) * lattice_constant,
                            (y_cell as f64 + fractional[1]) * lattice_constant,
                            (z_cell as f64 + fractional[2]) * lattice_constant,
                        ],
                        type_index,
                    ));
                }
            }
        }
    }
    let central_position = atoms[central_index].0;
    atoms
        .into_iter()
        .enumerate()
        .filter(|(index, _)| *index != central_index)
        .map(|(_, (position, type_index))| SnapNeighbor {
            displacement: std::array::from_fn(|axis| {
                let displacement = position[axis] - central_position[axis];
                displacement - box_length * (displacement / box_length).round()
            }),
            type_index,
        })
        .collect()
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

    let silicon_quadratic = SnapModelFiles::load(
        &directory.join("Si_Zuo_JPCA2020.quadratic.snapcoeff"),
        &directory.join("Si_Zuo_JPCA2020.quadratic.snapparam"),
        &["Si".to_string()],
    )
    .unwrap();
    assert!(silicon_quadratic.parameters.quadraticflag);
    assert_eq!(silicon_quadratic.coefficients.ncoeff, 1596);
    let total_energy = 64.0
        * silicon_quadratic
            .atomic_energy(0, &diamond_supercell_neighbors())
            .unwrap();
    assert!((total_energy - -347.062_817_967_136_6).abs() < 1.0e-9);

    let inp = SnapModelFiles::load(
        &directory.join("InP_JCPA2020.snapcoeff"),
        &directory.join("InP_JCPA2020.snapparam"),
        &["In".to_string(), "P".to_string()],
    )
    .unwrap();
    assert!(inp.parameters.chemflag);
    assert_eq!(inp.type_to_element, vec![0, 1]);
    assert_eq!(inp.descriptor_count(), 240);
    // LAMMPS 22 Jul 2025 Update 4, pure pair_style snap (without the ZBL
    // overlay in the distributed example), 3x3x3 zincblende InP.
    let in_energy = inp
        .atomic_energy(0, &zincblende_supercell_neighbors(0))
        .unwrap();
    let p_energy = inp
        .atomic_energy(1, &zincblende_supercell_neighbors(4))
        .unwrap();
    let total_energy = 108.0 * (in_energy + p_energy);
    assert!((total_energy - -1_229.893_944_643_025).abs() < 1.0e-8);
}
