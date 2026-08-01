use std::path::{Path, PathBuf};

use rsmith::ml_potential::snap_native::{SnapModelFiles, SnapNeighbor};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct ReferenceFixture {
    schema_version: u32,
    generator: ReferenceGenerator,
    model: ReferenceModel,
    cell: ReferenceCell,
    configurations: Vec<ReferenceConfiguration>,
}

#[derive(Debug, Deserialize)]
struct ReferenceGenerator {
    program: String,
    version: String,
    input: String,
}

#[derive(Debug, Deserialize)]
struct ReferenceModel {
    descriptor_count: usize,
    maximum_cutoff_angstrom: f64,
}

#[derive(Debug, Deserialize)]
struct ReferenceCell {
    atom_count: usize,
    box_lengths_angstrom: [f64; 3],
}

#[derive(Debug, Deserialize)]
struct ReferenceConfiguration {
    name: String,
    total_energy_ev: f64,
    all_atoms_equivalent: bool,
    representative_atoms: Vec<ReferenceAtom>,
}

#[derive(Debug, Deserialize)]
struct ReferenceAtom {
    id: usize,
    position_angstrom: [f64; 3],
    descriptors: Vec<f64>,
    energy_ev: f64,
}

fn test_data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/snap")
        .join(name)
}

fn load_fixture() -> ReferenceFixture {
    serde_json::from_str(
        &std::fs::read_to_string(test_data("linear_two_element_si.reference.json")).unwrap(),
    )
    .unwrap()
}

fn assert_close(actual: f64, expected: f64) {
    let tolerance = 1.0e-12 * expected.abs().max(1.0) + 1.0e-15;
    assert!(
        (actual - expected).abs() <= tolerance,
        "expected {expected:.16e}, got {actual:.16e}"
    );
}

fn diamond_supercell_positions() -> Vec<[f64; 3]> {
    let lattice_constant = 5.43;
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
    let mut positions = Vec::with_capacity(64);
    for z_cell in 0..2 {
        for y_cell in 0..2 {
            for x_cell in 0..2 {
                for fractional in basis {
                    positions.push([
                        (x_cell as f64 + fractional[0]) * lattice_constant,
                        (y_cell as f64 + fractional[1]) * lattice_constant,
                        (z_cell as f64 + fractional[2]) * lattice_constant,
                    ]);
                }
            }
        }
    }
    positions
}

fn minimum_image_displacement(
    central: [f64; 3],
    neighbor: [f64; 3],
    box_lengths: [f64; 3],
) -> [f64; 3] {
    std::array::from_fn(|axis| {
        let displacement = neighbor[axis] - central[axis];
        displacement - box_lengths[axis] * (displacement / box_lengths[axis]).round()
    })
}

fn neighbors_for_atom(
    atom_index: usize,
    positions: &[[f64; 3]],
    box_lengths: [f64; 3],
) -> Vec<SnapNeighbor> {
    positions
        .iter()
        .enumerate()
        .filter(|(neighbor_index, _)| *neighbor_index != atom_index)
        .map(|(_, position)| SnapNeighbor {
            displacement: minimum_image_displacement(positions[atom_index], *position, box_lengths),
            type_index: 0,
        })
        .collect()
}

#[test]
fn reference_fixture_matches_the_synthetic_model_contract() {
    let fixture = load_fixture();
    let model = SnapModelFiles::load(
        &test_data("linear_two_element.snapcoeff"),
        &test_data("linear_two_element.snapparam"),
        &["Si".to_string()],
    )
    .unwrap();

    assert_eq!(fixture.schema_version, 1);
    assert_eq!(fixture.generator.program, "LAMMPS");
    assert!(!fixture.generator.version.is_empty());
    assert_eq!(fixture.generator.input, "linear_two_element_si.lammps.in");
    assert!(test_data(&fixture.generator.input).is_file());
    assert_eq!(
        fixture.model.descriptor_count,
        model.coefficients.ncoeff - 1
    );
    assert_eq!(
        fixture
            .configurations
            .iter()
            .map(|configuration| configuration.name.as_str())
            .collect::<Vec<_>>(),
        ["diamond_equilibrium", "diamond_atom_1_displaced"]
    );
    assert!(fixture
        .configurations
        .iter()
        .all(|configuration| configuration.total_energy_ev.is_finite()));

    let maximum_model_cutoff = model.parameters.rcutfac
        * model
            .coefficients
            .elements
            .iter()
            .map(|element| 2.0 * element.radius)
            .fold(0.0_f64, f64::max);
    assert_close(fixture.model.maximum_cutoff_angstrom, maximum_model_cutoff);
    assert!(fixture
        .cell
        .box_lengths_angstrom
        .iter()
        .all(|length| *length > 2.0 * maximum_model_cutoff));
}

#[test]
fn reference_atomic_energies_are_the_linear_descriptor_contraction() {
    let fixture = load_fixture();
    let model = SnapModelFiles::load(
        &test_data("linear_two_element.snapcoeff"),
        &test_data("linear_two_element.snapparam"),
        &["Si".to_string()],
    )
    .unwrap();
    let coefficients = &model.coefficients.elements[0].coefficients;

    for configuration in &fixture.configurations {
        assert!(!configuration.name.is_empty());
        for atom in &configuration.representative_atoms {
            assert!((1..=fixture.cell.atom_count).contains(&atom.id));
            assert!(atom.position_angstrom.iter().all(|value| value.is_finite()));
            assert_eq!(atom.descriptors.len(), fixture.model.descriptor_count);
            let energy = coefficients[0]
                + coefficients[1..]
                    .iter()
                    .zip(&atom.descriptors)
                    .map(|(coefficient, descriptor)| coefficient * descriptor)
                    .sum::<f64>();
            assert_close(energy, atom.energy_ev);
        }
    }

    let equilibrium = fixture
        .configurations
        .iter()
        .find(|configuration| configuration.all_atoms_equivalent)
        .unwrap();
    assert_eq!(equilibrium.representative_atoms.len(), 1);
    assert_close(
        equilibrium.representative_atoms[0].energy_ev * fixture.cell.atom_count as f64,
        equilibrium.total_energy_ev,
    );

    let displaced = fixture
        .configurations
        .iter()
        .find(|configuration| configuration.name == "diamond_atom_1_displaced")
        .unwrap();
    assert!(!displaced.all_atoms_equivalent);
    // Atom 1 is the moved center, 2/3/5 sample distinct perturbed neighbor
    // environments, and 6 is outside the moved atom's cutoff and remains an
    // equilibrium environment. The selected subset does not sum to the total.
    assert_eq!(
        displaced
            .representative_atoms
            .iter()
            .map(|atom| atom.id)
            .collect::<Vec<_>>(),
        [1, 2, 3, 5, 6]
    );
}

#[test]
fn native_descriptors_and_energies_match_lammps() {
    let fixture = load_fixture();
    let model = SnapModelFiles::load(
        &test_data("linear_two_element.snapcoeff"),
        &test_data("linear_two_element.snapparam"),
        &["Si".to_string()],
    )
    .unwrap();

    for configuration in &fixture.configurations {
        let mut positions = diamond_supercell_positions();
        if configuration.name == "diamond_atom_1_displaced" {
            positions[0] = configuration.representative_atoms[0].position_angstrom;
        }

        for expected_atom in &configuration.representative_atoms {
            let atom_index = expected_atom.id - 1;
            let neighbors =
                neighbors_for_atom(atom_index, &positions, fixture.cell.box_lengths_angstrom);
            let descriptors = model.atomic_descriptors(0, &neighbors).unwrap();
            for (actual, expected) in descriptors.iter().zip(&expected_atom.descriptors) {
                assert_close(*actual, *expected);
            }
            assert_close(
                model.atomic_energy(0, &neighbors).unwrap(),
                expected_atom.energy_ev,
            );
        }

        let total_energy = (0..positions.len())
            .map(|atom_index| {
                let neighbors =
                    neighbors_for_atom(atom_index, &positions, fixture.cell.box_lengths_angstrom);
                model.atomic_energy(0, &neighbors).unwrap()
            })
            .sum::<f64>();
        assert_close(total_energy, configuration.total_energy_ev);
    }
}

#[test]
fn native_descriptors_are_rotation_invariant() {
    let model = SnapModelFiles::load(
        &test_data("linear_two_element.snapcoeff"),
        &test_data("linear_two_element.snapparam"),
        &["Si".to_string()],
    )
    .unwrap();
    let positions = diamond_supercell_positions();
    let neighbors = neighbors_for_atom(0, &positions, [10.86; 3]);
    let rotated = neighbors
        .iter()
        .map(|neighbor| SnapNeighbor {
            displacement: [
                neighbor.displacement[1],
                neighbor.displacement[2],
                neighbor.displacement[0],
            ],
            type_index: neighbor.type_index,
        })
        .collect::<Vec<_>>();

    let expected = model.atomic_descriptors(0, &neighbors).unwrap();
    let actual = model.atomic_descriptors(0, &rotated).unwrap();
    for (actual, expected) in actual.iter().zip(expected) {
        assert_close(*actual, expected);
    }
}
