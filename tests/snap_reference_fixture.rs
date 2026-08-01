use std::path::{Path, PathBuf};

use rsmith::ml_potential::snap_native::SnapModelFiles;
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
