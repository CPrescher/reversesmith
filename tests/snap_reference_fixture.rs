use std::collections::HashMap;
use std::path::{Path, PathBuf};

use rsmith::atoms::{Atom, Configuration};
use rsmith::cells::CellList;
use rsmith::config::{MlBackend, MlPotentialConfig};
use rsmith::energy::EnergyModel;
use rsmith::ml_potential::snap_native::{SnapModelFiles, SnapNativeModel, SnapNeighbor};
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

fn load_chemical_fixture() -> ReferenceFixture {
    serde_json::from_str(
        &std::fs::read_to_string(test_data("chemical_two_element_inp.reference.json")).unwrap(),
    )
    .unwrap()
}

fn load_quadratic_fixture() -> ReferenceFixture {
    serde_json::from_str(
        &std::fs::read_to_string(test_data("quadratic_si.reference.json")).unwrap(),
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
    let type_indices = vec![0; positions.len()];
    typed_neighbors_for_atom(atom_index, positions, &type_indices, box_lengths)
}

fn typed_neighbors_for_atom(
    atom_index: usize,
    positions: &[[f64; 3]],
    type_indices: &[usize],
    box_lengths: [f64; 3],
) -> Vec<SnapNeighbor> {
    positions
        .iter()
        .enumerate()
        .filter(|(neighbor_index, _)| *neighbor_index != atom_index)
        .map(|(neighbor_index, position)| SnapNeighbor {
            displacement: minimum_image_displacement(positions[atom_index], *position, box_lengths),
            type_index: type_indices[neighbor_index],
        })
        .collect()
}

fn zincblende_inp_positions() -> Vec<[f64; 3]> {
    let lattice_constant = 5.87;
    [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.25, 0.75, 0.75],
        [0.75, 0.25, 0.75],
        [0.75, 0.75, 0.25],
    ]
    .into_iter()
    .map(|fractional| fractional.map(|coordinate| coordinate * lattice_constant))
    .collect()
}

fn inp_configuration(positions: &[[f64; 3]]) -> Configuration {
    let type_indices = [0, 0, 0, 0, 1, 1, 1, 1];
    Configuration {
        atoms: positions
            .iter()
            .zip(type_indices)
            .map(|(position, type_id)| Atom {
                position: *position,
                species: if type_id == 0 { "In" } else { "P" }.to_string(),
                type_id,
            })
            .collect(),
        box_lengths: [5.87; 3],
        species: vec!["In".to_string(), "P".to_string()],
        composition: HashMap::from([("In".to_string(), 4), ("P".to_string(), 4)]),
    }
}

fn silicon_configuration(positions: &[[f64; 3]]) -> Configuration {
    silicon_configuration_in_box(positions, [10.86; 3])
}

fn silicon_configuration_in_box(positions: &[[f64; 3]], box_lengths: [f64; 3]) -> Configuration {
    Configuration {
        atoms: positions
            .iter()
            .map(|position| Atom {
                position: *position,
                species: "Si".to_string(),
                type_id: 0,
            })
            .collect(),
        box_lengths,
        species: vec!["Si".to_string()],
        composition: HashMap::from([("Si".to_string(), positions.len())]),
    }
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

#[test]
fn chemical_descriptors_energies_and_local_trials_match_lammps() {
    let fixture = load_chemical_fixture();
    let model_files = SnapModelFiles::load(
        &test_data("chemical_two_element.snapcoeff"),
        &test_data("chemical_two_element.snapparam"),
        &["In".to_string(), "P".to_string()],
    )
    .unwrap();
    assert_eq!(model_files.descriptor_count(), 40);
    assert_eq!(fixture.model.descriptor_count, 40);
    let type_indices = [0, 0, 0, 0, 1, 1, 1, 1];

    for expected_configuration in &fixture.configurations {
        let mut positions = zincblende_inp_positions();
        if expected_configuration.name == "zincblende_atom_1_displaced" {
            positions[0] = expected_configuration.representative_atoms[0].position_angstrom;
        }

        let expected_atom = &expected_configuration.representative_atoms[0];
        let atom_index = expected_atom.id - 1;
        let neighbors = typed_neighbors_for_atom(
            atom_index,
            &positions,
            &type_indices,
            fixture.cell.box_lengths_angstrom,
        );
        let descriptors = model_files
            .atomic_descriptors(type_indices[atom_index], &neighbors)
            .unwrap();
        assert_eq!(descriptors.len(), expected_atom.descriptors.len());
        for (actual, expected) in descriptors.iter().zip(&expected_atom.descriptors) {
            assert_close(*actual, *expected);
        }
        assert_close(
            model_files
                .atomic_energy(type_indices[atom_index], &neighbors)
                .unwrap(),
            expected_atom.energy_ev,
        );

        let total_energy = positions
            .iter()
            .enumerate()
            .map(|(index, _)| {
                let neighbors = typed_neighbors_for_atom(
                    index,
                    &positions,
                    &type_indices,
                    fixture.cell.box_lengths_angstrom,
                );
                model_files
                    .atomic_energy(type_indices[index], &neighbors)
                    .unwrap()
            })
            .sum::<f64>();
        assert_close(total_energy, expected_configuration.total_energy_ev);
    }

    let positions = zincblende_inp_positions();
    let mut configuration = inp_configuration(&positions);
    let mut runtime = SnapNativeModel::new(model_files, &configuration, 1.0).unwrap();
    let rdf_cell_list = CellList::new(&positions, &configuration.box_lengths, 2.6);
    let equilibrium_energy = fixture.configurations[0].total_energy_ev;
    let displaced_energy = fixture.configurations[1].total_energy_ev;
    assert_close(runtime.cached_total_energy(), equilibrium_energy);

    let old_position = configuration.atoms[0].position;
    let new_position = fixture.configurations[1].representative_atoms[0].position_angstrom;
    configuration.atoms[0].position = new_position;
    let delta = runtime.energy_delta_atom(
        &configuration,
        0,
        &old_position,
        &new_position,
        &rdf_cell_list,
        0,
        0,
    );
    assert_close(equilibrium_energy + delta, displaced_energy);
    runtime.accept_move(0, &new_position);
    assert_close(runtime.cached_total_energy(), displaced_energy);

    let center_self_model = SnapModelFiles::load(
        &test_data("chemical_two_element.snapcoeff"),
        &test_data("chemical_two_element_center_self.snapparam"),
        &["In".to_string(), "P".to_string()],
    )
    .unwrap();
    // LAMMPS 22 Jul 2025 Update 4, generated by
    // chemical_two_element_center_self.lammps.in.
    for (configuration_index, expected_total) in [-12.003_907_971_852_996, -12.004_551_450_816_276]
        .into_iter()
        .enumerate()
    {
        let mut positions = zincblende_inp_positions();
        if configuration_index == 1 {
            positions[0] = fixture.configurations[1].representative_atoms[0].position_angstrom;
        }
        let total = positions
            .iter()
            .enumerate()
            .map(|(index, _)| {
                let neighbors = typed_neighbors_for_atom(
                    index,
                    &positions,
                    &type_indices,
                    fixture.cell.box_lengths_angstrom,
                );
                center_self_model
                    .atomic_energy(type_indices[index], &neighbors)
                    .unwrap()
            })
            .sum::<f64>();
        assert_close(total, expected_total);
    }
}

#[test]
fn quadratic_energy_contraction_and_local_trials_match_lammps() {
    let fixture = load_quadratic_fixture();
    let model_files = SnapModelFiles::load(
        &test_data("quadratic_si.snapcoeff"),
        &test_data("quadratic_si.snapparam"),
        &["Si".to_string()],
    )
    .unwrap();
    assert!(model_files.parameters.quadraticflag);
    assert_eq!(model_files.descriptor_count(), 5);
    assert_eq!(model_files.coefficients.ncoeff, 21);

    for expected_configuration in &fixture.configurations {
        let mut positions = diamond_supercell_positions();
        if expected_configuration.name == "diamond_atom_1_displaced" {
            positions[0] = expected_configuration.representative_atoms[0].position_angstrom;
        }
        let expected_atom = &expected_configuration.representative_atoms[0];
        let atom_index = expected_atom.id - 1;
        let neighbors =
            neighbors_for_atom(atom_index, &positions, fixture.cell.box_lengths_angstrom);
        let descriptors = model_files.atomic_descriptors(0, &neighbors).unwrap();
        for (actual, expected) in descriptors.iter().zip(&expected_atom.descriptors) {
            assert_close(*actual, *expected);
        }
        assert_close(
            model_files.atomic_energy(0, &neighbors).unwrap(),
            expected_atom.energy_ev,
        );

        let total_energy = positions
            .iter()
            .enumerate()
            .map(|(index, _)| {
                model_files
                    .atomic_energy(
                        0,
                        &neighbors_for_atom(index, &positions, fixture.cell.box_lengths_angstrom),
                    )
                    .unwrap()
            })
            .sum::<f64>();
        assert_close(total_energy, expected_configuration.total_energy_ev);
    }

    let positions = diamond_supercell_positions();
    let mut configuration = silicon_configuration(&positions);
    let mut runtime = SnapNativeModel::new(model_files, &configuration, 1.0).unwrap();
    let rdf_cell_list = CellList::new(&positions, &configuration.box_lengths, 4.0);
    let equilibrium_energy = fixture.configurations[0].total_energy_ev;
    let displaced_energy = fixture.configurations[1].total_energy_ev;
    let old_position = configuration.atoms[0].position;
    let new_position = fixture.configurations[1].representative_atoms[0].position_angstrom;
    configuration.atoms[0].position = new_position;
    let delta = runtime.energy_delta_atom(
        &configuration,
        0,
        &old_position,
        &new_position,
        &rdf_cell_list,
        0,
        0,
    );
    assert_close(equilibrium_energy + delta, displaced_energy);
    runtime.accept_move(0, &new_position);
    assert_close(runtime.cached_total_energy(), displaced_energy);
}

#[test]
fn quadratic_contraction_composes_with_chemical_descriptors() {
    let positions = zincblende_inp_positions();
    let type_indices = [0, 0, 0, 0, 1, 1, 1, 1];
    let neighbors = typed_neighbors_for_atom(0, &positions, &type_indices, [5.87; 3]);
    let mut model = SnapModelFiles::load(
        &test_data("chemical_two_element.snapcoeff"),
        &test_data("chemical_two_element.snapparam"),
        &["In".to_string(), "P".to_string()],
    )
    .unwrap();
    let descriptors = model.atomic_descriptors(0, &neighbors).unwrap();
    let linear_energy = model.atomic_energy(0, &neighbors).unwrap();
    let quadratic_count = descriptors.len() * (descriptors.len() + 1) / 2;
    model.parameters.quadraticflag = true;
    model.coefficients.ncoeff += quadratic_count;
    for element in &mut model.coefficients.elements {
        element
            .coefficients
            .extend(std::iter::repeat_n(0.0, quadratic_count));
    }

    let quadratic_start = descriptors.len() + 1;
    let descriptor_count = descriptors.len();
    let coefficients = &mut model.coefficients.elements[0].coefficients;
    coefficients[quadratic_start] = 2.0;
    coefficients[quadratic_start + 1] = 3.0;
    coefficients[quadratic_start + descriptor_count] = 4.0;
    let expected = linear_energy
        + descriptors[0].powi(2)
        + 3.0 * descriptors[0] * descriptors[1]
        + 2.0 * descriptors[1].powi(2);
    assert_close(model.atomic_energy(0, &neighbors).unwrap(), expected);
}

#[test]
fn cached_local_trials_match_lammps_and_are_transactional() {
    let fixture = load_fixture();
    let positions = diamond_supercell_positions();
    let mut configuration = silicon_configuration(&positions);
    let model_files = SnapModelFiles::load(
        &test_data("linear_two_element.snapcoeff"),
        &test_data("linear_two_element.snapparam"),
        &["Si".to_string()],
    )
    .unwrap();
    let mut model = SnapNativeModel::new(model_files, &configuration, 0.75).unwrap();
    let rdf_cell_list = CellList::new(&positions, &configuration.box_lengths, 4.0);
    let equilibrium_energy = fixture.configurations[0].total_energy_ev;
    let displaced_energy = fixture.configurations[1].total_energy_ev;
    assert_close(
        model.total_energy(&configuration, &rdf_cell_list),
        equilibrium_energy,
    );

    let old_position = configuration.atoms[0].position;
    let new_position = fixture.configurations[1].representative_atoms[0].position_angstrom;
    configuration.atoms[0].position = new_position;
    let delta = model.energy_delta_atom(
        &configuration,
        0,
        &old_position,
        &new_position,
        &rdf_cell_list,
        0,
        0,
    );
    assert_close(equilibrium_energy + delta, displaced_energy);
    model.reject_move(0, &old_position);
    configuration.atoms[0].position = old_position;
    assert_close(model.cached_total_energy(), equilibrium_energy);

    configuration.atoms[0].position = new_position;
    let delta = model.energy_delta_atom(
        &configuration,
        0,
        &old_position,
        &new_position,
        &rdf_cell_list,
        0,
        0,
    );
    model.accept_move(0, &new_position);
    assert_close(equilibrium_energy + delta, displaced_energy);
    assert_close(model.cached_total_energy(), displaced_energy);

    configuration.atoms[0].position = old_position;
    let reverse_delta = model.energy_delta_atom(
        &configuration,
        0,
        &new_position,
        &old_position,
        &rdf_cell_list,
        0,
        0,
    );
    model.accept_move(0, &old_position);
    assert_close(displaced_energy + reverse_delta, equilibrium_energy);
    assert_close(model.cached_total_energy(), equilibrium_energy);
}

#[test]
fn config_constructor_loads_relative_fit_snap_files() {
    let positions = diamond_supercell_positions();
    let mut configuration = silicon_configuration(&positions);
    let config = MlPotentialConfig {
        backend: MlBackend::SnapNative,
        model: None,
        coefficient_file: Some("linear_two_element.snapcoeff".to_string()),
        parameter_file: Some("linear_two_element.snapparam".to_string()),
        init_args: None,
        weight: Some(0.75),
        cutoff: None,
        device: None,
        torch_threads: None,
        python: None,
        worker: None,
    };
    let base_directory = test_data("");
    let mut model = SnapNativeModel::from_config(&config, &configuration, &base_directory).unwrap();
    let rdf_cell_list = CellList::new(&positions, &configuration.box_lengths, 4.0);
    let equilibrium_energy = load_fixture().configurations[0].total_energy_ev;

    assert_close(
        model.total_energy(&configuration, &rdf_cell_list),
        equilibrium_energy,
    );
    assert_close(model.weight(), 0.75);

    let old_position = configuration.atoms[0].position;
    let new_position = [
        old_position[0] + 0.08,
        old_position[1] - 0.03,
        old_position[2] + 0.04,
    ];
    configuration.atoms[0].position = new_position;
    let accepted_delta = model.energy_delta_atom(
        &configuration,
        0,
        &old_position,
        &new_position,
        &rdf_cell_list,
        0,
        0,
    );
    model.accept_move(0, &new_position);

    let rebuilt = SnapNativeModel::from_config(&config, &configuration, &base_directory).unwrap();
    assert_close(model.cached_total_energy(), rebuilt.cached_total_energy());
    assert_close(
        equilibrium_energy + accepted_delta,
        rebuilt.cached_total_energy(),
    );

    let rejected_position = [
        new_position[0] - 0.02,
        new_position[1] + 0.06,
        new_position[2],
    ];
    configuration.atoms[0].position = rejected_position;
    let rejected_delta = model.energy_delta_atom(
        &configuration,
        0,
        &new_position,
        &rejected_position,
        &rdf_cell_list,
        0,
        0,
    );
    assert!(rejected_delta.is_finite());
    model.reject_move(0, &new_position);
    configuration.atoms[0].position = new_position;
    assert_close(model.cached_total_energy(), rebuilt.cached_total_energy());
}

#[test]
fn local_trial_delta_matches_a_fresh_full_rebuild_in_a_larger_cell() {
    let mut positions = Vec::new();
    for z in 0..5 {
        for y in 0..5 {
            for x in 0..5 {
                positions.push([
                    0.5 + 3.4 * x as f64,
                    0.5 + 3.4 * y as f64,
                    0.5 + 3.4 * z as f64,
                ]);
            }
        }
    }
    let mut configuration = silicon_configuration_in_box(&positions, [17.0; 3]);
    let model_files = SnapModelFiles::load(
        &test_data("linear_two_element.snapcoeff"),
        &test_data("linear_two_element.snapparam"),
        &["Si".to_string()],
    )
    .unwrap();
    let mut local_model = SnapNativeModel::new(model_files, &configuration, 1.0).unwrap();
    let rdf_cell_list = CellList::new(&positions, &configuration.box_lengths, 4.0);
    let initial_energy = local_model.cached_total_energy();
    let old_position = configuration.atoms[0].position;
    let new_position = [0.82, 0.37, 0.71];
    configuration.atoms[0].position = new_position;

    let local_delta = local_model.energy_delta_atom(
        &configuration,
        0,
        &old_position,
        &new_position,
        &rdf_cell_list,
        0,
        0,
    );
    let rebuilt_files = SnapModelFiles::load(
        &test_data("linear_two_element.snapcoeff"),
        &test_data("linear_two_element.snapparam"),
        &["Si".to_string()],
    )
    .unwrap();
    let rebuilt_model = SnapNativeModel::new(rebuilt_files, &configuration, 1.0).unwrap();
    assert_close(
        initial_energy + local_delta,
        rebuilt_model.cached_total_energy(),
    );

    local_model.accept_move(0, &new_position);
    assert_close(
        local_model.cached_total_energy(),
        rebuilt_model.cached_total_energy(),
    );
}

#[test]
#[ignore = "manual release-mode local-delta performance comparison"]
fn benchmark_local_trials_against_full_rebuilds() {
    let mut positions = Vec::new();
    for z in 0..5 {
        for y in 0..5 {
            for x in 0..5 {
                positions.push([
                    0.5 + 3.4 * x as f64,
                    0.5 + 3.4 * y as f64,
                    0.5 + 3.4 * z as f64,
                ]);
            }
        }
    }
    let mut configuration = silicon_configuration_in_box(&positions, [17.0; 3]);
    let model_files = SnapModelFiles::load(
        &test_data("linear_two_element.snapcoeff"),
        &test_data("linear_two_element.snapparam"),
        &["Si".to_string()],
    )
    .unwrap();
    let mut local_model = SnapNativeModel::new(model_files.clone(), &configuration, 1.0).unwrap();
    let rdf_cell_list = CellList::new(&positions, &configuration.box_lengths, 4.0);
    let old_position = configuration.atoms[0].position;
    let trial_count = 200;

    let local_start = std::time::Instant::now();
    let mut local_checksum = 0.0;
    for trial in 0..trial_count {
        let offset = 0.001 * (trial + 1) as f64;
        let new_position = [
            old_position[0] + offset,
            old_position[1] - 0.5 * offset,
            old_position[2] + 0.25 * offset,
        ];
        configuration.atoms[0].position = new_position;
        local_checksum += local_model.energy_delta_atom(
            &configuration,
            0,
            &old_position,
            &new_position,
            &rdf_cell_list,
            0,
            0,
        );
        local_model.reject_move(0, &old_position);
    }
    let local_elapsed = local_start.elapsed();

    let rebuild_start = std::time::Instant::now();
    let mut rebuild_checksum = 0.0;
    for trial in 0..trial_count {
        let offset = 0.001 * (trial + 1) as f64;
        configuration.atoms[0].position = [
            old_position[0] + offset,
            old_position[1] - 0.5 * offset,
            old_position[2] + 0.25 * offset,
        ];
        let rebuilt = SnapNativeModel::new(model_files.clone(), &configuration, 1.0).unwrap();
        rebuild_checksum += rebuilt.cached_total_energy();
    }
    let rebuild_elapsed = rebuild_start.elapsed();

    assert!(local_checksum.is_finite());
    assert!(rebuild_checksum.is_finite());
    eprintln!(
        "{trial_count} trials: local={local_elapsed:?}, full rebuild={rebuild_elapsed:?}, speedup={:.1}x",
        rebuild_elapsed.as_secs_f64() / local_elapsed.as_secs_f64()
    );
}
