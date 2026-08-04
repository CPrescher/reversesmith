use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;

use rsmith::atoms::{Atom, Configuration};
use rsmith::cells::CellList;
use rsmith::config::{MlBackend, MlPotentialConfig};
use rsmith::energy::EnergyModel;
use rsmith::ml_potential::pace_native::{PaceModel, PaceNativeModel, PaceNeighbor};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/pace")
        .join(name)
}

fn copper_configuration(distance: f64) -> Configuration {
    let atoms = vec![
        Atom {
            position: [5.0, 5.0, 5.0],
            species: "Cu".to_string(),
            type_id: 0,
        },
        Atom {
            position: [5.0 + distance, 5.0, 5.0],
            species: "Cu".to_string(),
            type_id: 0,
        },
    ];
    Configuration {
        atoms,
        box_lengths: [20.0; 3],
        species: vec!["Cu".to_string()],
        composition: HashMap::from([("Cu".to_string(), 2)]),
    }
}

fn copper_cluster_configuration() -> Configuration {
    let positions = [
        [5.0, 5.0, 5.0],
        [7.2, 5.0, 5.0],
        [5.4, 7.1, 5.3],
        [7.0, 6.8, 6.5],
    ];
    Configuration {
        atoms: positions
            .into_iter()
            .map(|position| Atom {
                position,
                species: "Cu".to_string(),
                type_id: 0,
            })
            .collect(),
        box_lengths: [20.0; 3],
        species: vec!["Cu".to_string()],
        composition: HashMap::from([("Cu".to_string(), 4)]),
    }
}

fn tantalum_cluster_configuration() -> Configuration {
    let mut configuration = copper_cluster_configuration();
    for atom in &mut configuration.atoms {
        atom.species = "Ta".to_string();
    }
    configuration.species = vec!["Ta".to_string()];
    configuration.composition = HashMap::from([("Ta".to_string(), 4)]);
    configuration
}

fn binary_configuration() -> Configuration {
    Configuration {
        atoms: vec![
            Atom {
                position: [5.0, 5.0, 5.0],
                species: "H".to_string(),
                type_id: 0,
            },
            Atom {
                position: [7.0, 5.0, 5.0],
                species: "O".to_string(),
                type_id: 1,
            },
            Atom {
                position: [5.0, 7.2, 5.3],
                species: "H".to_string(),
                type_id: 0,
            },
        ],
        box_lengths: [20.0; 3],
        species: vec!["H".to_string(), "O".to_string()],
        composition: HashMap::from([("H".to_string(), 2), ("O".to_string(), 1)]),
    }
}

#[test]
fn rank_one_chebpow_matches_lammps_oracle() {
    let model = PaceModel::load(&fixture("probe_chebpow.yace"), &["Cu".to_string()]).unwrap();
    // Frozen from LAMMPS 22 Jul 2025 Update 4 + ML-PACE 2023.11.25.
    for (distance, expected_atomic_energy) in [
        (0.5, 0.81),
        (1.0, 0.6400000000000001),
        (2.0, 0.36),
        (4.25, 0.022499999999999964),
        (4.75, 0.0025000000000000577),
        (5.1, 0.0),
    ] {
        let energy = model
            .atomic_energy(
                0,
                &[PaceNeighbor {
                    displacement: [distance, 0.0, 0.0],
                    type_index: 0,
                }],
            )
            .unwrap();
        assert!(
            (energy - expected_atomic_energy).abs() < 2.0e-13,
            "r={distance}: native={energy:.16e}, LAMMPS={expected_atomic_energy:.16e}"
        );
    }
}

#[test]
fn sbessel_rank_one_matches_python_ace_oracle() {
    // Frozen from official ICAMS python-ace commit
    // b143ac3c5d18c55d8a1f0701fae855b2638536fe.
    for (distance, oracle_total_energy) in [
        (0.5, 1.364_228_739_562_039),
        (1.0, 1.2030019100150908),
        (2.0, 0.7042495846821467),
        (4.25, 0.013172351831062924),
        (4.75, 0.0004587334497791528),
        (5.1, 0.0),
    ] {
        let configuration = copper_configuration(distance);
        let model =
            PaceModel::load(&fixture("probe_sbessel.yace"), &configuration.species).unwrap();
        let runtime = PaceNativeModel::new(model, &configuration, 1.0).unwrap();
        let error = (runtime.cached_total_energy() - oracle_total_energy).abs();
        assert!(
            error < 2.0e-11,
            "r={distance}: native={:.16e}, python-ace={oracle_total_energy:.16e}, error={error:.3e}",
            runtime.cached_total_energy()
        );
    }
}

#[test]
fn core_repulsion_and_e0_match_lammps_oracle() {
    for (distance, lammps_total) in [
        (0.5, 24.167_511_919_302_637),
        (1.0, 2.948_238_273_135_793),
        (1.234, 1.160_833_562_657_103_1),
        (2.0, 0.502_195_631_404_510_3),
    ] {
        let configuration = copper_configuration(distance);
        let model = PaceModel::load(&fixture("probe_core.yace"), &configuration.species).unwrap();
        let runtime = PaceNativeModel::new(model, &configuration, 1.0).unwrap();
        let error = (runtime.cached_total_energy() - lammps_total).abs();
        assert!(
            error < 2.0e-11,
            "r={distance}: native={:.16e}, LAMMPS={lammps_total:.16e}, error={error:.3e}",
            runtime.cached_total_energy()
        );
    }
}

#[test]
fn multi_element_channels_and_directed_bonds_match_lammps_oracle() {
    let configuration = binary_configuration();
    let model = PaceModel::load(&fixture("probe_binary.yace"), &configuration.species).unwrap();
    assert_eq!(model.cutoff_for_types(0, 1).unwrap(), 5.0);
    assert_eq!(model.cutoff_for_types(1, 0).unwrap(), 4.5);
    let runtime = PaceNativeModel::new(model, &configuration, 1.0).unwrap();
    let lammps_energy = 0.566_947_475_847_808_1;
    let error = (runtime.cached_total_energy() - lammps_energy).abs();
    assert!(
        error < 2.0e-12,
        "native={:.16e}, LAMMPS={lammps_energy:.16e}, error={error:.3e}",
        runtime.cached_total_energy()
    );
}

#[test]
fn accepted_and_rejected_moves_keep_the_local_cache_transactional() {
    let mut configuration = copper_configuration(2.0);
    let model = PaceModel::load(&fixture("probe_chebpow.yace"), &configuration.species).unwrap();
    let mut runtime = PaceNativeModel::new(model, &configuration, 0.75).unwrap();
    let positions = configuration
        .atoms
        .iter()
        .map(|atom| atom.position)
        .collect::<Vec<_>>();
    let cell_list = CellList::new(&positions, &configuration.box_lengths, 5.0);
    let initial = runtime.total_energy(&configuration, &cell_list);
    assert!((initial - 0.72).abs() < 2.0e-13);
    assert_eq!(runtime.weight(), 0.75);
    assert_eq!(runtime.cutoff(), 5.0);

    let old = configuration.atoms[1].position;
    let accepted = [7.4, 5.0, 5.0];
    configuration.atoms[1].position = accepted;
    let delta = runtime.energy_delta_atom(&configuration, 1, &old, &accepted, &cell_list, 0, 0);
    runtime.accept_move(1, &accepted);
    let rebuilt_model =
        PaceModel::load(&fixture("probe_chebpow.yace"), &configuration.species).unwrap();
    let rebuilt = PaceNativeModel::new(rebuilt_model, &configuration, 0.75).unwrap();
    assert!((initial + delta - rebuilt.cached_total_energy()).abs() < 2.0e-13);
    assert!((runtime.cached_total_energy() - rebuilt.cached_total_energy()).abs() < 2.0e-13);

    let rejected = [8.0, 5.0, 5.0];
    configuration.atoms[1].position = rejected;
    let before_rejection = runtime.cached_total_energy();
    let _ = runtime.energy_delta_atom(&configuration, 1, &accepted, &rejected, &cell_list, 0, 0);
    runtime.reject_move(1, &accepted);
    assert_eq!(runtime.cached_total_energy(), before_rejection);
}

#[test]
fn config_constructor_resolves_a_relative_yace_model() {
    let configuration = copper_configuration(2.0);
    let config = MlPotentialConfig {
        backend: MlBackend::PaceNative,
        model: Some("probe_chebpow.yace".to_string()),
        coefficient_file: None,
        parameter_file: None,
        init_args: None,
        weight: Some(0.125),
        cutoff: None,
        delta: None,
        device: None,
        torch_threads: None,
        dtype: None,
        compile_mode: None,
        python: None,
        worker: None,
    };
    let base_directory = fixture("");
    let mut runtime =
        PaceNativeModel::from_config(&config, &configuration, &base_directory).unwrap();
    let positions = configuration
        .atoms
        .iter()
        .map(|atom| atom.position)
        .collect::<Vec<_>>();
    let cell_list = CellList::new(&positions, &configuration.box_lengths, 5.0);
    assert!((runtime.total_energy(&configuration, &cell_list) - 0.72).abs() < 2.0e-13);
    assert_eq!(runtime.weight(), 0.125);
}

#[test]
fn nonlinear_rank_two_copper_model_matches_lammps_oracle() {
    let mut configuration = copper_cluster_configuration();
    let model = PaceModel::load(&fixture("probe_rank2.yace"), &configuration.species).unwrap();
    let mut runtime = PaceNativeModel::new(model, &configuration, 1.0).unwrap();
    // LAMMPS 22 Jul 2025 Update 4, pair_style pace product, ML-PACE 2023.11.25.
    let lammps_energy = 6.466_977_752_567_922;
    let error = (runtime.cached_total_energy() - lammps_energy).abs();
    assert!(
        error < 2.0e-11,
        "native={:.16e}, LAMMPS={lammps_energy:.16e}, error={error:.3e}",
        runtime.cached_total_energy()
    );

    let positions = configuration
        .atoms
        .iter()
        .map(|atom| atom.position)
        .collect::<Vec<_>>();
    let cell_list = CellList::new(&positions, &configuration.box_lengths, 3.9);
    let old_position = configuration.atoms[0].position;
    let new_position = [5.13, 4.93, 5.08];
    configuration.atoms[0].position = new_position;
    let delta = runtime.energy_delta_atom(
        &configuration,
        0,
        &old_position,
        &new_position,
        &cell_list,
        0,
        0,
    );
    let rebuilt_model =
        PaceModel::load(&fixture("probe_rank2.yace"), &configuration.species).unwrap();
    let rebuilt = PaceNativeModel::new(rebuilt_model, &configuration, 1.0).unwrap();
    assert!((lammps_energy + delta - rebuilt.cached_total_energy()).abs() < 2.0e-11);
    runtime.accept_move(0, &new_position);
    assert!((runtime.cached_total_energy() - rebuilt.cached_total_energy()).abs() < 2.0e-11);
}

#[test]
fn rank_four_chebexpcos_model_matches_lammps_oracle() {
    let configuration = tantalum_cluster_configuration();
    let model = PaceModel::load(&fixture("probe_rank4.yace"), &configuration.species).unwrap();
    let runtime = PaceNativeModel::new(model, &configuration, 1.0).unwrap();
    // LAMMPS 22 Jul 2025 Update 4, pair_style pace product, ML-PACE 2023.11.25.
    let lammps_energy = 8.868_259_922_646_295;
    let error = (runtime.cached_total_energy() - lammps_energy).abs();
    assert!(
        error < 2.0e-10,
        "native={:.16e}, LAMMPS={lammps_energy:.16e}, error={error:.3e}",
        runtime.cached_total_energy()
    );
}

#[test]
#[ignore = "requires a LAMMPS binary with ML-PACE built in or RSMITH_LAMMPS_PACE_PLUGIN_DIR"]
fn live_lammps_product_evaluator_matches_frozen_copper_reference() {
    let executable = std::env::var("RSMITH_LAMMPS").unwrap_or_else(|_| "lmp_serial".to_string());
    let input = fixture("lammps_cu_cluster.in");
    let model = fixture("probe_rank2.yace");
    let mut command = Command::new(executable);
    command
        .arg("-log")
        .arg("none")
        .arg("-in")
        .arg(input)
        .arg("-var")
        .arg("model")
        .arg(model);
    if let Ok(plugin_directory) = std::env::var("RSMITH_LAMMPS_PACE_PLUGIN_DIR") {
        command.env("LAMMPS_PLUGIN_PATH", plugin_directory);
    }
    let output = command.output().expect("failed to launch LAMMPS");
    assert!(
        output.status.success(),
        "LAMMPS failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).unwrap();
    let marker = stdout
        .lines()
        .find(|line| line.starts_with("PACE_CU_CLUSTER"))
        .expect("LAMMPS output did not contain PACE_CU_CLUSTER marker");
    let energy = marker
        .split_whitespace()
        .find_map(|field| field.strip_prefix("pe="))
        .unwrap()
        .parse::<f64>()
        .unwrap();
    assert!((energy - 6.466_977_752_567_922).abs() < 2.0e-11);
}
