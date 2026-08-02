use std::collections::HashMap;
use std::env;
use std::path::Path;

use rsmith::atoms::{Atom, Configuration};
use rsmith::cells::CellList;
use rsmith::config::{MaceDtype, MlBackend, MlEnergyDelta, MlPotentialConfig};
use rsmith::energy::EnergyModel;
use rsmith::ml_potential::MacePythonModel;

fn diamond_si_config() -> Configuration {
    diamond_si_supercell(1)
}

fn diamond_si_supercell(cells: usize) -> Configuration {
    let a = 5.43;
    let frac = [
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.0, 0.5, 0.5],
        [0.25, 0.75, 0.75],
        [0.5, 0.0, 0.5],
        [0.75, 0.25, 0.75],
        [0.5, 0.5, 0.0],
        [0.75, 0.75, 0.25],
    ];
    let mut atoms = Vec::with_capacity(8 * cells.pow(3));
    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                for f in frac {
                    atoms.push(Atom {
                        position: [
                            a * (x as f64 + f[0]),
                            a * (y as f64 + f[1]),
                            a * (z as f64 + f[2]),
                        ],
                        species: "Si".to_string(),
                        type_id: 0,
                    });
                }
            }
        }
    }
    let mut composition = HashMap::new();
    composition.insert("Si".to_string(), atoms.len());
    Configuration {
        atoms,
        box_lengths: [a * cells as f64; 3],
        species: vec!["Si".to_string()],
        composition,
    }
}

fn mace_test_cfg(model_path: String, delta: Option<MlEnergyDelta>) -> MlPotentialConfig {
    let python = env::var("RSMITH_MACE_TEST_PYTHON").unwrap_or_else(|_| "python3".to_string());
    let device = env::var("RSMITH_MACE_TEST_DEVICE").unwrap_or_else(|_| "cpu".to_string());
    let torch_threads = env::var("RSMITH_MACE_TEST_TORCH_THREADS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .or(Some(1));
    let dtype = match env::var("RSMITH_MACE_TEST_DTYPE").as_deref() {
        Ok("float32") => Some(MaceDtype::Float32),
        Ok("float64") => Some(MaceDtype::Float64),
        Ok(value) => panic!("unsupported RSMITH_MACE_TEST_DTYPE: {value}"),
        Err(_) => None,
    };

    MlPotentialConfig {
        backend: MlBackend::MacePython,
        model: Some(model_path),
        coefficient_file: None,
        parameter_file: None,
        init_args: None,
        weight: Some(0.001),
        cutoff: Some(5.0),
        delta,
        device: Some(device),
        torch_threads,
        dtype,
        compile_mode: None,
        python: Some(python),
        worker: None,
    }
}

#[test]
fn mace_python_backend_runs_real_model_when_configured() {
    let Ok(model_path) = env::var("RSMITH_MACE_TEST_MODEL") else {
        eprintln!("skipping real MACE test because RSMITH_MACE_TEST_MODEL is not set");
        return;
    };
    if !Path::new(&model_path).exists() {
        panic!("RSMITH_MACE_TEST_MODEL does not exist: {model_path}");
    }

    let mut config = diamond_si_config();
    let positions = config
        .atoms
        .iter()
        .map(|atom| atom.position)
        .collect::<Vec<_>>();
    let cell_list = CellList::new(&positions, &config.box_lengths, 5.0);
    let cfg = mace_test_cfg(model_path, None);

    let mut model = MacePythonModel::from_config(&cfg, &config, Path::new(".")).unwrap();
    let initial = model.total_energy(&config, &cell_list);
    assert!(initial.is_finite(), "initial MACE energy is not finite");

    let atom_idx = 0;
    let old_pos = config.atoms[atom_idx].position;
    let new_pos = [old_pos[0] + 0.05, old_pos[1] + 0.02, old_pos[2] + 0.01];
    config.atoms[atom_idx].position = new_pos;
    let delta = model.energy_delta_atom(&config, atom_idx, &old_pos, &new_pos, &cell_list, 0, 0);
    assert!(delta.is_finite(), "MACE move delta is not finite");

    let moved = model.total_energy(&config, &cell_list);
    assert!(
        (moved - initial - delta).abs() < 1.0e-5,
        "delta does not match full energy difference: initial={initial}, moved={moved}, delta={delta}"
    );

    config.atoms[atom_idx].position = old_pos;
    model.reject_move(atom_idx, &old_pos);
    let reverted = model.total_energy(&config, &cell_list);
    assert!(
        (reverted - initial).abs() < 1.0e-5,
        "reject did not restore worker energy: initial={initial}, reverted={reverted}"
    );
}

#[test]
fn mace_python_local_delta_matches_full_delta_when_configured() {
    let Ok(model_path) = env::var("RSMITH_MACE_TEST_MODEL") else {
        eprintln!("skipping real MACE local delta test because RSMITH_MACE_TEST_MODEL is not set");
        return;
    };
    if !Path::new(&model_path).exists() {
        panic!("RSMITH_MACE_TEST_MODEL does not exist: {model_path}");
    }

    let mut config = diamond_si_supercell(2);
    let positions = config
        .atoms
        .iter()
        .map(|atom| atom.position)
        .collect::<Vec<_>>();
    let cell_list = CellList::new(&positions, &config.box_lengths, 5.0);
    let atom_idx = 0;
    let old_pos = config.atoms[atom_idx].position;
    let new_pos = [old_pos[0] + 0.03, old_pos[1] + 0.02, old_pos[2] + 0.01];

    let mut full_model = MacePythonModel::from_config(
        &mace_test_cfg(model_path.clone(), None),
        &config,
        Path::new("."),
    )
    .unwrap();
    let mut local_model = MacePythonModel::from_config(
        &mace_test_cfg(model_path.clone(), Some(MlEnergyDelta::Local)),
        &config,
        Path::new("."),
    )
    .unwrap();
    let mut incremental_model = MacePythonModel::from_config(
        &mace_test_cfg(model_path, Some(MlEnergyDelta::Incremental)),
        &config,
        Path::new("."),
    )
    .unwrap();

    let initial_full = full_model.total_energy(&config, &cell_list);
    let initial_local = local_model.total_energy(&config, &cell_list);
    let initial_incremental = incremental_model.total_energy(&config, &cell_list);
    assert!(
        (initial_local - initial_full).abs() < 1.0e-5,
        "local cache has the wrong initial energy: local={initial_local}, full={initial_full}"
    );
    assert!(
        (initial_incremental - initial_full).abs() < 1.0e-5,
        "incremental cache has the wrong initial energy: incremental={initial_incremental}, full={initial_full}"
    );

    let full_delta =
        full_model.energy_delta_atom(&config, atom_idx, &old_pos, &new_pos, &cell_list, 0, 0);
    let local_delta =
        local_model.energy_delta_atom(&config, atom_idx, &old_pos, &new_pos, &cell_list, 0, 0);
    let incremental_delta = incremental_model
        .energy_delta_atom(&config, atom_idx, &old_pos, &new_pos, &cell_list, 0, 0);
    assert!(
        (local_delta - full_delta).abs() < 1.0e-3,
        "local MACE delta differs from full delta: local={local_delta}, full={full_delta}"
    );
    assert!(
        (incremental_delta - full_delta).abs() < 1.0e-3,
        "incremental MACE delta differs from full delta: incremental={incremental_delta}, full={full_delta}"
    );

    full_model.reject_move(atom_idx, &old_pos);
    local_model.reject_move(atom_idx, &old_pos);
    incremental_model.reject_move(atom_idx, &old_pos);
    assert!(
        (local_model.total_energy(&config, &cell_list) - initial_local).abs() < 1.0e-5,
        "reject changed the accepted local energy cache"
    );
    assert!(
        (incremental_model.total_energy(&config, &cell_list) - initial_incremental).abs() < 1.0e-5,
        "reject changed the accepted incremental energy cache"
    );

    let accepted_pos = [old_pos[0] + 0.02, old_pos[1] - 0.01, old_pos[2] + 0.015];
    let accepted_full_delta =
        full_model.energy_delta_atom(&config, atom_idx, &old_pos, &accepted_pos, &cell_list, 0, 0);
    let accepted_local_delta =
        local_model.energy_delta_atom(&config, atom_idx, &old_pos, &accepted_pos, &cell_list, 0, 0);
    let accepted_incremental_delta = incremental_model.energy_delta_atom(
        &config,
        atom_idx,
        &old_pos,
        &accepted_pos,
        &cell_list,
        0,
        0,
    );
    assert!(
        (accepted_local_delta - accepted_full_delta).abs() < 1.0e-3,
        "accepted local MACE delta differs from full delta: local={accepted_local_delta}, full={accepted_full_delta}"
    );
    assert!(
        (accepted_incremental_delta - accepted_full_delta).abs() < 1.0e-3,
        "accepted incremental MACE delta differs from full delta: incremental={accepted_incremental_delta}, full={accepted_full_delta}"
    );
    config.atoms[atom_idx].position = accepted_pos;
    full_model.accept_move(atom_idx, &accepted_pos);
    local_model.accept_move(atom_idx, &accepted_pos);
    incremental_model.accept_move(atom_idx, &accepted_pos);
    let accepted_full = full_model.total_energy(&config, &cell_list);
    let accepted_local = local_model.total_energy(&config, &cell_list);
    let accepted_incremental = incremental_model.total_energy(&config, &cell_list);
    assert!(
        (accepted_local - accepted_full).abs() < 1.0e-3,
        "accepted local energy cache differs from full energy: local={accepted_local}, full={accepted_full}"
    );
    assert!(
        (accepted_incremental - accepted_full).abs() < 1.0e-3,
        "accepted incremental energy cache differs from full energy: incremental={accepted_incremental}, full={accepted_full}"
    );

    let second_atom = 1;
    let second_old = config.atoms[second_atom].position;
    let second_new = [
        second_old[0] - 0.01,
        second_old[1] + 0.02,
        second_old[2] + 0.01,
    ];
    let second_full_delta = full_model.energy_delta_atom(
        &config,
        second_atom,
        &second_old,
        &second_new,
        &cell_list,
        0,
        0,
    );
    let second_local_delta = local_model.energy_delta_atom(
        &config,
        second_atom,
        &second_old,
        &second_new,
        &cell_list,
        0,
        0,
    );
    let second_incremental_delta = incremental_model.energy_delta_atom(
        &config,
        second_atom,
        &second_old,
        &second_new,
        &cell_list,
        0,
        0,
    );
    assert!(
        (second_local_delta - second_full_delta).abs() < 1.0e-3,
        "local cache is stale after accept: local={second_local_delta}, full={second_full_delta}"
    );
    assert!(
        (second_incremental_delta - second_full_delta).abs() < 1.0e-3,
        "incremental cache is stale after accept: incremental={second_incremental_delta}, full={second_full_delta}"
    );
    full_model.reject_move(second_atom, &second_old);
    local_model.reject_move(second_atom, &second_old);
    incremental_model.reject_move(second_atom, &second_old);
}
