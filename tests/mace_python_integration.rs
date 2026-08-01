use std::collections::HashMap;
use std::env;
use std::path::Path;

use rsmith::atoms::{Atom, Configuration};
use rsmith::cells::CellList;
use rsmith::config::{MlBackend, MlPotentialConfig};
use rsmith::energy::EnergyModel;
use rsmith::ml_potential::MacePythonModel;

fn diamond_si_config() -> Configuration {
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
    let atoms = frac
        .iter()
        .map(|f| Atom {
            position: [a * f[0], a * f[1], a * f[2]],
            species: "Si".to_string(),
            type_id: 0,
        })
        .collect::<Vec<_>>();
    let mut composition = HashMap::new();
    composition.insert("Si".to_string(), atoms.len());
    Configuration {
        atoms,
        box_lengths: [a, a, a],
        species: vec!["Si".to_string()],
        composition,
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
    let python = env::var("RSMITH_MACE_TEST_PYTHON").unwrap_or_else(|_| "python3".to_string());
    let device = env::var("RSMITH_MACE_TEST_DEVICE").unwrap_or_else(|_| "cpu".to_string());
    let torch_threads = env::var("RSMITH_MACE_TEST_TORCH_THREADS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .or(Some(1));

    let cfg = MlPotentialConfig {
        backend: MlBackend::MacePython,
        model: model_path,
        init_args: None,
        weight: Some(0.001),
        cutoff: 5.0,
        device: Some(device),
        torch_threads,
        python: Some(python),
        worker: None,
    };

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
