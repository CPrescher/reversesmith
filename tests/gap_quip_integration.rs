#![cfg(feature = "gap-quip")]

use std::collections::HashMap;
use std::path::Path;

use rsmith::atoms::{Atom, Configuration};
use rsmith::cells::CellList;
use rsmith::config::{MlBackend, MlPotentialConfig};
use rsmith::energy::EnergyModel;
use rsmith::ml_potential::GapQuipModel;

fn small_config() -> Configuration {
    let atoms = vec![
        Atom {
            position: [1.0, 1.0, 1.0],
            species: "Si".to_string(),
            type_id: 0,
        },
        Atom {
            position: [2.6, 1.0, 1.0],
            species: "O".to_string(),
            type_id: 1,
        },
    ];
    let mut composition = HashMap::new();
    composition.insert("Si".to_string(), 1);
    composition.insert("O".to_string(), 1);
    Configuration {
        atoms,
        box_lengths: [10.0, 10.0, 10.0],
        species: vec!["Si".to_string(), "O".to_string()],
        composition,
    }
}

#[test]
fn gap_quip_backend_initializes_when_test_model_is_available() {
    let Ok(model) = std::env::var("RSMITH_GAP_TEST_MODEL") else {
        eprintln!("skipping GAP/QUIP integration test: RSMITH_GAP_TEST_MODEL is not set");
        return;
    };
    let init_args = std::env::var("RSMITH_GAP_TEST_INIT_ARGS").ok();

    let config = small_config();
    let ml_cfg = MlPotentialConfig {
        backend: MlBackend::GapQuip,
        model,
        init_args,
        weight: Some(0.001),
        cutoff: 5.0,
    };
    let mut backend = GapQuipModel::from_config(&ml_cfg, &config, Path::new("."))
        .expect("GAP/QUIP backend should initialize with the configured test model");

    assert_eq!(backend.label(), "GAP/QUIP");
    assert!(backend.cutoff() > 0.0);

    let positions: Vec<[f64; 3]> = config.atoms.iter().map(|atom| atom.position).collect();
    let cell_list = CellList::new(&positions, &config.box_lengths, backend.cutoff());
    let energy = backend.total_energy(&config, &cell_list);
    assert!(energy.is_finite(), "GAP/QUIP energy should be finite");
}
