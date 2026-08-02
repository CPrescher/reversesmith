#![cfg(feature = "gap-quip")]

use std::collections::HashMap;
use std::path::Path;

use rsmith::atoms::{Atom, Configuration};
use rsmith::cells::CellList;
use rsmith::config::{MlBackend, MlEnergyDelta, MlPotentialConfig};
use rsmith::energy::EnergyModel;
use rsmith::ml_potential::GapQuipModel;

fn diamond_silicon(cells: usize) -> Configuration {
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
    let mut atoms = Vec::new();
    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                for fractional in basis {
                    atoms.push(Atom {
                        position: [
                            (x as f64 + fractional[0]) * lattice_constant,
                            (y as f64 + fractional[1]) * lattice_constant,
                            (z as f64 + fractional[2]) * lattice_constant,
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
        box_lengths: [lattice_constant * cells as f64; 3],
        species: vec!["Si".to_string()],
        composition,
    }
}

#[test]
fn local_gap_delta_matches_full_and_is_transactional_when_test_model_is_available() {
    let Ok(model_path) = std::env::var("RSMITH_GAP_TEST_MODEL") else {
        eprintln!("skipping GAP/QUIP integration test: RSMITH_GAP_TEST_MODEL is not set");
        return;
    };
    let init_args = std::env::var("RSMITH_GAP_TEST_INIT_ARGS").ok();
    std::thread::Builder::new()
        .name("gap-quip-integration".to_string())
        .stack_size(64 * 1024 * 1024)
        .spawn(move || {
            let config = diamond_silicon(3);
            let make_cfg = |delta| MlPotentialConfig {
                backend: MlBackend::GapQuip,
                model: Some(model_path.clone()),
                coefficient_file: None,
                parameter_file: None,
                init_args: init_args.clone(),
                weight: Some(1.0),
                cutoff: Some(5.0),
                delta: Some(delta),
                device: None,
                torch_threads: None,
                dtype: None,
                compile_mode: None,
                python: None,
                worker: None,
            };
            let positions = config
                .atoms
                .iter()
                .map(|atom| atom.position)
                .collect::<Vec<_>>();
            let cell_list = CellList::new(&positions, &config.box_lengths, 5.0);
            let atom = 0;
            let old_position = config.atoms[atom].position;
            let new_position = [0.017, config.box_lengths[1] - 0.009, 0.004];

            let mut full =
                GapQuipModel::from_config(&make_cfg(MlEnergyDelta::Full), &config, Path::new("."))
                    .expect("full GAP model should initialize");
            let full_delta = full.energy_delta_atom(
                &config,
                atom,
                &old_position,
                &new_position,
                &cell_list,
                0,
                0,
            );
            full.reject_move(atom, &old_position);

            let mut local =
                GapQuipModel::from_config(&make_cfg(MlEnergyDelta::Local), &config, Path::new("."))
                    .expect("local GAP model should initialize");
            let local_delta = local.energy_delta_atom(
                &config,
                atom,
                &old_position,
                &new_position,
                &cell_list,
                0,
                0,
            );
            let cluster_sizes = local
                .last_local_cluster_sizes()
                .expect("local GAP should report cluster sizes");
            local.reject_move(atom, &old_position);
            let repeated_delta = local.energy_delta_atom(
                &config,
                atom,
                &old_position,
                &new_position,
                &cell_list,
                0,
                0,
            );
            local.reject_move(atom, &old_position);

            let accepted_delta = local.energy_delta_atom(
                &config,
                atom,
                &old_position,
                &new_position,
                &cell_list,
                0,
                0,
            );
            local.accept_move(atom, &new_position);
            let reverse_delta = local.energy_delta_atom(
                &config,
                atom,
                &new_position,
                &old_position,
                &cell_list,
                0,
                0,
            );
            local.reject_move(atom, &new_position);

            assert!(
                (local_delta - full_delta).abs() < 1.0e-8,
                "local GAP delta {local_delta} differs from full delta {full_delta}"
            );
            assert!((repeated_delta - local_delta).abs() < 1.0e-12);
            assert!((accepted_delta + reverse_delta).abs() < 1.0e-8);
            assert!(cluster_sizes.0 < config.atoms.len());
            assert!(cluster_sizes.0 <= cluster_sizes.1);
        })
        .expect("GAP/QUIP integration thread should start")
        .join()
        .expect("GAP/QUIP integration thread should complete");
}
