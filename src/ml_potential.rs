//! Machine-learning potential support.

use crate::atoms::Configuration;

#[cfg(feature = "gap-quip")]
pub mod gap_quip;
pub mod mace_python;
pub mod snap_native;

#[cfg(feature = "gap-quip")]
pub use gap_quip::GapQuipModel;
pub use mace_python::MacePythonModel;
pub use snap_native::SnapNativeModel;

#[cfg(not(feature = "gap-quip"))]
use crate::config::MlPotentialConfig;

#[cfg(not(feature = "gap-quip"))]
pub fn gap_quip_disabled_error(_cfg: &MlPotentialConfig) -> Box<dyn std::error::Error> {
    "GAP/QUIP support was requested, but rsmith was built without the 'gap-quip' feature".into()
}

/// Return the atoms whose local GAP environments can change when `atom_idx` is moved.
///
/// This includes the moved atom and any atom within `cutoff` of either the old
/// or proposed new position under periodic boundary conditions.
pub fn affected_atoms(
    config: &Configuration,
    atom_idx: usize,
    old_pos: &[f64; 3],
    new_pos: &[f64; 3],
    cutoff: f64,
) -> Vec<usize> {
    let cutoff2 = cutoff * cutoff;
    let mut affected = Vec::new();
    affected.push(atom_idx);

    for (j, atom) in config.atoms.iter().enumerate() {
        if j == atom_idx {
            continue;
        }

        let old_r2 = distance2_pbc(&atom.position, old_pos, &config.box_lengths);
        let new_r2 = distance2_pbc(&atom.position, new_pos, &config.box_lengths);
        if old_r2 < cutoff2 || new_r2 < cutoff2 {
            affected.push(j);
        }
    }

    affected.sort_unstable();
    affected.dedup();
    affected
}

fn distance2_pbc(a: &[f64; 3], b: &[f64; 3], box_lengths: &[f64; 3]) -> f64 {
    let mut r2 = 0.0;
    for d in 0..3 {
        let mut delta = a[d] - b[d];
        let l = box_lengths[d];
        delta -= l * (delta / l).round();
        r2 += delta * delta;
    }
    r2
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use crate::atoms::{Atom, Configuration};

    use super::affected_atoms;

    fn test_config() -> Configuration {
        let atoms = vec![
            Atom {
                position: [0.2, 0.2, 0.2],
                species: "Si".to_string(),
                type_id: 0,
            },
            Atom {
                position: [9.8, 0.2, 0.2],
                species: "O".to_string(),
                type_id: 1,
            },
            Atom {
                position: [5.0, 5.0, 5.0],
                species: "O".to_string(),
                type_id: 1,
            },
            Atom {
                position: [2.1, 0.2, 0.2],
                species: "O".to_string(),
                type_id: 1,
            },
        ];
        let mut composition = HashMap::new();
        composition.insert("Si".to_string(), 1);
        composition.insert("O".to_string(), 3);
        Configuration {
            atoms,
            box_lengths: [10.0, 10.0, 10.0],
            species: vec!["Si".to_string(), "O".to_string()],
            composition,
        }
    }

    #[test]
    fn affected_atoms_use_periodic_boundaries() {
        let config = test_config();
        let affected = affected_atoms(&config, 0, &[0.2, 0.2, 0.2], &[0.5, 0.2, 0.2], 0.6);
        assert_eq!(affected, vec![0, 1]);
    }

    #[test]
    fn affected_atoms_include_old_and_new_neighbor_shells() {
        let config = test_config();
        let affected = affected_atoms(&config, 0, &[0.2, 0.2, 0.2], &[2.0, 0.2, 0.2], 0.35);
        assert_eq!(affected, vec![0, 3]);
    }
}
