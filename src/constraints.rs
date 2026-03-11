use std::collections::HashMap;

use crate::atoms::Configuration;

/// Minimum interatomic distance constraints per pair type.
/// Key: pair label like "Ca-O", Value: minimum distance in Angstrom.
pub type MinDistanceConstraints = HashMap<String, f64>;

/// Coordination number constraint for a specific pair.
#[derive(Debug, Clone)]
pub struct CoordinationConstraint {
    pub pair: String, // e.g. "Si-O"
    pub min: usize,
    pub max: usize,
    pub cutoff: f64,
}

/// All constraints.
#[derive(Debug, Clone)]
pub struct Constraints {
    pub min_distances: MinDistanceConstraints,
    pub coordination: Vec<CoordinationConstraint>,
}

impl Constraints {
    pub fn new() -> Self {
        Constraints {
            min_distances: HashMap::new(),
            coordination: Vec::new(),
        }
    }

    /// Get minimum distance for a pair given type ids.
    /// Returns 0.0 if no constraint is set for this pair.
    pub fn min_distance_for_pair(&self, config: &Configuration, type_a: usize, type_b: usize) -> f64 {
        let sa = &config.species[type_a];
        let sb = &config.species[type_b];

        // Try both orderings
        let key1 = format!("{}-{}", sa, sb);
        let key2 = format!("{}-{}", sb, sa);

        self.min_distances
            .get(&key1)
            .or_else(|| self.min_distances.get(&key2))
            .copied()
            .unwrap_or(0.0)
    }
}

/// Check if moving atom `atom_idx` to `new_pos` violates minimum distance constraints.
/// Returns true if the move is VALID (no violations).
pub fn check_min_distances(
    config: &Configuration,
    atom_idx: usize,
    new_pos: &[f64; 3],
    constraints: &Constraints,
) -> bool {
    let ti = config.atoms[atom_idx].type_id;

    for (j, atom_j) in config.atoms.iter().enumerate() {
        if j == atom_idx {
            continue;
        }

        let min_d = constraints.min_distance_for_pair(config, ti, atom_j.type_id);
        if min_d <= 0.0 {
            continue;
        }

        // Compute distance from new_pos to atom j
        let mut r2 = 0.0f64;
        for d in 0..3 {
            let mut delta = atom_j.position[d] - new_pos[d];
            let l = config.box_lengths[d];
            delta -= l * (delta / l).round();
            r2 += delta * delta;
        }

        if r2 < min_d * min_d {
            return false;
        }
    }

    true
}

/// Check coordination number constraints after a proposed move.
/// Only checks constraints involving the moved atom's species.
/// Returns true if the move is VALID.
pub fn check_coordination(
    config: &Configuration,
    atom_idx: usize,
    new_pos: &[f64; 3],
    constraints: &Constraints,
) -> bool {
    let moved_species = &config.atoms[atom_idx].species;

    for cc in &constraints.coordination {
        let parts: Vec<&str> = cc.pair.split('-').collect();
        if parts.len() != 2 {
            continue;
        }
        let (sp_a, sp_b) = (parts[0], parts[1]);

        // Only check if moved atom is involved in this constraint
        if moved_species != sp_a && moved_species != sp_b {
            continue;
        }

        // For each atom of type sp_a, count neighbours of type sp_b within cutoff
        // We only need to re-check atoms that could be affected by the move
        let cutoff2 = cc.cutoff * cc.cutoff;

        // Check atoms of type sp_a: their coordination with sp_b
        for (i, atom_i) in config.atoms.iter().enumerate() {
            if &atom_i.species != sp_a {
                continue;
            }

            // Only check if this atom is the moved atom, or if it's a neighbour of the moved atom
            let pos_i = if i == atom_idx {
                *new_pos
            } else {
                atom_i.position
            };

            // Check if this atom is close enough to the moved atom to be affected
            if i != atom_idx {
                // Is this atom within cutoff of either old or new position?
                let old_r2 = min_image_r2(
                    &config.atoms[atom_idx].position,
                    &atom_i.position,
                    &config.box_lengths,
                );
                let new_r2 =
                    min_image_r2(new_pos, &atom_i.position, &config.box_lengths);
                if old_r2 > cutoff2 * 4.0 && new_r2 > cutoff2 * 4.0 {
                    continue; // This atom is far from the moved atom
                }
            }

            // Count sp_b neighbours of atom i
            let mut count = 0usize;
            for (j, atom_j) in config.atoms.iter().enumerate() {
                if j == i {
                    continue;
                }
                if &atom_j.species != sp_b {
                    continue;
                }

                let pos_j = if j == atom_idx {
                    *new_pos
                } else {
                    atom_j.position
                };

                let r2 = min_image_r2(&pos_i, &pos_j, &config.box_lengths);
                if r2 < cutoff2 {
                    count += 1;
                }
            }

            if count < cc.min || count > cc.max {
                return false;
            }
        }
    }

    true
}

fn min_image_r2(a: &[f64; 3], b: &[f64; 3], box_lengths: &[f64; 3]) -> f64 {
    let mut r2 = 0.0f64;
    for d in 0..3 {
        let mut delta = b[d] - a[d];
        let l = box_lengths[d];
        delta -= l * (delta / l).round();
        r2 += delta * delta;
    }
    r2
}
