use std::collections::HashSet;
use std::path::{Path, PathBuf};

use crate::atoms::Configuration;
use crate::cells::CellList;
use crate::config::MlPotentialConfig;
use crate::energy::EnergyModel;

use super::{SnapModelFiles, SnapNeighbor};

const POSITION_TOLERANCE_ANGSTROM: f64 = 1.0e-10;

/// Stateful native SNAP evaluator with a dedicated neighbor list and
/// accepted-state per-atom energy cache.
pub struct SnapNativeModel {
    model: SnapModelFiles,
    weight: f64,
    maximum_cutoff: f64,
    box_lengths: [f64; 3],
    type_indices: Vec<usize>,
    accepted_positions: Vec<[f64; 3]>,
    cell_list: SnapCellList,
    atom_energies: Vec<f64>,
    total_energy: f64,
    pending_trial: Option<PendingTrial>,
}

struct PendingTrial {
    atom_index: usize,
    old_position: [f64; 3],
    new_position: [f64; 3],
    energy_updates: Vec<(usize, f64)>,
    delta: f64,
}

impl SnapNativeModel {
    pub fn from_config(
        config: &MlPotentialConfig,
        configuration: &Configuration,
        base_directory: &Path,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let coefficient_file = config
            .coefficient_file
            .as_deref()
            .ok_or("[ml_potential] snap_native requires coefficient_file")?;
        let parameter_file = config
            .parameter_file
            .as_deref()
            .ok_or("[ml_potential] snap_native requires parameter_file")?;
        let coefficient_path = resolve_path(base_directory, coefficient_file);
        let parameter_path = resolve_path(base_directory, parameter_file);
        let model =
            SnapModelFiles::load(&coefficient_path, &parameter_path, &configuration.species)?;
        Ok(Self::new(
            model,
            configuration,
            config.weight.unwrap_or(0.001),
        )?)
    }

    pub fn new(
        model: SnapModelFiles,
        configuration: &Configuration,
        weight: f64,
    ) -> Result<Self, String> {
        if !weight.is_finite() || weight < 0.0 {
            return Err("native SNAP weight must be finite and non-negative".to_string());
        }
        if configuration.atoms.is_empty() {
            return Err("native SNAP requires at least one atom".to_string());
        }
        if configuration
            .box_lengths
            .iter()
            .any(|length| !length.is_finite() || *length <= 0.0)
        {
            return Err("native SNAP box lengths must be finite and positive".to_string());
        }
        if configuration.atoms.iter().any(|atom| {
            atom.position
                .iter()
                .any(|coordinate| !coordinate.is_finite())
        }) {
            return Err("native SNAP atom positions must be finite".to_string());
        }

        let type_indices = configuration
            .atoms
            .iter()
            .map(|atom| {
                model.element_index_for_type(atom.type_id)?;
                Ok(atom.type_id)
            })
            .collect::<Result<Vec<_>, String>>()?;
        let maximum_cutoff = model.maximum_cutoff_for_types(&type_indices)?;
        for (axis, length) in configuration.box_lengths.iter().enumerate() {
            if 2.0 * maximum_cutoff >= *length {
                return Err(format!(
                    "native SNAP minimum-image evaluation requires box length {axis} ({length}) to exceed twice the maximum cutoff ({})",
                    2.0 * maximum_cutoff
                ));
            }
        }

        let accepted_positions = configuration
            .atoms
            .iter()
            .map(|atom| wrap_position(atom.position, configuration.box_lengths))
            .collect::<Vec<_>>();
        let cell_list = SnapCellList::new(
            &accepted_positions,
            configuration.box_lengths,
            maximum_cutoff,
        );
        let mut evaluator = Self {
            model,
            weight,
            maximum_cutoff,
            box_lengths: configuration.box_lengths,
            type_indices,
            accepted_positions,
            cell_list,
            atom_energies: vec![0.0; configuration.atoms.len()],
            total_energy: 0.0,
            pending_trial: None,
        };
        for atom_index in 0..evaluator.accepted_positions.len() {
            evaluator.atom_energies[atom_index] = evaluator.evaluate_atom(atom_index, None)?;
        }
        evaluator.total_energy = evaluator.atom_energies.iter().sum();
        Ok(evaluator)
    }

    pub fn cached_total_energy(&self) -> f64 {
        self.total_energy
    }

    fn evaluate_atom(
        &self,
        atom_index: usize,
        trial: Option<(usize, [f64; 3])>,
    ) -> Result<f64, String> {
        let central_position = match trial {
            Some((moved_index, new_position)) if moved_index == atom_index => new_position,
            _ => self.accepted_positions[atom_index],
        };
        let mut neighbors = Vec::new();
        for neighbor_index in self.cell_list.candidate_atoms(central_position) {
            if neighbor_index == atom_index
                || trial.is_some_and(|(moved_index, _)| moved_index == neighbor_index)
            {
                continue;
            }
            let displacement = minimum_image_displacement(
                central_position,
                self.accepted_positions[neighbor_index],
                self.box_lengths,
            );
            // This is a coarse cell-list filter. `atomic_energy` applies the
            // exact species-pair cutoff, and `maximum_cutoff` was constructed
            // as the maximum over every configured pair.
            if squared_norm(displacement) <= self.maximum_cutoff.powi(2) {
                neighbors.push(SnapNeighbor {
                    displacement,
                    type_index: self.type_indices[neighbor_index],
                });
            }
        }

        if let Some((moved_index, new_position)) = trial {
            if moved_index != atom_index {
                let displacement =
                    minimum_image_displacement(central_position, new_position, self.box_lengths);
                if squared_norm(displacement) <= self.maximum_cutoff.powi(2) {
                    neighbors.push(SnapNeighbor {
                        displacement,
                        type_index: self.type_indices[moved_index],
                    });
                }
            }
        }

        self.model
            .atomic_energy(self.type_indices[atom_index], &neighbors)
    }

    fn affected_atoms(
        &self,
        atom_index: usize,
        old_position: [f64; 3],
        new_position: [f64; 3],
    ) -> Result<Vec<usize>, String> {
        let mut affected = vec![false; self.accepted_positions.len()];
        affected[atom_index] = true;
        let moved_type = self.type_indices[atom_index];

        for (position, candidates) in [
            (old_position, self.cell_list.candidate_atoms(old_position)),
            (new_position, self.cell_list.candidate_atoms(new_position)),
        ] {
            for candidate in candidates {
                if candidate == atom_index {
                    continue;
                }
                let cutoff = self
                    .model
                    .cutoff_for_types(moved_type, self.type_indices[candidate])?;
                // SNAP pair cutoffs are symmetric because they are
                // rcutfac * (radius_moved + radius_candidate). One check is
                // therefore sufficient for the candidate central atom.
                let displacement = minimum_image_displacement(
                    self.accepted_positions[candidate],
                    position,
                    self.box_lengths,
                );
                if squared_norm(displacement) <= cutoff * cutoff {
                    affected[candidate] = true;
                }
            }
        }

        Ok(affected
            .into_iter()
            .enumerate()
            .filter_map(|(index, is_affected)| is_affected.then_some(index))
            .collect())
    }

    fn require_no_pending_trial(&self) {
        assert!(
            self.pending_trial.is_none(),
            "native SNAP received a new trial before the previous trial was accepted or rejected"
        );
    }

    #[cfg(debug_assertions)]
    fn debug_assert_cache_consistent(&self) {
        let cached_sum = self.atom_energies.iter().sum::<f64>();
        let scale = self
            .atom_energies
            .iter()
            .map(|energy| energy.abs())
            .sum::<f64>()
            .max(1.0);
        debug_assert!(
            (cached_sum - self.total_energy).abs() <= 1.0e-10 * scale,
            "native SNAP total energy cache became inconsistent"
        );
    }
}

impl EnergyModel for SnapNativeModel {
    fn label(&self) -> &str {
        "SNAP/native"
    }

    fn weight(&self) -> f64 {
        self.weight
    }

    fn cutoff(&self) -> f64 {
        self.maximum_cutoff
    }

    fn total_energy(&mut self, config: &Configuration, _cell_list: &CellList) -> f64 {
        self.require_no_pending_trial();
        assert_eq!(
            config.atoms.len(),
            self.accepted_positions.len(),
            "native SNAP configuration atom count changed after initialization"
        );
        self.total_energy
    }

    fn energy_delta_atom(
        &mut self,
        _config: &Configuration,
        atom_idx: usize,
        old_pos: &[f64; 3],
        new_pos: &[f64; 3],
        _cell_list: &CellList,
        _old_cell: usize,
        _new_cell: usize,
    ) -> f64 {
        self.require_no_pending_trial();
        assert!(atom_idx < self.accepted_positions.len());
        assert!(
            old_pos.iter().chain(new_pos).all(|value| value.is_finite()),
            "native SNAP trial positions must be finite"
        );
        // Trial positions may use any periodic image. Wrapping before all
        // cell and minimum-image operations keeps the accepted/trial states
        // comparable at box boundaries.
        let old_position = wrap_position(*old_pos, self.box_lengths);
        let new_position = wrap_position(*new_pos, self.box_lengths);
        let accepted_error = squared_norm(minimum_image_displacement(
            self.accepted_positions[atom_idx],
            old_position,
            self.box_lengths,
        ));
        assert!(
            accepted_error <= POSITION_TOLERANCE_ANGSTROM.powi(2),
            "native SNAP old trial position does not match its accepted state"
        );

        let affected = self
            .affected_atoms(atom_idx, old_position, new_position)
            .unwrap_or_else(|error| panic!("native SNAP failed to find affected atoms: {error}"));
        let old_energy = affected
            .iter()
            .map(|index| self.atom_energies[*index])
            .sum::<f64>();
        let energy_updates = affected
            .into_iter()
            .map(|index| {
                let energy = self
                    .evaluate_atom(index, Some((atom_idx, new_position)))
                    .unwrap_or_else(|error| {
                        panic!("native SNAP failed to evaluate trial environment: {error}")
                    });
                (index, energy)
            })
            .collect::<Vec<_>>();
        let new_energy = energy_updates.iter().map(|(_, energy)| energy).sum::<f64>();
        let delta = new_energy - old_energy;
        self.pending_trial = Some(PendingTrial {
            atom_index: atom_idx,
            old_position,
            new_position,
            energy_updates,
            delta,
        });
        delta
    }

    fn accept_move(&mut self, atom_idx: usize, new_pos: &[f64; 3]) {
        let pending = self
            .pending_trial
            .take()
            .expect("native SNAP accept_move called without a pending trial");
        assert_eq!(pending.atom_index, atom_idx);
        let new_position = wrap_position(*new_pos, self.box_lengths);
        assert!(
            squared_norm(minimum_image_displacement(
                pending.new_position,
                new_position,
                self.box_lengths
            )) <= POSITION_TOLERANCE_ANGSTROM.powi(2),
            "native SNAP accepted position differs from its pending trial"
        );

        self.accepted_positions[atom_idx] = pending.new_position;
        self.cell_list.move_atom(atom_idx, pending.new_position);
        for (index, energy) in pending.energy_updates {
            self.atom_energies[index] = energy;
        }
        self.total_energy += pending.delta;
        #[cfg(debug_assertions)]
        self.debug_assert_cache_consistent();
    }

    fn reject_move(&mut self, atom_idx: usize, old_pos: &[f64; 3]) {
        let pending = self
            .pending_trial
            .take()
            .expect("native SNAP reject_move called without a pending trial");
        assert_eq!(pending.atom_index, atom_idx);
        let old_position = wrap_position(*old_pos, self.box_lengths);
        assert!(
            squared_norm(minimum_image_displacement(
                pending.old_position,
                old_position,
                self.box_lengths
            )) <= POSITION_TOLERANCE_ANGSTROM.powi(2),
            "native SNAP rejected position differs from its accepted state"
        );
    }
}

struct SnapCellList {
    cell_counts: [usize; 3],
    box_lengths: [f64; 3],
    buckets: Vec<Vec<usize>>,
    atom_cells: Vec<usize>,
}

impl SnapCellList {
    fn new(positions: &[[f64; 3]], box_lengths: [f64; 3], cutoff: f64) -> Self {
        let cell_counts =
            std::array::from_fn(|axis| (box_lengths[axis] / cutoff).floor().max(1.0) as usize);
        let bucket_count = cell_counts.iter().product();
        let mut list = Self {
            cell_counts,
            box_lengths,
            buckets: vec![Vec::new(); bucket_count],
            atom_cells: vec![0; positions.len()],
        };
        for (atom_index, position) in positions.iter().enumerate() {
            let cell = list.cell_for_position(*position);
            list.atom_cells[atom_index] = cell;
            list.buckets[cell].push(atom_index);
        }
        list
    }

    fn candidate_atoms(&self, position: [f64; 3]) -> Vec<usize> {
        let center = self.cell_coordinates(self.cell_for_position(position));
        let mut seen_cells = HashSet::with_capacity(27);
        let mut cells = Vec::with_capacity(27);
        for z_offset in -1..=1 {
            for y_offset in -1..=1 {
                for x_offset in -1..=1 {
                    let coordinates = [
                        periodic_cell(center[0], x_offset, self.cell_counts[0]),
                        periodic_cell(center[1], y_offset, self.cell_counts[1]),
                        periodic_cell(center[2], z_offset, self.cell_counts[2]),
                    ];
                    let cell = self.linear_cell(coordinates);
                    if seen_cells.insert(cell) {
                        cells.push(cell);
                    }
                }
            }
        }
        cells
            .into_iter()
            .flat_map(|cell| self.buckets[cell].iter().copied())
            .collect()
    }

    fn move_atom(&mut self, atom_index: usize, new_position: [f64; 3]) {
        let old_cell = self.atom_cells[atom_index];
        let new_cell = self.cell_for_position(new_position);
        if old_cell == new_cell {
            return;
        }
        let offset = self.buckets[old_cell]
            .iter()
            .position(|candidate| *candidate == atom_index)
            .expect("native SNAP cell list lost an atom");
        self.buckets[old_cell].swap_remove(offset);
        self.buckets[new_cell].push(atom_index);
        self.atom_cells[atom_index] = new_cell;
    }

    fn cell_for_position(&self, position: [f64; 3]) -> usize {
        let coordinates = std::array::from_fn(|axis| {
            let wrapped = position[axis].rem_euclid(self.box_lengths[axis]);
            ((wrapped / self.box_lengths[axis] * self.cell_counts[axis] as f64).floor() as usize)
                .min(self.cell_counts[axis] - 1)
        });
        self.linear_cell(coordinates)
    }

    fn cell_coordinates(&self, cell: usize) -> [usize; 3] {
        let xy = self.cell_counts[0] * self.cell_counts[1];
        let z = cell / xy;
        let remainder = cell % xy;
        let y = remainder / self.cell_counts[0];
        let x = remainder % self.cell_counts[0];
        [x, y, z]
    }

    fn linear_cell(&self, coordinates: [usize; 3]) -> usize {
        coordinates[2] * self.cell_counts[0] * self.cell_counts[1]
            + coordinates[1] * self.cell_counts[0]
            + coordinates[0]
    }
}

fn periodic_cell(cell: usize, offset: isize, cell_count: usize) -> usize {
    (cell as isize + offset).rem_euclid(cell_count as isize) as usize
}

fn wrap_position(mut position: [f64; 3], box_lengths: [f64; 3]) -> [f64; 3] {
    for axis in 0..3 {
        position[axis] = position[axis].rem_euclid(box_lengths[axis]);
    }
    position
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

fn squared_norm(vector: [f64; 3]) -> f64 {
    vector.iter().map(|component| component * component).sum()
}

fn resolve_path(base_directory: &Path, path: &str) -> PathBuf {
    let path = Path::new(path);
    if path.is_absolute() {
        path.to_path_buf()
    } else {
        base_directory.join(path)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_cell_periodic_candidate_search_does_not_duplicate_atoms() {
        let positions = [[0.1, 0.1, 0.1], [9.9, 9.9, 9.9], [5.0, 5.0, 5.0]];
        let list = SnapCellList::new(&positions, [10.0; 3], 4.9);
        assert_eq!(list.cell_counts, [2; 3]);
        let mut candidates = list.candidate_atoms([0.0, 0.0, 0.0]);
        candidates.sort_unstable();
        assert_eq!(candidates, [0, 1, 2]);
    }

    #[test]
    fn moving_an_atom_updates_the_dedicated_cell_list() {
        let positions = [[0.1, 0.1, 0.1], [5.1, 0.1, 0.1]];
        let mut list = SnapCellList::new(&positions, [12.0; 3], 3.0);
        let old_cell = list.atom_cells[0];
        list.move_atom(0, [8.1, 0.1, 0.1]);
        assert_ne!(list.atom_cells[0], old_cell);
        assert!(!list.buckets[old_cell].contains(&0));
        assert!(list.buckets[list.atom_cells[0]].contains(&0));
    }

    #[test]
    fn candidate_search_excludes_distant_cells_in_large_boxes() {
        let positions = [
            [1.0, 1.0, 1.0],
            [2.0, 1.0, 1.0],
            [14.0, 14.0, 14.0],
            [25.0, 25.0, 25.0],
        ];
        let list = SnapCellList::new(&positions, [30.0; 3], 4.0);
        let mut candidates = list.candidate_atoms([1.0, 1.0, 1.0]);
        candidates.sort_unstable();
        assert_eq!(candidates, [0, 1]);
    }
}
