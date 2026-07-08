//! Energy regularizers for RMC acceptance.
//!
//! The RMC loop only needs energies and trial-move energy differences. Pair
//! potentials and machine-learning potentials can therefore share this small
//! interface while keeping their own state and evaluation strategy.

use crate::atoms::Configuration;
use crate::cells::CellList;

pub trait EnergyModel {
    /// Human-readable backend name for logging.
    fn label(&self) -> &str;

    /// Weight multiplying the raw energy contribution in the RMC cost.
    fn weight(&self) -> f64;

    /// Interaction/locality cutoff in Angstrom.
    fn cutoff(&self) -> f64;

    /// Total energy of the current configuration.
    fn total_energy(&mut self, config: &Configuration, cell_list: &CellList) -> f64;

    /// Energy change for a trial move of one atom.
    ///
    /// The RMC loop may already have placed `atom_idx` at `new_pos` in
    /// `config`, so implementations must use the explicit old/new positions
    /// rather than assuming `config.atoms[atom_idx]` is the old position.
    fn energy_delta_atom(
        &mut self,
        config: &Configuration,
        atom_idx: usize,
        old_pos: &[f64; 3],
        new_pos: &[f64; 3],
        cell_list: &CellList,
        old_cell: usize,
        new_cell: usize,
    ) -> f64;

    /// Commit any backend state after an accepted move.
    fn accept_move(&mut self, _atom_idx: usize, _new_pos: &[f64; 3]) {}

    /// Revert any backend state after a rejected move.
    fn reject_move(&mut self, _atom_idx: usize, _old_pos: &[f64; 3]) {}
}
