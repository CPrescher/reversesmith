//! GAP/QUIP energy backend.
//!
//! This module expects a small C ABI shim to be available at link time. The
//! shim owns the QUIP/libAtoms calculator and exposes per-atom GAP energies.

use std::ffi::CString;
use std::os::raw::{c_char, c_double, c_int};
use std::path::Path;

use crate::atoms::Configuration;
use crate::cells::CellList;
use crate::config::MlPotentialConfig;
use crate::energy::EnergyModel;
use crate::ml_potential::affected_atoms;

#[repr(C)]
struct GapQuipHandle {
    _private: [u8; 0],
}

extern "C" {
    fn rsmith_gap_quip_create(
        model_path: *const c_char,
        init_args: *const c_char,
    ) -> *mut GapQuipHandle;
    fn rsmith_gap_quip_destroy(handle: *mut GapQuipHandle);
    fn rsmith_gap_quip_set_structure(
        handle: *mut GapQuipHandle,
        n_atoms: usize,
        species: *const *const c_char,
        positions: *const c_double,
        box_lengths: *const c_double,
    ) -> c_int;
    fn rsmith_gap_quip_move_atom(
        handle: *mut GapQuipHandle,
        atom_idx: usize,
        position: *const c_double,
    ) -> c_int;
    fn rsmith_gap_quip_per_atom_energy(
        handle: *mut GapQuipHandle,
        n_indices: usize,
        indices: *const usize,
        out_energy: *mut c_double,
    ) -> c_int;
}

pub struct GapQuipModel {
    handle: *mut GapQuipHandle,
    weight: f64,
    cutoff: f64,
    _species_names: Vec<CString>,
}

impl GapQuipModel {
    pub fn from_config(
        cfg: &MlPotentialConfig,
        config: &Configuration,
        base_dir: &Path,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        if cfg.cutoff <= 0.0 {
            return Err("[ml_potential] cutoff must be greater than 0".into());
        }

        let model_path = if Path::new(&cfg.model).is_absolute() {
            Path::new(&cfg.model).to_path_buf()
        } else {
            base_dir.join(&cfg.model)
        };
        if !model_path.exists() {
            return Err(format!("GAP model file not found: {}", model_path.display()).into());
        }

        let model_c = CString::new(model_path.to_string_lossy().as_bytes())?;
        let init_args = cfg.init_args.as_deref().map(CString::new).transpose()?;
        let init_args_ptr = init_args.as_ref().map_or(std::ptr::null(), |s| s.as_ptr());
        let handle = unsafe { rsmith_gap_quip_create(model_c.as_ptr(), init_args_ptr) };
        if handle.is_null() {
            return Err("failed to create GAP/QUIP calculator".into());
        }

        let species_names = config
            .atoms
            .iter()
            .map(|atom| CString::new(atom.species.as_str()))
            .collect::<Result<Vec<_>, _>>()?;
        let species_ptrs: Vec<*const c_char> = species_names.iter().map(|s| s.as_ptr()).collect();
        let positions = flatten_positions(config);
        let rc = unsafe {
            rsmith_gap_quip_set_structure(
                handle,
                config.atoms.len(),
                species_ptrs.as_ptr(),
                positions.as_ptr(),
                config.box_lengths.as_ptr(),
            )
        };
        if rc != 0 {
            unsafe {
                rsmith_gap_quip_destroy(handle);
            }
            return Err(format!("failed to initialize GAP/QUIP structure (code {rc})").into());
        }

        Ok(GapQuipModel {
            handle,
            weight: cfg.weight.unwrap_or(0.001),
            cutoff: cfg.cutoff,
            _species_names: species_names,
        })
    }

    fn move_atom_or_panic(&mut self, atom_idx: usize, pos: &[f64; 3]) {
        let rc = unsafe { rsmith_gap_quip_move_atom(self.handle, atom_idx, pos.as_ptr()) };
        if rc != 0 {
            panic!("GAP/QUIP failed to move atom {atom_idx} (code {rc})");
        }
    }

    fn sum_per_atom_energy_or_panic(&mut self, indices: &[usize]) -> f64 {
        let mut energy = 0.0f64;
        let rc = unsafe {
            rsmith_gap_quip_per_atom_energy(
                self.handle,
                indices.len(),
                indices.as_ptr(),
                &mut energy as *mut f64,
            )
        };
        if rc != 0 {
            panic!("GAP/QUIP failed to compute per-atom energy (code {rc})");
        }
        energy
    }
}

impl EnergyModel for GapQuipModel {
    fn label(&self) -> &str {
        "GAP/QUIP"
    }

    fn weight(&self) -> f64 {
        self.weight
    }

    fn cutoff(&self) -> f64 {
        self.cutoff
    }

    fn total_energy(&mut self, config: &Configuration, _cell_list: &CellList) -> f64 {
        let indices: Vec<usize> = (0..config.atoms.len()).collect();
        self.sum_per_atom_energy_or_panic(&indices)
    }

    fn energy_delta_atom(
        &mut self,
        config: &Configuration,
        atom_idx: usize,
        old_pos: &[f64; 3],
        new_pos: &[f64; 3],
        _cell_list: &CellList,
        _old_cell: usize,
        _new_cell: usize,
    ) -> f64 {
        let affected = affected_atoms(config, atom_idx, old_pos, new_pos, self.cutoff);
        let old_energy = self.sum_per_atom_energy_or_panic(&affected);
        self.move_atom_or_panic(atom_idx, new_pos);
        let new_energy = self.sum_per_atom_energy_or_panic(&affected);
        new_energy - old_energy
    }

    fn accept_move(&mut self, _atom_idx: usize, _new_pos: &[f64; 3]) {
        // The calculator was already moved during energy_delta_atom.
    }

    fn reject_move(&mut self, atom_idx: usize, old_pos: &[f64; 3]) {
        self.move_atom_or_panic(atom_idx, old_pos);
    }
}

impl Drop for GapQuipModel {
    fn drop(&mut self) {
        if !self.handle.is_null() {
            unsafe {
                rsmith_gap_quip_destroy(self.handle);
            }
        }
    }
}

fn flatten_positions(config: &Configuration) -> Vec<f64> {
    let mut positions = Vec::with_capacity(config.atoms.len() * 3);
    for atom in &config.atoms {
        positions.extend_from_slice(&atom.position);
    }
    positions
}
