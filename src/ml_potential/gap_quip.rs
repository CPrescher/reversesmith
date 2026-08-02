//! GAP/QUIP energy backend.
//!
//! This module expects a small C ABI shim to be available at link time. The
//! shim owns the QUIP/libAtoms calculator and exposes per-atom GAP energies.

use std::collections::HashSet;
use std::ffi::CString;
use std::os::raw::{c_char, c_double, c_int};
use std::path::Path;

use crate::atoms::Configuration;
use crate::cells::CellList;
use crate::config::{MlEnergyDelta, MlPotentialConfig};
use crate::energy::EnergyModel;

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
    fn rsmith_gap_quip_set_local_cluster(
        handle: *mut GapQuipHandle,
        n_atoms: usize,
        n_central: usize,
        atom_ids: *const usize,
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
    fn rsmith_gap_quip_per_atom_energies(
        handle: *mut GapQuipHandle,
        n_indices: usize,
        indices: *const usize,
        out_energies: *mut c_double,
    ) -> c_int;
}

pub struct GapQuipModel {
    handle: *mut GapQuipHandle,
    weight: f64,
    cutoff: f64,
    delta: MlEnergyDelta,
    _species_names: Vec<CString>,
    local_builder: Option<GapLocalClusterBuilder>,
    accepted_local_energies: Option<Vec<f64>>,
    pending_local: Option<PendingLocalDelta>,
    last_local_cluster_sizes: Option<(usize, usize)>,
}

struct PendingLocalDelta {
    central_atoms: Vec<usize>,
    local_energies: Vec<f64>,
}

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
struct PeriodicImage {
    atom: usize,
    shift: [i32; 3],
}

struct GapLocalCluster {
    images: Vec<PeriodicImage>,
    n_central: usize,
}

struct GapLocalClusterBuilder {
    positions: Vec<[f64; 3]>,
    box_lengths: [f64; 3],
    cutoff: f64,
}

impl GapLocalClusterBuilder {
    fn new(config: &Configuration, cutoff: f64) -> Result<Self, String> {
        if config.atoms.is_empty() {
            return Err("local GAP requires at least one atom".to_string());
        }
        if config
            .box_lengths
            .iter()
            .any(|length| !length.is_finite() || *length <= 0.0)
        {
            return Err("local GAP box lengths must be finite and positive".to_string());
        }
        Ok(Self {
            positions: config
                .atoms
                .iter()
                .map(|atom| wrap_position(atom.position, config.box_lengths))
                .collect(),
            box_lengths: config.box_lengths,
            cutoff,
        })
    }

    fn trial_cluster(
        &self,
        moved_atom: usize,
        old_position: [f64; 3],
        new_position: [f64; 3],
    ) -> Result<(Vec<usize>, GapLocalCluster), String> {
        if moved_atom >= self.positions.len() {
            return Err(format!("local GAP atom index {moved_atom} is out of range"));
        }
        let old_position = wrap_position(old_position, self.box_lengths);
        let new_position = wrap_position(new_position, self.box_lengths);
        if squared_norm(minimum_image_displacement(
            self.positions[moved_atom],
            old_position,
            self.box_lengths,
        )) > 1.0e-20
        {
            return Err("local GAP old position does not match its accepted state".to_string());
        }

        let cutoff2 = self.cutoff * self.cutoff;
        let central_atoms = self
            .positions
            .iter()
            .enumerate()
            .filter_map(|(atom, position)| {
                (atom == moved_atom
                    || squared_norm(minimum_image_displacement(
                        *position,
                        old_position,
                        self.box_lengths,
                    )) < cutoff2
                    || squared_norm(minimum_image_displacement(
                        *position,
                        new_position,
                        self.box_lengths,
                    )) < cutoff2)
                    .then_some(atom)
            })
            .collect::<Vec<_>>();

        let central_images = central_atoms
            .iter()
            .map(|&atom| PeriodicImage {
                atom,
                shift: [0; 3],
            })
            .collect::<Vec<_>>();
        let mut context = HashSet::new();
        for central in &central_images {
            let central_base = if central.atom == moved_atom {
                new_position
            } else {
                self.positions[central.atom]
            };
            let central_position = std::array::from_fn(|axis| {
                central_base[axis] + central.shift[axis] as f64 * self.box_lengths[axis]
            });
            for atom in 0..self.positions.len() {
                let position = if atom == moved_atom {
                    new_position
                } else {
                    self.positions[atom]
                };
                self.insert_images_within(&mut context, atom, position, central_position);
            }
        }
        for central in &central_images {
            context.remove(central);
        }
        let mut context = context.into_iter().collect::<Vec<_>>();
        context.sort_unstable();
        let n_central = central_images.len();
        let mut images = central_images;
        images.extend(context);
        Ok((central_atoms, GapLocalCluster { images, n_central }))
    }

    fn insert_images_within(
        &self,
        images: &mut HashSet<PeriodicImage>,
        atom: usize,
        atom_position: [f64; 3],
        central_position: [f64; 3],
    ) {
        let minimum_shift: [i32; 3] = std::array::from_fn(|axis| {
            ((central_position[axis] - atom_position[axis] - self.cutoff) / self.box_lengths[axis])
                .ceil() as i32
        });
        let maximum_shift: [i32; 3] = std::array::from_fn(|axis| {
            ((central_position[axis] - atom_position[axis] + self.cutoff) / self.box_lengths[axis])
                .floor() as i32
        });
        for shift_x in minimum_shift[0]..=maximum_shift[0] {
            for shift_y in minimum_shift[1]..=maximum_shift[1] {
                for shift_z in minimum_shift[2]..=maximum_shift[2] {
                    let shift = [shift_x, shift_y, shift_z];
                    let displacement = std::array::from_fn(|axis| {
                        atom_position[axis] + shift[axis] as f64 * self.box_lengths[axis]
                            - central_position[axis]
                    });
                    if squared_norm(displacement) <= self.cutoff * self.cutoff + 1.0e-12 {
                        images.insert(PeriodicImage { atom, shift });
                    }
                }
            }
        }
    }

    fn cluster_arrays(
        &self,
        cluster: &GapLocalCluster,
        moved_atom: usize,
        new_position: [f64; 3],
    ) -> (Vec<usize>, Vec<f64>) {
        let mut atoms = Vec::with_capacity(cluster.images.len());
        let mut positions = Vec::with_capacity(3 * cluster.images.len());
        for image in &cluster.images {
            atoms.push(image.atom);
            let base = if image.atom == moved_atom {
                new_position
            } else {
                self.positions[image.atom]
            };
            for axis in 0..3 {
                positions.push(base[axis] + image.shift[axis] as f64 * self.box_lengths[axis]);
            }
        }
        (atoms, positions)
    }

    fn accept_move(&mut self, atom: usize, new_position: [f64; 3]) {
        self.positions[atom] = wrap_position(new_position, self.box_lengths);
    }
}

impl GapQuipModel {
    pub fn from_config(
        cfg: &MlPotentialConfig,
        config: &Configuration,
        base_dir: &Path,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let cutoff = cfg
            .cutoff
            .ok_or("[ml_potential] cutoff is required for GAP/QUIP")?;
        if cutoff <= 0.0 {
            return Err("[ml_potential] cutoff must be greater than 0".into());
        }
        let delta = cfg.delta.unwrap_or(MlEnergyDelta::Full);
        if matches!(delta, MlEnergyDelta::Incremental) {
            return Err("[ml_potential] GAP/QUIP supports delta = full or local".into());
        }
        let model = cfg
            .model
            .as_deref()
            .ok_or("[ml_potential] model is required for GAP/QUIP")?;

        let model_path = if Path::new(model).is_absolute() {
            Path::new(model).to_path_buf()
        } else {
            base_dir.join(model)
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

        let mut model = GapQuipModel {
            handle,
            weight: cfg.weight.unwrap_or(0.001),
            cutoff,
            delta,
            _species_names: species_names,
            local_builder: None,
            accepted_local_energies: None,
            pending_local: None,
            last_local_cluster_sizes: None,
        };
        if matches!(delta, MlEnergyDelta::Local) {
            let indices = (0..config.atoms.len()).collect::<Vec<_>>();
            model.accepted_local_energies = Some(model.per_atom_energies_or_panic(&indices));
            model.local_builder = Some(GapLocalClusterBuilder::new(config, cutoff)?);
        }
        Ok(model)
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

    fn per_atom_energies_or_panic(&mut self, indices: &[usize]) -> Vec<f64> {
        let mut energies = vec![0.0; indices.len()];
        let rc = unsafe {
            rsmith_gap_quip_per_atom_energies(
                self.handle,
                indices.len(),
                indices.as_ptr(),
                energies.as_mut_ptr(),
            )
        };
        if rc != 0 {
            panic!("GAP/QUIP failed to compute per-atom energies (code {rc})");
        }
        energies
    }

    fn local_delta_or_panic(
        &mut self,
        atom_idx: usize,
        old_pos: &[f64; 3],
        new_pos: &[f64; 3],
    ) -> f64 {
        assert!(
            self.pending_local.is_none(),
            "local GAP trial already pending"
        );
        let (central_atoms, cluster, cluster_atoms, cluster_positions, box_lengths) = {
            let builder = self
                .local_builder
                .as_ref()
                .expect("local GAP delta requested without a cluster builder");
            let new_position = wrap_position(*new_pos, builder.box_lengths);
            let (central_atoms, cluster) = builder
                .trial_cluster(atom_idx, *old_pos, new_position)
                .unwrap_or_else(|error| {
                    panic!("rsmith failed to build local GAP cluster: {error}")
                });
            let (cluster_atoms, cluster_positions) =
                builder.cluster_arrays(&cluster, atom_idx, new_position);
            (
                central_atoms,
                cluster,
                cluster_atoms,
                cluster_positions,
                builder.box_lengths,
            )
        };
        let species = cluster_atoms
            .iter()
            .map(|&atom| self._species_names[atom].as_ptr())
            .collect::<Vec<_>>();
        let rc = unsafe {
            rsmith_gap_quip_set_local_cluster(
                self.handle,
                cluster_atoms.len(),
                cluster.n_central,
                cluster_atoms.as_ptr(),
                species.as_ptr(),
                cluster_positions.as_ptr(),
                box_lengths.as_ptr(),
            )
        };
        if rc != 0 {
            panic!("GAP/QUIP failed to set local cluster (code {rc})");
        }
        let indices = (0..cluster.n_central).collect::<Vec<_>>();
        let new_local_energies = self.per_atom_energies_or_panic(&indices);
        let old_energy = central_atoms
            .iter()
            .map(|&atom| {
                self.accepted_local_energies
                    .as_ref()
                    .expect("local GAP energy cache is missing")[atom]
            })
            .sum::<f64>();
        let new_energy = new_local_energies.iter().sum::<f64>();
        self.last_local_cluster_sizes = Some((cluster.n_central, cluster_atoms.len()));
        self.pending_local = Some(PendingLocalDelta {
            central_atoms,
            local_energies: new_local_energies,
        });
        new_energy - old_energy
    }

    pub fn last_local_cluster_sizes(&self) -> Option<(usize, usize)> {
        self.last_local_cluster_sizes
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
        if let Some(energies) = &self.accepted_local_energies {
            energies.iter().sum()
        } else {
            let indices: Vec<usize> = (0..config.atoms.len()).collect();
            self.sum_per_atom_energy_or_panic(&indices)
        }
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
        match self.delta {
            MlEnergyDelta::Full => {
                let indices: Vec<usize> = (0..config.atoms.len()).collect();
                let old_energy = self.sum_per_atom_energy_or_panic(&indices);
                self.move_atom_or_panic(atom_idx, new_pos);
                let new_energy = self.sum_per_atom_energy_or_panic(&indices);
                new_energy - old_energy
            }
            MlEnergyDelta::Local => self.local_delta_or_panic(atom_idx, old_pos, new_pos),
            MlEnergyDelta::Incremental => unreachable!(),
        }
    }

    fn accept_move(&mut self, atom_idx: usize, new_pos: &[f64; 3]) {
        if matches!(self.delta, MlEnergyDelta::Local) {
            let pending = self
                .pending_local
                .take()
                .expect("accepted local GAP move without a pending trial");
            let cache = self
                .accepted_local_energies
                .as_mut()
                .expect("local GAP energy cache is missing");
            for (atom, energy) in pending
                .central_atoms
                .into_iter()
                .zip(pending.local_energies)
            {
                cache[atom] = energy;
            }
            self.local_builder
                .as_mut()
                .expect("local GAP cluster builder is missing")
                .accept_move(atom_idx, *new_pos);
        }
    }

    fn reject_move(&mut self, atom_idx: usize, old_pos: &[f64; 3]) {
        if matches!(self.delta, MlEnergyDelta::Local) {
            self.pending_local = None;
        } else {
            self.move_atom_or_panic(atom_idx, old_pos);
        }
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

fn wrap_position(mut position: [f64; 3], box_lengths: [f64; 3]) -> [f64; 3] {
    for axis in 0..3 {
        position[axis] = position[axis].rem_euclid(box_lengths[axis]);
    }
    position
}

fn minimum_image_displacement(
    first: [f64; 3],
    second: [f64; 3],
    box_lengths: [f64; 3],
) -> [f64; 3] {
    std::array::from_fn(|axis| {
        let displacement = second[axis] - first[axis];
        displacement - box_lengths[axis] * (displacement / box_lengths[axis]).round()
    })
}

fn squared_norm(vector: [f64; 3]) -> f64 {
    vector.iter().map(|component| component * component).sum()
}
