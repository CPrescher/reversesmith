//! MACE energy backend through a small Python worker.
//!
//! The default backend evaluates full-system energies for correctness. For
//! ordinary short-range MACE models, `delta = "local"` builds a bounded
//! explicit-image cluster with a Rust cell list and asks the worker to sum the
//! affected per-atom energies.

use std::collections::HashSet;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};

use serde_json::{json, Value};

use crate::atoms::Configuration;
use crate::cells::CellList;
use crate::config::{MlEnergyDelta, MlPotentialConfig};
use crate::energy::EnergyModel;

const MACE_WORKER: &str = include_str!("../../python/rsmith_mace_worker.py");

pub struct MacePythonModel {
    child: Child,
    stdin: std::io::BufWriter<ChildStdin>,
    stdout: BufReader<ChildStdout>,
    weight: f64,
    cutoff: f64,
    delta: MlEnergyDelta,
    local_cluster_builder: Option<MaceLocalClusterBuilder>,
    last_local_cluster_sizes: Option<(usize, usize)>,
}

#[derive(Debug, Clone, serde::Serialize)]
struct MaceClusterSpec {
    atoms: Vec<usize>,
    shifts: Vec<[i32; 3]>,
    central_indices: Vec<usize>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct PeriodicImage {
    atom: usize,
    shift: [i32; 3],
}

/// Builds the explicit periodic-image clusters needed by local MACE without
/// rescanning atom/central/image triples in Python.
struct MaceLocalClusterBuilder {
    positions: Vec<[f64; 3]>,
    box_lengths: [f64; 3],
    radius: f64,
    cell_list: MaceCellList,
}

impl MaceLocalClusterBuilder {
    fn new(configuration: &Configuration, radius: f64) -> Result<Self, String> {
        if !radius.is_finite() || radius <= 0.0 {
            return Err("local MACE radius must be finite and positive".to_string());
        }
        if configuration.atoms.is_empty() {
            return Err("local MACE requires at least one atom".to_string());
        }
        if configuration
            .box_lengths
            .iter()
            .any(|length| !length.is_finite() || *length <= 0.0)
        {
            return Err("local MACE box lengths must be finite and positive".to_string());
        }
        let positions = configuration
            .atoms
            .iter()
            .map(|atom| wrap_position(atom.position, configuration.box_lengths))
            .collect::<Vec<_>>();
        let cell_list = MaceCellList::new(&positions, configuration.box_lengths, radius);
        Ok(Self {
            positions,
            box_lengths: configuration.box_lengths,
            radius,
            cell_list,
        })
    }

    fn trial_cluster(
        &self,
        moved_atom: usize,
        old_position: [f64; 3],
        new_position: [f64; 3],
    ) -> Result<(Vec<usize>, MaceClusterSpec), String> {
        if moved_atom >= self.positions.len() {
            return Err(format!(
                "local MACE atom index {moved_atom} is out of range"
            ));
        }
        if old_position
            .iter()
            .chain(new_position.iter())
            .any(|coordinate| !coordinate.is_finite())
        {
            return Err("local MACE trial positions must be finite".to_string());
        }
        let old_position = wrap_position(old_position, self.box_lengths);
        let new_position = wrap_position(new_position, self.box_lengths);
        let accepted_error = squared_norm(minimum_image_displacement(
            self.positions[moved_atom],
            old_position,
            self.box_lengths,
        ));
        if accepted_error > 1.0e-20 {
            return Err("local MACE old position does not match its accepted state".to_string());
        }

        let central_atoms = self.affected_atoms(moved_atom, old_position, new_position);
        let after_cluster = self.build_cluster(&central_atoms, moved_atom, new_position);
        Ok((central_atoms, after_cluster))
    }

    fn affected_atoms(
        &self,
        moved_atom: usize,
        old_position: [f64; 3],
        new_position: [f64; 3],
    ) -> Vec<usize> {
        let radius_squared = self.radius * self.radius;
        self.positions
            .iter()
            .enumerate()
            .filter_map(|(atom, position)| {
                if atom == moved_atom
                    || squared_norm(minimum_image_displacement(
                        *position,
                        old_position,
                        self.box_lengths,
                    )) <= radius_squared
                    || squared_norm(minimum_image_displacement(
                        *position,
                        new_position,
                        self.box_lengths,
                    )) <= radius_squared
                {
                    Some(atom)
                } else {
                    None
                }
            })
            .collect()
    }

    fn build_cluster(
        &self,
        central_atoms: &[usize],
        moved_atom: usize,
        moved_position: [f64; 3],
    ) -> MaceClusterSpec {
        let mut images = HashSet::new();
        let needs_multiple_images = self
            .box_lengths
            .iter()
            .any(|length| 2.0 * self.radius >= *length);

        let central_images = central_atoms
            .iter()
            .map(|&central_atom| {
                let atom_position = if central_atom == moved_atom {
                    moved_position
                } else {
                    self.positions[central_atom]
                };
                PeriodicImage {
                    atom: central_atom,
                    shift: image_shift_near(atom_position, moved_position, self.box_lengths),
                }
            })
            .collect::<Vec<_>>();

        for central_image in &central_images {
            let central_atom = central_image.atom;
            let central_base_position = if central_atom == moved_atom {
                moved_position
            } else {
                self.positions[central_atom]
            };
            let central_position = std::array::from_fn(|axis| {
                central_base_position[axis]
                    + central_image.shift[axis] as f64 * self.box_lengths[axis]
            });
            if needs_multiple_images {
                for atom in 0..self.positions.len() {
                    let atom_position = if atom == moved_atom {
                        moved_position
                    } else {
                        self.positions[atom]
                    };
                    self.insert_all_images_within(
                        &mut images,
                        atom,
                        atom_position,
                        central_position,
                    );
                }
            } else {
                for atom in self.cell_list.candidate_atoms(central_position) {
                    if atom == moved_atom {
                        continue;
                    }
                    self.insert_minimum_image(
                        &mut images,
                        atom,
                        self.positions[atom],
                        central_position,
                    );
                }
                self.insert_minimum_image(
                    &mut images,
                    moved_atom,
                    moved_position,
                    central_position,
                );
            }
        }

        let mut images = images.into_iter().collect::<Vec<_>>();
        images.sort_unstable();
        debug_assert!(central_images
            .iter()
            .all(|central| images.contains(central)));
        let central_indices = central_images
            .iter()
            .map(|central| {
                images
                    .binary_search(central)
                    .expect("local MACE cluster missed a central periodic image")
            })
            .collect();
        MaceClusterSpec {
            atoms: images.iter().map(|image| image.atom).collect(),
            shifts: images.iter().map(|image| image.shift).collect(),
            central_indices,
        }
    }

    fn insert_minimum_image(
        &self,
        images: &mut HashSet<PeriodicImage>,
        atom: usize,
        atom_position: [f64; 3],
        central_position: [f64; 3],
    ) {
        let mut displacement = [0.0; 3];
        let mut shift = [0; 3];
        for axis in 0..3 {
            let raw = atom_position[axis] - central_position[axis];
            shift[axis] = -(raw / self.box_lengths[axis]).round() as i32;
            displacement[axis] = raw + shift[axis] as f64 * self.box_lengths[axis];
        }
        if squared_norm(displacement) <= self.radius * self.radius + 1.0e-12 {
            images.insert(PeriodicImage { atom, shift });
        }
    }

    fn insert_all_images_within(
        &self,
        images: &mut HashSet<PeriodicImage>,
        atom: usize,
        atom_position: [f64; 3],
        central_position: [f64; 3],
    ) {
        let minimum_shift: [i32; 3] = std::array::from_fn(|axis| {
            ((central_position[axis] - atom_position[axis] - self.radius) / self.box_lengths[axis])
                .ceil() as i32
        });
        let maximum_shift: [i32; 3] = std::array::from_fn(|axis| {
            ((central_position[axis] - atom_position[axis] + self.radius) / self.box_lengths[axis])
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
                    if squared_norm(displacement) <= self.radius * self.radius + 1.0e-12 {
                        images.insert(PeriodicImage { atom, shift });
                    }
                }
            }
        }
    }

    fn accept_move(&mut self, atom: usize, new_position: [f64; 3]) {
        let new_position = wrap_position(new_position, self.box_lengths);
        self.positions[atom] = new_position;
        self.cell_list.move_atom(atom, new_position);
    }
}

struct MaceCellList {
    cell_counts: [usize; 3],
    box_lengths: [f64; 3],
    buckets: Vec<Vec<usize>>,
    atom_cells: Vec<usize>,
}

impl MaceCellList {
    fn new(positions: &[[f64; 3]], box_lengths: [f64; 3], radius: f64) -> Self {
        let cell_counts =
            std::array::from_fn(|axis| (box_lengths[axis] / radius).floor().max(1.0) as usize);
        let bucket_count = cell_counts.iter().product();
        let mut list = Self {
            cell_counts,
            box_lengths,
            buckets: vec![Vec::new(); bucket_count],
            atom_cells: vec![0; positions.len()],
        };
        for (atom, position) in positions.iter().enumerate() {
            let cell = list.cell_for_position(*position);
            list.atom_cells[atom] = cell;
            list.buckets[cell].push(atom);
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

    fn move_atom(&mut self, atom: usize, new_position: [f64; 3]) {
        let old_cell = self.atom_cells[atom];
        let new_cell = self.cell_for_position(new_position);
        if old_cell == new_cell {
            return;
        }
        let offset = self.buckets[old_cell]
            .iter()
            .position(|candidate| *candidate == atom)
            .expect("local MACE cell list lost an atom");
        self.buckets[old_cell].swap_remove(offset);
        self.buckets[new_cell].push(atom);
        self.atom_cells[atom] = new_cell;
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
    first: [f64; 3],
    second: [f64; 3],
    box_lengths: [f64; 3],
) -> [f64; 3] {
    std::array::from_fn(|axis| {
        let displacement = second[axis] - first[axis];
        displacement - box_lengths[axis] * (displacement / box_lengths[axis]).round()
    })
}

fn image_shift_near(position: [f64; 3], reference: [f64; 3], box_lengths: [f64; 3]) -> [i32; 3] {
    std::array::from_fn(|axis| {
        -((position[axis] - reference[axis]) / box_lengths[axis]).round() as i32
    })
}

fn squared_norm(vector: [f64; 3]) -> f64 {
    vector.iter().map(|component| component * component).sum()
}

impl MacePythonModel {
    pub fn from_config(
        cfg: &MlPotentialConfig,
        config: &Configuration,
        base_dir: &Path,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let cutoff = cfg
            .cutoff
            .ok_or("[ml_potential] cutoff is required for MACE/Python")?;
        if cutoff <= 0.0 {
            return Err("[ml_potential] cutoff must be greater than 0".into());
        }
        if cfg.torch_threads == Some(0) {
            return Err("[ml_potential] torch_threads must be greater than 0".into());
        }

        let model = cfg
            .model
            .as_deref()
            .ok_or("[ml_potential] model is required for MACE/Python")?;
        let model_path = resolve_path(base_dir, model);
        if !model_path.exists() {
            return Err(format!("MACE model file not found: {}", model_path.display()).into());
        }

        let python = cfg.python.as_deref().unwrap_or("python3");
        let worker_path = cfg
            .worker
            .as_ref()
            .map(|worker| resolve_path(base_dir, worker));
        if let Some(worker_path) = &worker_path {
            if !worker_path.exists() {
                return Err(format!(
                    "MACE Python worker file not found: {}",
                    worker_path.display()
                )
                .into());
            }
        }

        let mut child = spawn_worker(python, worker_path.as_deref())?;
        let child_stdin = child
            .stdin
            .take()
            .ok_or("failed to open MACE worker stdin")?;
        let child_stdout = child
            .stdout
            .take()
            .ok_or("failed to open MACE worker stdout")?;
        let mut model = MacePythonModel {
            child,
            stdin: std::io::BufWriter::new(child_stdin),
            stdout: BufReader::new(child_stdout),
            weight: cfg.weight.unwrap_or(0.001),
            cutoff,
            delta: cfg.delta.unwrap_or(MlEnergyDelta::Full),
            local_cluster_builder: None,
            last_local_cluster_sizes: None,
        };

        let init = json!({
            "cmd": "init",
            "model": model_path.to_string_lossy(),
            "device": cfg.device.as_deref().unwrap_or("cpu"),
            "torch_threads": cfg.torch_threads,
            "delta": delta_label(cfg.delta.unwrap_or(MlEnergyDelta::Full)),
            "species": config.atoms.iter().map(|atom| atom.species.as_str()).collect::<Vec<_>>(),
            "positions": config.atoms.iter().map(|atom| wrap_position(atom.position, config.box_lengths)).collect::<Vec<_>>(),
            "box": config.box_lengths,
        });
        model.request(init)?;

        if matches!(model.delta, MlEnergyDelta::Local) {
            let metadata = model.request(json!({"cmd": "metadata"}))?;
            let local_radius = metadata
                .get("local_context_radius")
                .and_then(Value::as_f64)
                .ok_or("MACE Python worker metadata did not contain local_context_radius")?;
            model.local_cluster_builder = Some(MaceLocalClusterBuilder::new(config, local_radius)?);
        }

        Ok(model)
    }

    fn request(&mut self, request: Value) -> Result<Value, Box<dyn std::error::Error>> {
        serde_json::to_writer(&mut self.stdin, &request)?;
        self.stdin.write_all(b"\n")?;
        self.stdin.flush()?;

        let mut line = String::new();
        let n = self.stdout.read_line(&mut line)?;
        if n == 0 {
            return Err("MACE Python worker exited without a response".into());
        }

        let response: Value = serde_json::from_str(&line)?;
        if response.get("ok").and_then(Value::as_bool) == Some(true) {
            Ok(response)
        } else {
            let error = response
                .get("error")
                .and_then(Value::as_str)
                .unwrap_or("unknown MACE Python worker error");
            Err(error.to_string().into())
        }
    }

    fn total_energy_or_panic(&mut self) -> f64 {
        let response = self
            .request(json!({"cmd": "energy"}))
            .unwrap_or_else(|e| panic!("MACE Python energy evaluation failed: {e}"));
        response
            .get("energy")
            .and_then(Value::as_f64)
            .unwrap_or_else(|| panic!("MACE Python worker response did not contain numeric energy"))
    }

    fn move_atom_or_panic(&mut self, atom_idx: usize, pos: &[f64; 3]) {
        self.request(json!({
            "cmd": "move",
            "atom": atom_idx,
            "position": pos,
        }))
        .unwrap_or_else(|e| panic!("MACE Python failed to move atom {atom_idx}: {e}"));
    }

    fn local_delta_or_panic(
        &mut self,
        atom_idx: usize,
        old_pos: &[f64; 3],
        new_pos: &[f64; 3],
    ) -> f64 {
        let (old_position, new_position, central_atoms, after_cluster) = {
            let builder = self
                .local_cluster_builder
                .as_ref()
                .expect("local MACE delta requested without a cluster builder");
            let old_position = wrap_position(*old_pos, builder.box_lengths);
            let new_position = wrap_position(*new_pos, builder.box_lengths);
            let (central_atoms, after_cluster) = builder
                .trial_cluster(atom_idx, old_position, new_position)
                .unwrap_or_else(|error| {
                    panic!("rsmith failed to build local MACE trial clusters: {error}")
                });
            (old_position, new_position, central_atoms, after_cluster)
        };
        let response = self
            .request(json!({
                "cmd": "local_delta",
                "atom": atom_idx,
                "old_position": old_position,
                "new_position": new_position,
                "central_atoms": central_atoms,
                "after_cluster": after_cluster,
            }))
            .unwrap_or_else(|e| panic!("MACE Python local delta failed: {e}"));
        self.last_local_cluster_sizes = Some((
            response
                .get("central_atoms")
                .and_then(Value::as_u64)
                .unwrap_or(0) as usize,
            response
                .get("cluster_atoms")
                .and_then(Value::as_u64)
                .unwrap_or(0) as usize,
        ));
        response
            .get("delta")
            .and_then(Value::as_f64)
            .unwrap_or_else(|| panic!("MACE Python worker response did not contain numeric delta"))
    }

    /// Diagnostics from the most recent local trial: affected central atoms
    /// and atoms (including periodic images) in the evaluated cluster.
    pub fn last_local_cluster_sizes(&self) -> Option<(usize, usize)> {
        self.last_local_cluster_sizes
    }
}

impl EnergyModel for MacePythonModel {
    fn label(&self) -> &str {
        "MACE/Python"
    }

    fn weight(&self) -> f64 {
        self.weight
    }

    fn cutoff(&self) -> f64 {
        self.cutoff
    }

    fn total_energy(&mut self, _config: &Configuration, _cell_list: &CellList) -> f64 {
        self.total_energy_or_panic()
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
        match self.delta {
            MlEnergyDelta::Full => {
                let old_energy = self.total_energy_or_panic();
                self.move_atom_or_panic(atom_idx, new_pos);
                let new_energy = self.total_energy_or_panic();
                new_energy - old_energy
            }
            MlEnergyDelta::Local => self.local_delta_or_panic(atom_idx, old_pos, new_pos),
        }
    }

    fn accept_move(&mut self, atom_idx: usize, new_pos: &[f64; 3]) {
        // The worker state was already moved during energy_delta_atom.
        if matches!(self.delta, MlEnergyDelta::Local) {
            self.request(json!({"cmd": "accept_local"}))
                .unwrap_or_else(|error| panic!("MACE Python local accept failed: {error}"));
        }
        if let Some(builder) = &mut self.local_cluster_builder {
            builder.accept_move(atom_idx, *new_pos);
        }
    }

    fn reject_move(&mut self, atom_idx: usize, old_pos: &[f64; 3]) {
        match self.delta {
            MlEnergyDelta::Full => self.move_atom_or_panic(atom_idx, old_pos),
            MlEnergyDelta::Local => {
                let old_position = self
                    .local_cluster_builder
                    .as_ref()
                    .map(|builder| wrap_position(*old_pos, builder.box_lengths))
                    .unwrap_or(*old_pos);
                self.request(json!({
                    "cmd": "reject_local",
                    "atom": atom_idx,
                    "old_position": old_position,
                }))
                .unwrap_or_else(|error| panic!("MACE Python local reject failed: {error}"));
            }
        }
    }
}

impl Drop for MacePythonModel {
    fn drop(&mut self) {
        let _ = self.request(json!({"cmd": "shutdown"}));
        let _ = self.child.wait();
    }
}

fn spawn_worker(
    python: &str,
    worker_path: Option<&Path>,
) -> Result<Child, Box<dyn std::error::Error>> {
    let mut cmd = Command::new(python);
    cmd.arg("-u");
    if let Some(worker_path) = worker_path {
        cmd.arg(worker_path);
    } else {
        cmd.arg("-c").arg(MACE_WORKER);
    }
    Ok(cmd
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit())
        .spawn()?)
}

fn resolve_path(base_dir: &Path, path: &str) -> PathBuf {
    let path = Path::new(path);
    if path.is_absolute() {
        path.to_path_buf()
    } else {
        base_dir.join(path)
    }
}

fn delta_label(delta: MlEnergyDelta) -> &'static str {
    match delta {
        MlEnergyDelta::Full => "full",
        MlEnergyDelta::Local => "local",
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::fs;
    use std::process::Command;

    use crate::atoms::{Atom, Configuration};
    use crate::cells::CellList;
    use crate::config::{MlBackend, MlEnergyDelta, MlPotentialConfig};
    use crate::energy::EnergyModel;

    use super::{MaceLocalClusterBuilder, MacePythonModel, PeriodicImage};

    fn python_available() -> bool {
        Command::new("python3")
            .arg("--version")
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status()
            .is_ok_and(|status| status.success())
    }

    fn small_config() -> Configuration {
        let atoms = vec![
            Atom {
                position: [1.0, 0.0, 0.0],
                species: "Si".to_string(),
                type_id: 0,
            },
            Atom {
                position: [0.0, 2.0, 0.0],
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

    fn silicon_config(positions: Vec<[f64; 3]>, box_lengths: [f64; 3]) -> Configuration {
        let atom_count = positions.len();
        Configuration {
            atoms: positions
                .into_iter()
                .map(|position| Atom {
                    position,
                    species: "Si".to_string(),
                    type_id: 0,
                })
                .collect(),
            box_lengths,
            species: vec!["Si".to_string()],
            composition: HashMap::from([("Si".to_string(), atom_count)]),
        }
    }

    fn cluster_images(
        cluster: &super::MaceClusterSpec,
    ) -> std::collections::HashSet<PeriodicImage> {
        cluster
            .atoms
            .iter()
            .copied()
            .zip(cluster.shifts.iter().copied())
            .map(|(atom, shift)| PeriodicImage { atom, shift })
            .collect()
    }

    #[test]
    fn rust_local_cluster_tracks_periodic_image_shifts() {
        let configuration = silicon_config(
            vec![[0.5, 0.5, 0.5], [10.0, 10.0, 10.0], [19.5, 0.5, 0.5]],
            [20.0; 3],
        );
        let builder = MaceLocalClusterBuilder::new(&configuration, 2.0).unwrap();
        let (central, cluster) = builder
            .trial_cluster(0, [0.5, 0.5, 0.5], [0.6, 0.5, 0.5])
            .unwrap();
        assert_eq!(central, vec![0, 2]);
        let images = cluster_images(&cluster);
        assert!(images.contains(&PeriodicImage {
            atom: 0,
            shift: [0, 0, 0]
        }));
        assert!(images.contains(&PeriodicImage {
            atom: 2,
            shift: [-1, 0, 0]
        }));
        assert!(!images.iter().any(|image| image.atom == 1));
    }

    #[test]
    fn rust_local_cluster_keeps_multiple_images_in_small_boxes() {
        let configuration = silicon_config(vec![[0.0, 0.0, 0.0]], [5.0; 3]);
        let builder = MaceLocalClusterBuilder::new(&configuration, 6.0).unwrap();
        let (_, cluster) = builder
            .trial_cluster(0, [0.0, 0.0, 0.0], [0.1, 0.0, 0.0])
            .unwrap();
        let images = cluster_images(&cluster);
        assert!(images.contains(&PeriodicImage {
            atom: 0,
            shift: [-1, 0, 0]
        }));
        assert!(images.contains(&PeriodicImage {
            atom: 0,
            shift: [0, 0, 0]
        }));
        assert!(images.contains(&PeriodicImage {
            atom: 0,
            shift: [1, 0, 0]
        }));
        assert!(!images.contains(&PeriodicImage {
            atom: 0,
            shift: [1, 1, 0]
        }));
    }

    #[test]
    fn rust_local_cluster_updates_accepted_positions() {
        let configuration = silicon_config(vec![[1.0, 1.0, 1.0]], [20.0; 3]);
        let mut builder = MaceLocalClusterBuilder::new(&configuration, 3.0).unwrap();
        builder.accept_move(0, [4.0, 1.0, 1.0]);
        assert!(builder
            .trial_cluster(0, [4.0, 1.0, 1.0], [4.1, 1.0, 1.0])
            .is_ok());
        assert!(builder
            .trial_cluster(0, [1.0, 1.0, 1.0], [1.1, 1.0, 1.0])
            .is_err());
    }

    #[test]
    fn mace_python_backend_updates_and_reverts_worker_state() {
        if !python_available() {
            eprintln!("skipping MACE Python mock test because python3 is unavailable");
            return;
        }

        let temp_dir =
            std::env::temp_dir().join(format!("rsmith-mace-python-test-{}", std::process::id()));
        fs::create_dir_all(&temp_dir).unwrap();
        let worker_path = temp_dir.join("mock_worker.py");
        let model_path = temp_dir.join("mock.model");
        fs::write(&model_path, "mock").unwrap();
        fs::write(
            &worker_path,
            r#"
import json
import sys

positions = []

def energy():
    return sum(sum(x * x for x in pos) for pos in positions)

for line in sys.stdin:
    request = json.loads(line)
    cmd = request["cmd"]
    if cmd == "init":
        positions = [list(pos) for pos in request["positions"]]
        response = {"ok": True}
    elif cmd == "metadata":
        response = {
            "ok": True,
            "r_max": 2.5,
            "num_interactions": 2,
            "local_central_radius": 5.0,
            "local_context_radius": 5.0,
        }
    elif cmd == "energy":
        response = {"ok": True, "energy": energy()}
    elif cmd == "move":
        positions[int(request["atom"])] = list(request["position"])
        response = {"ok": True}
    elif cmd == "local_delta":
        old = energy()
        positions[int(request["atom"])] = list(request["new_position"])
        response = {"ok": True, "delta": energy() - old}
    elif cmd == "accept_local":
        response = {"ok": True}
    elif cmd == "reject_local":
        positions[int(request["atom"])] = list(request["old_position"])
        response = {"ok": True}
    elif cmd == "shutdown":
        print(json.dumps({"ok": True}), flush=True)
        break
    else:
        response = {"ok": False, "error": "unknown command"}
    print(json.dumps(response), flush=True)
"#,
        )
        .unwrap();

        let config = small_config();
        let cell_list = CellList::new(
            &config
                .atoms
                .iter()
                .map(|atom| atom.position)
                .collect::<Vec<_>>(),
            &config.box_lengths,
            4.0,
        );
        let cfg = MlPotentialConfig {
            backend: MlBackend::MacePython,
            model: Some(model_path.to_string_lossy().to_string()),
            coefficient_file: None,
            parameter_file: None,
            init_args: None,
            weight: Some(0.5),
            cutoff: Some(5.0),
            delta: None,
            device: Some("cpu".to_string()),
            torch_threads: Some(1),
            python: Some("python3".to_string()),
            worker: Some(worker_path.to_string_lossy().to_string()),
        };
        let mut model = MacePythonModel::from_config(&cfg, &config, &temp_dir).unwrap();

        assert_eq!(model.total_energy(&config, &cell_list), 5.0);
        let delta = model.energy_delta_atom(
            &config,
            0,
            &[1.0, 0.0, 0.0],
            &[2.0, 0.0, 0.0],
            &cell_list,
            0,
            0,
        );
        assert_eq!(delta, 3.0);
        assert_eq!(model.total_energy(&config, &cell_list), 8.0);

        model.reject_move(0, &[1.0, 0.0, 0.0]);
        assert_eq!(model.total_energy(&config, &cell_list), 5.0);
        drop(model);

        let mut local_cfg = cfg.clone();
        local_cfg.delta = Some(MlEnergyDelta::Local);
        let mut local_model = MacePythonModel::from_config(&local_cfg, &config, &temp_dir).unwrap();
        let local_delta = local_model.energy_delta_atom(
            &config,
            0,
            &[1.0, 0.0, 0.0],
            &[2.0, 0.0, 0.0],
            &cell_list,
            0,
            0,
        );
        assert_eq!(local_delta, 3.0);
        assert_eq!(local_model.total_energy(&config, &cell_list), 8.0);
        local_model.reject_move(0, &[1.0, 0.0, 0.0]);
        assert_eq!(local_model.total_energy(&config, &cell_list), 5.0);
    }
}
