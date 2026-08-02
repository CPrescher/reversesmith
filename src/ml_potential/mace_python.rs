//! MACE energy backend through a small Python worker.
//!
//! The default backend evaluates full-system energies for correctness. For
//! ordinary short-range MACE models, `delta = "local"` asks the worker to build
//! a bounded explicit-image cluster around the moved atom and sum the affected
//! per-atom energies.

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
        };

        let init = json!({
            "cmd": "init",
            "model": model_path.to_string_lossy(),
            "device": cfg.device.as_deref().unwrap_or("cpu"),
            "torch_threads": cfg.torch_threads,
            "delta": delta_label(cfg.delta.unwrap_or(MlEnergyDelta::Full)),
            "species": config.atoms.iter().map(|atom| atom.species.as_str()).collect::<Vec<_>>(),
            "positions": config.atoms.iter().map(|atom| atom.position).collect::<Vec<_>>(),
            "box": config.box_lengths,
        });
        model.request(init)?;

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
        let response = self
            .request(json!({
                "cmd": "local_delta",
                "atom": atom_idx,
                "old_position": old_pos,
                "new_position": new_pos,
            }))
            .unwrap_or_else(|e| panic!("MACE Python local delta failed: {e}"));
        response
            .get("delta")
            .and_then(Value::as_f64)
            .unwrap_or_else(|| panic!("MACE Python worker response did not contain numeric delta"))
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

    fn accept_move(&mut self, _atom_idx: usize, _new_pos: &[f64; 3]) {
        // The worker state was already moved during energy_delta_atom.
    }

    fn reject_move(&mut self, atom_idx: usize, old_pos: &[f64; 3]) {
        self.move_atom_or_panic(atom_idx, old_pos);
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

    use super::MacePythonModel;

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
    elif cmd == "energy":
        response = {"ok": True, "energy": energy()}
    elif cmd == "move":
        positions[int(request["atom"])] = list(request["position"])
        response = {"ok": True}
    elif cmd == "local_delta":
        old = energy()
        positions[int(request["atom"])] = list(request["new_position"])
        response = {"ok": True, "delta": energy() - old}
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
