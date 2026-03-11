use std::collections::HashMap;
use std::fs;
use std::path::Path;

use serde::Deserialize;

use crate::constraints::{Constraints, CoordinationConstraint};
use crate::rmc::RmcParams;

/// Top-level TOML configuration.
#[derive(Debug, Deserialize)]
pub struct Config {
    pub system: SystemConfig,
    pub data: DataConfig,
    pub rmc: RmcConfig,
    pub sq: Option<SqConfig>,
    pub constraints: Option<ConstraintsConfig>,
}

#[derive(Debug, Deserialize)]
pub struct SystemConfig {
    pub structure: String,
    pub format: String,
    pub types: Option<HashMap<String, String>>,
}

#[derive(Debug, Deserialize)]
pub struct DataConfig {
    pub xray_sq: Option<DatasetConfig>,
    pub neutron_sq: Option<DatasetConfig>,
    pub xray_gr: Option<DatasetConfig>,
}

#[derive(Debug, Deserialize)]
pub struct DatasetConfig {
    pub file: String,
    pub weight: Option<f64>,
    pub sigma: Option<f64>,
    pub fit_min: Option<f64>,
    pub fit_max: Option<f64>,
    /// Q_max used when deriving g(r) from S(Q). Controls Q cutoff in inverse FT.
    pub qmax: Option<f64>,
    /// Apply Lorch modification in Q-space for g(r) inverse FT (default: true).
    pub lorch: Option<bool>,
}

#[derive(Debug, Deserialize)]
pub struct RmcConfig {
    pub max_moves: Option<u64>,
    pub max_step: Option<f64>,
    pub checkpoint_every: Option<u64>,
    pub seed: Option<u64>,
    pub print_every: Option<u64>,
    pub target_acceptance: Option<f64>,
    pub adjust_step_every: Option<u64>,
    pub anneal_start: Option<f64>,
    pub anneal_end: Option<f64>,
    pub anneal_steps: Option<u64>,
    pub convergence_threshold: Option<f64>,
    pub convergence_window: Option<u64>,
}

#[derive(Debug, Deserialize)]
pub struct SqConfig {
    pub qmin: Option<f64>,
    pub qmax: Option<f64>,
    pub nq: Option<usize>,
    pub lorch: Option<bool>,
    pub rdf_cutoff: Option<f64>,
    pub rdf_nbins: Option<usize>,
}

#[derive(Debug, Deserialize)]
pub struct ConstraintsConfig {
    pub min_distance: Option<HashMap<String, f64>>,
    pub coordination: Option<Vec<CoordinationConfig>>,
}

#[derive(Debug, Deserialize)]
pub struct CoordinationConfig {
    pub pair: String,
    pub min: usize,
    pub max: usize,
    pub cutoff: f64,
}

impl Config {
    pub fn load(path: &Path) -> Result<Self, Box<dyn std::error::Error>> {
        let content = fs::read_to_string(path)?;
        let config: Config = toml::from_str(&content)?;
        Ok(config)
    }

    /// Get LAMMPS type map as HashMap<u32, String>.
    pub fn type_map(&self) -> HashMap<u32, String> {
        match &self.system.types {
            Some(types) => types
                .iter()
                .map(|(k, v)| (k.parse::<u32>().unwrap(), v.clone()))
                .collect(),
            None => {
                // Default CaSiO3 mapping
                let mut map = HashMap::new();
                map.insert(1, "Ca".to_string());
                map.insert(2, "Si".to_string());
                map.insert(3, "O".to_string());
                map
            }
        }
    }

    /// Build RMC params from config.
    pub fn rmc_params(&self) -> RmcParams {
        let mut p = RmcParams::default();

        if let Some(v) = self.rmc.max_moves {
            p.max_moves = v;
        }
        if let Some(v) = self.rmc.max_step {
            p.max_step = v;
        }
        if let Some(v) = self.rmc.checkpoint_every {
            p.checkpoint_every = v;
        }
        if let Some(v) = self.rmc.seed {
            p.seed = v;
        }
        if let Some(v) = self.rmc.print_every {
            p.print_every = v;
        }
        if let Some(v) = self.rmc.target_acceptance {
            p.target_acceptance = v;
        }
        if let Some(v) = self.rmc.adjust_step_every {
            p.adjust_step_every = v;
        }
        if let Some(v) = self.rmc.anneal_start {
            p.anneal_start = v;
        }
        if let Some(v) = self.rmc.anneal_end {
            p.anneal_end = v;
        }
        if let Some(v) = self.rmc.anneal_steps {
            p.anneal_steps = v;
        }
        if let Some(v) = self.rmc.convergence_threshold {
            p.convergence_threshold = v;
        }
        if let Some(v) = self.rmc.convergence_window {
            p.convergence_window = v;
        }

        if let Some(ref sq) = self.sq {
            let qmin = sq.qmin.unwrap_or(0.3);
            let qmax = sq.qmax.unwrap_or(20.0);
            let nq = sq.nq.unwrap_or(500);
            let dq = (qmax - qmin) / nq as f64;
            p.q_grid = (0..nq).map(|i| qmin + i as f64 * dq).collect();
            p.lorch = sq.lorch.unwrap_or(true);
            if let Some(v) = sq.rdf_cutoff {
                p.rdf_cutoff = v;
            }
            if let Some(v) = sq.rdf_nbins {
                p.rdf_nbins = v;
            }
        }

        p
    }

    /// Build constraints from config.
    pub fn constraints(&self) -> Constraints {
        let mut c = Constraints::new();

        if let Some(ref cc) = self.constraints {
            if let Some(ref md) = cc.min_distance {
                c.min_distances = md.clone();
            }
            if let Some(ref coords) = cc.coordination {
                for coord in coords {
                    c.coordination.push(CoordinationConstraint {
                        pair: coord.pair.clone(),
                        min: coord.min,
                        max: coord.max,
                        cutoff: coord.cutoff,
                    });
                }
            }
        }

        c
    }
}
