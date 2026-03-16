use std::collections::HashMap;
use std::fs;
use std::path::Path;

use serde::Deserialize;

use crate::constraints::{Constraints, CoordinationConstraint};
use crate::rmc::RmcParams;

/// Top-level TOML configuration.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    pub system: SystemConfig,
    pub data: DataConfig,
    pub rmc: RmcConfig,
    pub sq: Option<SqConfig>,
    pub constraints: Option<ConstraintsConfig>,
    pub analysis: Option<AnalysisConfig>,
    pub potential: Option<PotentialConfig>,
    pub epsr: Option<EpsrConfig>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SystemConfig {
    pub structure: String,
    pub format: String,
    pub types: Option<HashMap<String, String>>,
    /// Target mass density in g/cm^3.  When set, the box and atom positions
    /// are rescaled isotropically before any computation begins.
    pub density: Option<f64>,
    /// Write a VASP POSCAR file alongside the refined XYZ output (default: false).
    pub output_poscar: Option<bool>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DataConfig {
    pub xray_sq: Option<DatasetConfig>,
    pub neutron_sq: Option<DatasetConfig>,
    pub xray_gr: Option<DatasetConfig>,
    pub xray_fr: Option<DatasetConfig>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DatasetConfig {
    pub file: String,
    pub weight: Option<f64>,
    pub sigma: Option<f64>,
    /// Linear Q-dependent scaling for sigma: sigma(Q) *= 1 + sigma_alpha * Q.
    /// At Q=0 sigma is unchanged; at Q=20 with sigma_alpha=0.05 sigma doubles.
    /// Default: 0.0 (no Q-scaling). Only applies to S(Q) datasets, ignored for g(r).
    pub sigma_alpha: Option<f64>,
    pub fit_min: Option<f64>,
    pub fit_max: Option<f64>,
    /// Q_max used when deriving g(r) from S(Q). Controls Q cutoff in inverse FT.
    pub qmax: Option<f64>,
    /// Apply Lorch modification in Q-space for g(r) inverse FT (default: true).
    pub lorch: Option<bool>,
    /// Data convention: "sq" for S(Q) (Faber-Ziman, oscillates around 1),
    /// "iq" for i(Q) = S(Q) - 1 (interference function, oscillates around 0),
    /// "fq" for F(Q) = Q*(S(Q) - 1) (reduced interference function).
    /// Default: "sq".
    pub convention: Option<String>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
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
#[serde(deny_unknown_fields)]
pub struct SqConfig {
    pub qmin: Option<f64>,
    pub qmax: Option<f64>,
    pub nq: Option<usize>,
    pub lorch: Option<bool>,
    pub rdf_cutoff: Option<f64>,
    pub rdf_nbins: Option<usize>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ConstraintsConfig {
    pub min_distance: Option<HashMap<String, f64>>,
    pub coordination: Option<Vec<CoordinationConfig>>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CoordinationConfig {
    pub pair: String,
    pub min: usize,
    pub max: usize,
    pub cutoff: f64,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AnalysisConfig {
    pub cutoffs: Option<HashMap<String, f64>>,
    pub angle_triplets: Option<Vec<String>>,
    pub angle_bins: Option<usize>,
}

/// EPSR (Empirical Potential Structure Refinement) configuration.
///
/// When present, wraps the RMC inner loop in an outer EPSR iteration that
/// refines a perturbation potential so the simulation naturally reproduces
/// the experimental data (Soper, 1996).
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EpsrConfig {
    /// EPSR mode: "hybrid" (chi2 + energy MC, default) or "pure" (energy-only MC
    /// with S(Q) computed after each epoch — the original Soper algorithm).
    pub mode: Option<String>,
    /// Number of outer EPSR iterations (default: 10).
    pub iterations: Option<usize>,
    /// Feedback factor for EP update: EP += feedback * kT * Δg(r) (default: 0.2).
    pub feedback: Option<f64>,
    /// Gaussian smoothing width in Å for EP (default: 0.02).
    pub smooth_sigma: Option<f64>,
    /// MC moves per EPSR epoch. Overrides [rmc] max_moves during EPSR (default: from [rmc]).
    pub moves_per_iteration: Option<u64>,
    /// kT in eV for the EP update (default: 0.025, ~300K).
    pub temperature: Option<f64>,
    /// Zero EP below this distance in Å (default: 1.0).
    pub min_r: Option<f64>,
    /// Stop if relative EP change (max |ΔEP| / max |EP|) falls below this for
    /// `convergence_window` consecutive iterations (default: 0.0 = run all iterations).
    pub convergence: Option<f64>,
    /// Number of consecutive iterations that must satisfy the convergence criterion
    /// before stopping (default: 3).
    pub convergence_window: Option<usize>,
    /// Directory containing previous `epsr_ep_{pair}.dat` files to seed the EP.
    /// When set, EP tables are loaded and accumulated on top (restart from previous run).
    pub ep_restart: Option<String>,
}

/// Configuration for pair potentials (hybrid RMC).
///
/// Parsed from the `[potential]` TOML section. Multiple potential types
/// can be combined: analytical forms are summed, tabulated replaces all.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PotentialConfig {
    pub weight: Option<f64>,
    pub cutoff: Option<f64>,
    pub buckingham: Option<Vec<BuckinghamConfig>>,
    pub pedone: Option<Vec<PedoneConfig>>,
    pub coulomb: Option<CoulombConfig>,
    pub tabulated: Option<Vec<TabulatedConfig>>,
}

/// Buckingham potential: V(r) = A exp(-r/rho) - C/r^6.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BuckinghamConfig {
    pub pair: String,
    #[serde(rename = "A")]
    pub a_param: f64,
    pub rho: f64,
    #[serde(rename = "C")]
    pub c_param: f64,
}

/// Pedone potential: V(r) = D0 [1 - exp(-alpha (r-r0))]^2 - D0 + C0/r^12.
///
/// Pedone et al. (2006), J. Phys. Chem. B, 110, 11780.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PedoneConfig {
    pub pair: String,
    #[serde(rename = "D0")]
    pub d0: f64,
    pub alpha: f64,
    pub r0: f64,
    #[serde(rename = "C0")]
    pub c0: f64,
}

/// Coulomb DSF electrostatics, matching LAMMPS `pair_style coul/dsf`.
///
/// V(r) = K qi qj [erfc(alpha r)/r - erfc(alpha rc)/rc + correction(r - rc)]
///
/// Automatically generates pair potentials for all species combinations
/// with nonzero charge products.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CoulombConfig {
    /// DSF damping parameter (1/A), matching LAMMPS coul/dsf alpha
    pub alpha: f64,
    /// Per-species charges, e.g. { "Ca" = 1.2, "Si" = 2.4, "O" = -1.2 }
    pub charges: HashMap<String, f64>,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TabulatedConfig {
    pub pair: String,
    pub file: String,
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

    /// Build analysis pair cutoffs by merging explicit [analysis.cutoffs] with
    /// fallback from [[constraints.coordination]] cutoffs.
    pub fn analysis_pairs(&self) -> HashMap<String, f64> {
        let mut pairs = HashMap::new();

        // Fallback: coordination constraint cutoffs
        if let Some(ref cc) = self.constraints {
            if let Some(ref coords) = cc.coordination {
                for coord in coords {
                    pairs.entry(coord.pair.clone()).or_insert(coord.cutoff);
                }
            }
        }

        // Explicit analysis cutoffs override
        if let Some(ref ac) = self.analysis {
            if let Some(ref cutoffs) = ac.cutoffs {
                for (pair, &cutoff) in cutoffs {
                    pairs.insert(pair.clone(), cutoff);
                }
            }
        }

        pairs
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
