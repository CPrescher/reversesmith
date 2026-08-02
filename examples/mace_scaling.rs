//! Manual CPU thread-scaling benchmark for the MACE/Python RMC backend.
//!
//! This measures the complete rejected-trial path used by RMC for either the
//! full-system or local-cluster MACE energy-delta strategy.

use std::collections::HashMap;
use std::env;
use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::time::Instant;

use rsmith::atoms::{Atom, Configuration};
use rsmith::cells::CellList;
use rsmith::config::{MlBackend, MlEnergyDelta, MlPotentialConfig};
use rsmith::energy::EnergyModel;
use rsmith::ml_potential::MacePythonModel;

const DEFAULT_CELLS: &[usize] = &[3, 5, 6];
const DEFAULT_THREADS: &[usize] = &[1, 2, 4, 8, 10, 14];

struct Arguments {
    model: PathBuf,
    python: String,
    device: String,
    delta: MlEnergyDelta,
    cells: Vec<usize>,
    threads: Vec<usize>,
    warmup: usize,
    trials: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let arguments = parse_arguments()?;
    println!(
        "delta,cells,atoms,threads,baseline_threads,init_ms,trials,central_atoms,cluster_atoms,mean_ms,median_ms,min_ms,max_ms,stddev_ms,speedup_vs_baseline,efficiency"
    );

    for &cells in &arguments.cells {
        let template = diamond_silicon(cells);
        let atom_count = template.atoms.len();
        let positions = template
            .atoms
            .iter()
            .map(|atom| atom.position)
            .collect::<Vec<_>>();
        let cell_list = CellList::new(&positions, &template.box_lengths, 5.0);
        let mut baseline_median = None;
        let baseline_threads = arguments.threads[0];

        for &threads in &arguments.threads {
            let config = MlPotentialConfig {
                backend: MlBackend::MacePython,
                model: Some(arguments.model.to_string_lossy().into_owned()),
                coefficient_file: None,
                parameter_file: None,
                init_args: None,
                weight: Some(1.0),
                cutoff: Some(5.0),
                delta: Some(arguments.delta),
                device: Some(arguments.device.clone()),
                torch_threads: Some(threads),
                python: Some(arguments.python.clone()),
                worker: None,
            };

            let mut structure = template.clone();
            let initialization_start = Instant::now();
            let mut model = MacePythonModel::from_config(&config, &structure, Path::new("."))?;
            let initialization_ms = initialization_start.elapsed().as_secs_f64() * 1.0e3;

            run_trials(
                &mut model,
                &mut structure,
                &cell_list,
                arguments.warmup,
                false,
            )?;
            let mut samples = run_trials(
                &mut model,
                &mut structure,
                &cell_list,
                arguments.trials,
                true,
            )?;
            let (central_atoms, cluster_atoms) = model.last_local_cluster_sizes().unwrap_or((0, 0));
            samples.sort_by(f64::total_cmp);

            let mean = samples.iter().sum::<f64>() / samples.len() as f64;
            let median = if samples.len().is_multiple_of(2) {
                let middle = samples.len() / 2;
                0.5 * (samples[middle - 1] + samples[middle])
            } else {
                samples[samples.len() / 2]
            };
            let variance = samples
                .iter()
                .map(|sample| (sample - mean).powi(2))
                .sum::<f64>()
                / samples.len() as f64;
            let standard_deviation = variance.sqrt();
            if threads == baseline_threads {
                baseline_median = Some(median);
            }
            let baseline = baseline_median.ok_or("failed to record the baseline thread count")?;
            let speedup = baseline / median;
            let efficiency = speedup / (threads as f64 / baseline_threads as f64);

            println!(
                "{},{cells},{atom_count},{threads},{baseline_threads},{initialization_ms:.3},{},{central_atoms},{cluster_atoms},{mean:.3},{median:.3},{:.3},{:.3},{standard_deviation:.3},{speedup:.3},{efficiency:.3}",
                delta_label(arguments.delta),
                arguments.trials,
                samples[0],
                samples[samples.len() - 1],
            );
        }
    }

    Ok(())
}

fn run_trials(
    model: &mut MacePythonModel,
    structure: &mut Configuration,
    cell_list: &CellList,
    count: usize,
    measure: bool,
) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
    let atom_index = 0;
    let old_position = structure.atoms[atom_index].position;
    let mut samples = Vec::with_capacity(count);
    let mut checksum = 0.0;

    for trial in 0..count {
        let offset = 0.001 * (trial + 1) as f64;
        let new_position = [
            old_position[0] + offset,
            old_position[1] - 0.5 * offset,
            old_position[2] + 0.25 * offset,
        ];
        structure.atoms[atom_index].position = new_position;
        let start = Instant::now();
        let delta = model.energy_delta_atom(
            structure,
            atom_index,
            &old_position,
            &new_position,
            cell_list,
            0,
            0,
        );
        model.reject_move(atom_index, &old_position);
        let elapsed_ms = start.elapsed().as_secs_f64() * 1.0e3;
        structure.atoms[atom_index].position = old_position;
        if !delta.is_finite() {
            return Err("MACE returned a non-finite trial energy difference".into());
        }
        checksum += delta;
        if measure {
            samples.push(elapsed_ms);
        }
    }
    black_box(checksum);
    Ok(samples)
}

fn diamond_silicon(cells: usize) -> Configuration {
    let lattice_constant = 5.43;
    let basis = [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.25, 0.75, 0.75],
        [0.75, 0.25, 0.75],
        [0.75, 0.75, 0.25],
    ];
    let mut atoms = Vec::with_capacity(8 * cells.pow(3));
    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                for fractional in basis {
                    atoms.push(Atom {
                        position: [
                            (x as f64 + fractional[0]) * lattice_constant,
                            (y as f64 + fractional[1]) * lattice_constant,
                            (z as f64 + fractional[2]) * lattice_constant,
                        ],
                        species: "Si".to_string(),
                        type_id: 0,
                    });
                }
            }
        }
    }
    Configuration {
        composition: HashMap::from([("Si".to_string(), atoms.len())]),
        atoms,
        box_lengths: [lattice_constant * cells as f64; 3],
        species: vec!["Si".to_string()],
    }
}

fn parse_arguments() -> Result<Arguments, Box<dyn std::error::Error>> {
    let mut model = env::var_os("RSMITH_MACE_TEST_MODEL").map(PathBuf::from);
    let mut python = env::var("RSMITH_MACE_TEST_PYTHON").unwrap_or_else(|_| {
        if Path::new(".venv/bin/python").is_file() {
            ".venv/bin/python".to_string()
        } else {
            "python3".to_string()
        }
    });
    let mut device = "cpu".to_string();
    let mut delta = MlEnergyDelta::Full;
    let mut cells = DEFAULT_CELLS.to_vec();
    let mut threads = DEFAULT_THREADS.to_vec();
    let mut warmup = 2;
    let mut trials = 5;

    let mut arguments = env::args().skip(1);
    while let Some(argument) = arguments.next() {
        match argument.as_str() {
            "--model" => model = Some(PathBuf::from(next_value(&mut arguments, "--model")?)),
            "--python" => python = next_value(&mut arguments, "--python")?,
            "--device" => device = next_value(&mut arguments, "--device")?,
            "--delta" => delta = parse_delta(&next_value(&mut arguments, "--delta")?)?,
            "--cells" => cells = parse_list(&next_value(&mut arguments, "--cells")?)?,
            "--threads" => threads = parse_list(&next_value(&mut arguments, "--threads")?)?,
            "--warmup" => warmup = next_value(&mut arguments, "--warmup")?.parse()?,
            "--trials" => trials = next_value(&mut arguments, "--trials")?.parse()?,
            "-h" | "--help" => {
                print_help();
                std::process::exit(0);
            }
            _ => return Err(format!("unknown argument: {argument}").into()),
        }
    }

    let model =
        model.ok_or("provide --model PATH or set RSMITH_MACE_TEST_MODEL to a MACE checkpoint")?;
    if !model.is_file() {
        return Err(format!("MACE model file does not exist: {}", model.display()).into());
    }
    if cells.is_empty() || cells.contains(&0) {
        return Err("--cells must contain positive integers".into());
    }
    if threads.is_empty() || threads.contains(&0) {
        return Err("--threads must contain positive integers".into());
    }
    if trials == 0 {
        return Err("--trials must be greater than zero".into());
    }

    Ok(Arguments {
        model,
        python,
        device,
        delta,
        cells,
        threads,
        warmup,
        trials,
    })
}

fn next_value(
    arguments: &mut impl Iterator<Item = String>,
    option: &str,
) -> Result<String, Box<dyn std::error::Error>> {
    arguments
        .next()
        .ok_or_else(|| format!("missing value after {option}").into())
}

fn parse_list(value: &str) -> Result<Vec<usize>, Box<dyn std::error::Error>> {
    value
        .split(',')
        .map(|item| Ok(item.trim().parse::<usize>()?))
        .collect()
}

fn parse_delta(value: &str) -> Result<MlEnergyDelta, Box<dyn std::error::Error>> {
    match value {
        "full" => Ok(MlEnergyDelta::Full),
        "local" => Ok(MlEnergyDelta::Local),
        _ => Err(format!("unknown delta mode {value:?}; expected full or local").into()),
    }
}

fn delta_label(delta: MlEnergyDelta) -> &'static str {
    match delta {
        MlEnergyDelta::Full => "full",
        MlEnergyDelta::Local => "local",
    }
}

fn print_help() {
    println!(
        "MACE/Python CPU thread-scaling benchmark\n\n\
Usage: cargo run --release --example mace_scaling -- --model PATH [OPTIONS]\n\n\
Options:\n\
  --model PATH       MACE checkpoint (or RSMITH_MACE_TEST_MODEL)\n\
  --python PATH      Python with MACE installed [default: .venv/bin/python]\n\
  --device DEVICE    MACE device [default: cpu]\n\
  --delta MODE       Energy delta: full or local [default: full]\n\
  --cells LIST       Diamond-Si supercell edges [default: 3,5,6]\n\
                     Atom count is 8*cells^3: 216,1000,1728\n\
  --threads LIST     PyTorch thread counts [default: 1,2,4,8,10,14]\n\
  --warmup N         Unmeasured trial moves per case [default: 2]\n\
  --trials N         Measured trial moves per case [default: 5]\n"
    );
}
