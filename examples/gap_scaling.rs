//! Manual atom-count scaling benchmark for the GAP/QUIP RMC backend.

use std::collections::HashMap;
use std::env;
use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::time::Instant;

use rsmith::atoms::{Atom, Configuration};
use rsmith::cells::CellList;
use rsmith::config::{MlBackend, MlPotentialConfig};
use rsmith::energy::EnergyModel;
use rsmith::ml_potential::GapQuipModel;

struct Arguments {
    model: PathBuf,
    init_args: Option<String>,
    cutoff: f64,
    cells: Vec<usize>,
    warmup: usize,
    trials: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let arguments = parse_arguments()?;
    let omp_threads = env::var("OMP_NUM_THREADS").unwrap_or_else(|_| "default".to_string());
    let openblas_threads =
        env::var("OPENBLAS_NUM_THREADS").unwrap_or_else(|_| "default".to_string());
    println!(
        "backend,omp_threads,openblas_threads,cells,atoms,init_ms,trials,mean_ms,median_ms,min_ms,max_ms,stddev_ms"
    );

    for cells in arguments.cells {
        let template = diamond_silicon(cells);
        let atom_count = template.atoms.len();
        let config = MlPotentialConfig {
            backend: MlBackend::GapQuip,
            model: Some(arguments.model.to_string_lossy().into_owned()),
            coefficient_file: None,
            parameter_file: None,
            init_args: arguments.init_args.clone(),
            weight: Some(1.0),
            cutoff: Some(arguments.cutoff),
            delta: None,
            device: None,
            torch_threads: None,
            dtype: None,
            compile_mode: None,
            python: None,
            worker: None,
        };
        let mut structure = template.clone();
        let initialization_start = Instant::now();
        let mut model = GapQuipModel::from_config(&config, &structure, Path::new("."))?;
        let initialization_ms = initialization_start.elapsed().as_secs_f64() * 1.0e3;
        let positions = structure
            .atoms
            .iter()
            .map(|atom| atom.position)
            .collect::<Vec<_>>();
        let cell_list = CellList::new(&positions, &structure.box_lengths, model.cutoff());

        run_trials(
            &mut model,
            &mut structure,
            &cell_list,
            arguments.warmup,
            false,
        );
        let mut samples = run_trials(
            &mut model,
            &mut structure,
            &cell_list,
            arguments.trials,
            true,
        );
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

        println!(
            "gap_quip,{omp_threads},{openblas_threads},{cells},{atom_count},{initialization_ms:.3},{},{mean:.6},{median:.6},{:.6},{:.6},{:.6}",
            arguments.trials,
            samples[0],
            samples[samples.len() - 1],
            variance.sqrt(),
        );
    }
    Ok(())
}

fn run_trials(
    model: &mut GapQuipModel,
    structure: &mut Configuration,
    cell_list: &CellList,
    count: usize,
    measure: bool,
) -> Vec<f64> {
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
        assert!(
            delta.is_finite(),
            "GAP/QUIP returned a non-finite trial energy difference"
        );
        checksum += delta;
        if measure {
            samples.push(elapsed_ms);
        }
    }
    black_box(checksum);
    samples
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
    let mut model = env::var_os("RSMITH_GAP_TEST_MODEL").map(PathBuf::from);
    let mut init_args = env::var("RSMITH_GAP_TEST_INIT_ARGS").ok();
    let mut cutoff = 5.0;
    let mut cells = vec![3, 5, 6, 8, 10];
    let mut warmup = 3;
    let mut trials = 7;
    let mut arguments = env::args().skip(1);
    while let Some(argument) = arguments.next() {
        match argument.as_str() {
            "--model" => model = Some(PathBuf::from(next_value(&mut arguments, "--model")?)),
            "--init-args" => init_args = Some(next_value(&mut arguments, "--init-args")?),
            "--cutoff" => cutoff = next_value(&mut arguments, "--cutoff")?.parse()?,
            "--cells" => cells = parse_list(&next_value(&mut arguments, "--cells")?)?,
            "--warmup" => warmup = next_value(&mut arguments, "--warmup")?.parse()?,
            "--trials" => trials = next_value(&mut arguments, "--trials")?.parse()?,
            "-h" | "--help" => {
                print_help();
                std::process::exit(0);
            }
            _ => return Err(format!("unknown argument: {argument}").into()),
        }
    }
    let model = model.ok_or("provide --model PATH or set RSMITH_GAP_TEST_MODEL")?;
    if !model.is_file() {
        return Err(format!("GAP model file does not exist: {}", model.display()).into());
    }
    if cutoff <= 0.0 || cells.is_empty() || cells.contains(&0) || trials == 0 {
        return Err("cutoff, cells, and trials must be positive".into());
    }
    Ok(Arguments {
        model,
        init_args,
        cutoff,
        cells,
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

fn print_help() {
    println!(
        "GAP/QUIP atom-count scaling benchmark\n\n\
Usage: cargo run --release --features gap-quip --example gap_scaling -- --model PATH [OPTIONS]\n\n\
Options:\n\
  --model PATH       GAP XML model (or RSMITH_GAP_TEST_MODEL)\n\
  --init-args TEXT   Optional QUIP initialization string\n\
  --cutoff FLOAT     Maximum GAP descriptor cutoff in A [default: 5.0]\n\
  --cells LIST       Diamond-Si supercell edges [default: 3,5,6,8,10]\n\
  --warmup N         Unmeasured rejected trial moves [default: 3]\n\
  --trials N         Measured rejected trial moves [default: 7]\n"
    );
}
