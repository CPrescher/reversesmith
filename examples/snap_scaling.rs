//! Manual atom-count scaling benchmark for the native SNAP RMC backend.

use std::collections::HashMap;
use std::env;
use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::time::Instant;

use rsmith::atoms::{Atom, Configuration};
use rsmith::cells::CellList;
use rsmith::config::{MlBackend, MlPotentialConfig};
use rsmith::energy::EnergyModel;
use rsmith::ml_potential::SnapNativeModel;

struct Arguments {
    coefficient_file: PathBuf,
    parameter_file: PathBuf,
    cells: Vec<usize>,
    warmup: usize,
    trials: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let arguments = parse_arguments()?;
    println!("backend,cells,atoms,init_ms,trials,mean_ms,median_ms,min_ms,max_ms,stddev_ms");

    for cells in arguments.cells {
        let template = diamond_silicon(cells);
        let atom_count = template.atoms.len();
        let config = MlPotentialConfig {
            backend: MlBackend::SnapNative,
            model: None,
            coefficient_file: Some(arguments.coefficient_file.to_string_lossy().into_owned()),
            parameter_file: Some(arguments.parameter_file.to_string_lossy().into_owned()),
            init_args: None,
            weight: Some(1.0),
            cutoff: None,
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
        let mut model = SnapNativeModel::from_config(&config, &structure, Path::new("."))?;
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
        )?;
        let mut samples = run_trials(
            &mut model,
            &mut structure,
            &cell_list,
            arguments.trials,
            true,
        )?;
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
            "snap_native,{cells},{atom_count},{initialization_ms:.3},{},{mean:.6},{median:.6},{:.6},{:.6},{:.6}",
            arguments.trials,
            samples[0],
            samples[samples.len() - 1],
            variance.sqrt(),
        );
    }
    Ok(())
}

fn run_trials(
    model: &mut SnapNativeModel,
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
            return Err("SNAP returned a non-finite trial energy difference".into());
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
    let mut coefficient_file = None;
    let mut parameter_file = None;
    let mut cells = vec![3, 5, 6, 8, 10];
    let mut warmup = 10;
    let mut trials = 100;
    let mut arguments = env::args().skip(1);
    while let Some(argument) = arguments.next() {
        match argument.as_str() {
            "--coeff" => {
                coefficient_file = Some(PathBuf::from(next_value(&mut arguments, "--coeff")?))
            }
            "--param" => {
                parameter_file = Some(PathBuf::from(next_value(&mut arguments, "--param")?))
            }
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
    let coefficient_file = coefficient_file.ok_or("--coeff PATH is required")?;
    let parameter_file = parameter_file.ok_or("--param PATH is required")?;
    if !coefficient_file.is_file() || !parameter_file.is_file() {
        return Err("SNAP coefficient or parameter file does not exist".into());
    }
    if cells.is_empty() || cells.contains(&0) || trials == 0 {
        return Err("cells and trials must be positive".into());
    }
    Ok(Arguments {
        coefficient_file,
        parameter_file,
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
        "SNAP/native atom-count scaling benchmark\n\n\
Usage: cargo run --release --example snap_scaling -- --coeff PATH --param PATH [OPTIONS]\n\n\
Options:\n\
  --coeff PATH    SNAP coefficient file\n\
  --param PATH    SNAP parameter file\n\
  --cells LIST    Diamond-Si supercell edges [default: 3,5,6,8,10]\n\
  --warmup N      Unmeasured trial moves [default: 10]\n\
  --trials N      Measured trial moves [default: 100]\n"
    );
}
