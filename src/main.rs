#![allow(clippy::needless_range_loop)]
#![allow(clippy::type_complexity)]

use std::path::{Path, PathBuf};
use std::process;

use rsmith::analyze;
use rsmith::config::Config;
use rsmith::epsr::{self, EpsrState};
use rsmith::io;
use rsmith::neutron;
use rsmith::potential::PotentialSet;
use rsmith::rdf;
use rsmith::rmc::{self, DataKind, ExperimentalData, ExperimentalGrData, RmcParams, SqConvention};
use rsmith::sq;
use rsmith::xray;
use rsmith::{log_eprintln, log_println};

fn parse_convention(s: Option<&str>) -> SqConvention {
    match s {
        Some("iq") => SqConvention::Iq,
        Some("fq") => SqConvention::Fq,
        Some("sq") | None => SqConvention::Sq,
        Some(other) => {
            eprintln!(
                "Warning: unknown convention '{}', using 'sq'. Valid: sq, iq, fq",
                other
            );
            SqConvention::Sq
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: rsmith <config.toml> [OPTIONS]");
        eprintln!("Options:");
        eprintln!("  --analyze [structure.xyz]   Run structural analysis");
        eprintln!("  --output-dir DIR            Write output to DIR instead of config directory");
        eprintln!("  --seed N                    Override RNG seed (default: random)");
        eprintln!("  --quiet                     Suppress terminal output (log file only)");
        eprintln!("  --resume                    Resume RMC from checkpoint.dat");
        process::exit(1);
    }

    let config_path = Path::new(&args[1]);
    let analyze_mode = args.iter().any(|a| a == "--analyze");
    let quiet_mode = args.iter().any(|a| a == "--quiet");
    let resume_mode = args.iter().any(|a| a == "--resume");

    // --seed N
    let cli_seed: Option<u64> = args
        .iter()
        .position(|a| a == "--seed")
        .and_then(|pos| args.get(pos + 1))
        .and_then(|s| s.parse().ok());

    // --output-dir DIR
    let cli_output_dir: Option<PathBuf> = args
        .iter()
        .position(|a| a == "--output-dir")
        .and_then(|pos| args.get(pos + 1))
        .map(PathBuf::from);

    // Optional structure path override after --analyze
    let analyze_structure: Option<String> = args
        .iter()
        .position(|a| a == "--analyze")
        .and_then(|pos| args.get(pos + 1))
        .filter(|a| !a.starts_with("--"))
        .cloned();

    // Resolve paths relative to config file location
    let config_dir = config_path
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .to_path_buf();

    // Output directory: --output-dir overrides config_dir
    let output_dir = cli_output_dir.unwrap_or_else(|| config_dir.clone());
    std::fs::create_dir_all(&output_dir).unwrap_or_else(|e| {
        eprintln!(
            "Error creating output directory {}: {}",
            output_dir.display(),
            e
        );
        process::exit(1);
    });

    // Init logging (use separate log file for --analyze to avoid overwriting RMC log)
    rsmith::logging::set_quiet(quiet_mode);
    let log_name = if analyze_mode {
        "rsmith_analyze.log"
    } else {
        "rsmith.log"
    };
    rsmith::logging::init_log_file_in(&output_dir, log_name);

    let cfg = Config::load(config_path).unwrap_or_else(|e| {
        log_eprintln!("Error loading config: {}", e);
        process::exit(1);
    });

    // Load structure (with optional override for --analyze <path>)
    let structure_path = if let Some(ref override_path) = analyze_structure {
        resolve_path(&config_dir, override_path)
    } else {
        resolve_path(&config_dir, &cfg.system.structure)
    };
    log_println!("Loading structure from {:?} ...", structure_path);

    // Auto-detect format from extension when path was overridden
    let format = if analyze_structure.is_some() {
        match structure_path.extension().and_then(|e| e.to_str()) {
            Some("xyz") => "xyz",
            Some("data") | Some("lmp") => "lammps",
            Some("poscar") | Some("vasp") => "poscar",
            _ if structure_path
                .file_name()
                .and_then(|f| f.to_str())
                .is_some_and(|n| {
                    n == "POSCAR"
                        || n == "CONTCAR"
                        || n.starts_with("POSCAR")
                        || n.starts_with("CONTCAR")
                }) =>
            {
                "poscar"
            }
            _ => cfg.system.format.as_str(),
        }
    } else {
        cfg.system.format.as_str()
    };

    let mut config = match format {
        "lammps" => {
            let type_map = cfg.type_map();
            io::read_lammps_data(&structure_path, &type_map).unwrap_or_else(|e| {
                log_eprintln!("Error reading LAMMPS data: {}", e);
                process::exit(1);
            })
        }
        "xyz" => io::read_xyz(&structure_path).unwrap_or_else(|e| {
            log_eprintln!("Error reading XYZ: {}", e);
            process::exit(1);
        }),
        "poscar" => io::read_poscar(&structure_path).unwrap_or_else(|e| {
            log_eprintln!("Error reading POSCAR: {}", e);
            process::exit(1);
        }),
        _ => {
            log_eprintln!("Unknown format: {}", format);
            process::exit(1);
        }
    };

    log_println!(
        "  {} atoms, {} species: {:?}",
        config.atoms.len(),
        config.species.len(),
        config.species
    );
    log_println!(
        "  Box: {:.4} x {:.4} x {:.4} A",
        config.box_lengths[0],
        config.box_lengths[1],
        config.box_lengths[2]
    );
    log_println!(
        "  Volume = {:.1} A^3, rho0 = {:.6} atoms/A^3",
        config.volume(),
        config.number_density()
    );

    // --- Optional density rescaling ---
    if let Some(target_density) = cfg.system.density {
        let current_density = config.mass_density();
        let scale = (current_density / target_density).cbrt();
        log_println!(
            "\nRescaling box to target density {:.4} g/cm^3 (current: {:.4} g/cm^3, scale factor: {:.6})",
            target_density, current_density, scale
        );
        for d in 0..3 {
            config.box_lengths[d] *= scale;
        }
        for atom in &mut config.atoms {
            for d in 0..3 {
                atom.position[d] *= scale;
            }
        }
        log_println!(
            "  New box: {:.4} x {:.4} x {:.4} A",
            config.box_lengths[0],
            config.box_lengths[1],
            config.box_lengths[2]
        );
        log_println!(
            "  New volume = {:.1} A^3, rho0 = {:.6} atoms/A^3, density = {:.4} g/cm^3",
            config.volume(),
            config.number_density(),
            config.mass_density()
        );
    }

    let mut params = cfg.rmc_params();

    // Seed priority: CLI --seed > TOML seed > random from entropy
    if let Some(seed) = cli_seed {
        params.seed = seed;
        log_println!("  RNG seed (CLI): {}", seed);
    } else if cfg.rmc.seed.is_some() {
        log_println!("  RNG seed (config): {}", params.seed);
    } else {
        let seed: u64 = rand::random();
        params.seed = seed;
        log_println!("  RNG seed (random): {}", seed);
    }

    // --- Resume from checkpoint ---
    let resume_state: Option<rmc::RmcState> = if resume_mode {
        let checkpoint_path = output_dir.join("checkpoint.dat");
        if !checkpoint_path.exists() {
            log_eprintln!(
                "Error: --resume specified but {:?} not found",
                checkpoint_path
            );
            process::exit(1);
        }
        log_println!("\nResuming from checkpoint {:?}", checkpoint_path);
        let (rs, ckpt_config) = io::read_checkpoint(&checkpoint_path, &config.species)
            .unwrap_or_else(|e| {
                log_eprintln!("Error reading checkpoint: {}", e);
                process::exit(1);
            });
        log_println!(
            "  Checkpoint: move {}/{}, accepted {}, chi2 = {:.6}, max_step = {:.4}",
            rs.move_count,
            params.max_moves,
            rs.accepted,
            rs.chi2,
            rs.max_step
        );
        if rs.move_count >= params.max_moves {
            log_eprintln!(
                "Checkpoint move_count ({}) >= max_moves ({}), nothing to do.",
                rs.move_count,
                params.max_moves
            );
            process::exit(0);
        }
        config = ckpt_config;
        Some(rs)
    } else {
        None
    };

    let rho0 = config.number_density();

    // --- Analyze mode ---
    if analyze_mode {
        if analyze_structure.is_some() {
            // Explicit path: analyze just that one structure
            run_analysis(&config, &cfg, &output_dir, "analysis");
        } else {
            // No explicit path: analyze starting structure, then refined if it exists
            run_analysis(&config, &cfg, &output_dir, "starting");
            let refined_path = output_dir.join("refined.xyz");
            if refined_path.exists() {
                log_println!("\n{}", "=".repeat(60));
                let refined = io::read_xyz(&refined_path).unwrap_or_else(|e| {
                    log_eprintln!("Error reading refined.xyz: {}", e);
                    process::exit(1);
                });
                log_println!("Loading refined structure from {:?} ...", refined_path);
                log_println!(
                    "  {} atoms, {} species: {:?}",
                    refined.atoms.len(),
                    refined.species.len(),
                    refined.species
                );
                run_analysis(&refined, &cfg, &output_dir, "refined");
            }
        }
        rsmith::logging::flush_log_file();
        return;
    }

    // --- Load experimental data ---
    let mut experiments: Vec<ExperimentalData> = Vec::new();

    if let Some(ref xray_cfg) = cfg.data.xray_sq {
        let path = resolve_path(&config_dir, &xray_cfg.file);
        log_println!("Loading X-ray S(Q) from {:?} ...", path);
        let (q, sq) = io::read_sq_data(&path).unwrap_or_else(|e| {
            log_eprintln!("Error reading X-ray S(Q): {}", e);
            process::exit(1);
        });
        let mut sigma = match xray_cfg.sigma {
            Some(val) => {
                log_println!("  Using constant sigma = {:.4}", val);
                vec![val; sq.len()]
            }
            None => {
                let s = rmc::estimate_sigma(&sq, 10);
                log_println!(
                    "  Estimated sigma from data: min={:.4}, max={:.4}",
                    s.iter().cloned().fold(f64::INFINITY, f64::min),
                    s.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
                );
                s
            }
        };
        let alpha = xray_cfg.sigma_alpha.unwrap_or(0.0);
        if alpha > 0.0 {
            for (i, s) in sigma.iter_mut().enumerate() {
                *s *= 1.0 + alpha * q[i];
            }
            log_println!(
                "  Applied Q-scaling (alpha={:.4}): sigma range {:.4} - {:.4}",
                alpha,
                sigma.iter().cloned().fold(f64::INFINITY, f64::min),
                sigma.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            );
        }
        let weight = xray_cfg.weight.unwrap_or(1.0);
        let fit_min = xray_cfg.fit_min.unwrap_or(0.0);
        let fit_max = xray_cfg.fit_max.unwrap_or(f64::INFINITY);
        log_println!(
            "  {} Q points, Q range: {:.2} - {:.2}, fit range: [{:.2}, {:.2}]",
            q.len(),
            q[0],
            q[q.len() - 1],
            fit_min,
            fit_max
        );
        let convention = parse_convention(xray_cfg.convention.as_deref());
        if convention != SqConvention::Sq {
            log_println!("  Convention: {:?}", convention);
        }
        experiments.push(ExperimentalData {
            q,
            sq,
            sigma,
            weight,
            kind: DataKind::Xray,
            fit_min,
            fit_max,
            convention,
        });
    }

    if let Some(ref neutron_cfg) = cfg.data.neutron_sq {
        let path = resolve_path(&config_dir, &neutron_cfg.file);
        log_println!("Loading neutron S(Q) from {:?} ...", path);
        let (q, sq) = io::read_sq_data(&path).unwrap_or_else(|e| {
            log_eprintln!("Error reading neutron S(Q): {}", e);
            process::exit(1);
        });
        let mut sigma = match neutron_cfg.sigma {
            Some(val) => {
                log_println!("  Using constant sigma = {:.4}", val);
                vec![val; sq.len()]
            }
            None => {
                let s = rmc::estimate_sigma(&sq, 10);
                log_println!(
                    "  Estimated sigma from data: min={:.4}, max={:.4}",
                    s.iter().cloned().fold(f64::INFINITY, f64::min),
                    s.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
                );
                s
            }
        };
        let alpha = neutron_cfg.sigma_alpha.unwrap_or(0.0);
        if alpha > 0.0 {
            for (i, s) in sigma.iter_mut().enumerate() {
                *s *= 1.0 + alpha * q[i];
            }
            log_println!(
                "  Applied Q-scaling (alpha={:.4}): sigma range {:.4} - {:.4}",
                alpha,
                sigma.iter().cloned().fold(f64::INFINITY, f64::min),
                sigma.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            );
        }
        let weight = neutron_cfg.weight.unwrap_or(1.0);
        let fit_min = neutron_cfg.fit_min.unwrap_or(0.0);
        let fit_max = neutron_cfg.fit_max.unwrap_or(f64::INFINITY);
        let convention = parse_convention(neutron_cfg.convention.as_deref());
        if convention != SqConvention::Sq {
            log_println!("  Convention: {:?}", convention);
        }
        experiments.push(ExperimentalData {
            q,
            sq,
            sigma,
            weight,
            kind: DataKind::Neutron,
            fit_min,
            fit_max,
            convention,
        });
    }

    // --- Load experimental g(r) data ---
    let mut gr_datasets: Vec<ExperimentalGrData> = Vec::new();

    if let Some(ref gr_cfg) = cfg.data.xray_gr {
        let path = resolve_path(&config_dir, &gr_cfg.file);
        log_println!("Loading X-ray g(r) from {:?} ...", path);
        let (r, gr) = io::read_sq_data(&path).unwrap_or_else(|e| {
            log_eprintln!("Error reading X-ray g(r): {}", e);
            process::exit(1);
        });
        let sigma = match gr_cfg.sigma {
            Some(val) => {
                log_println!("  Using constant sigma = {:.4}", val);
                vec![val; gr.len()]
            }
            None => {
                let s = rmc::estimate_sigma(&gr, 10);
                log_println!(
                    "  Estimated sigma from data: min={:.4}, max={:.4}",
                    s.iter().cloned().fold(f64::INFINITY, f64::min),
                    s.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
                );
                s
            }
        };
        let weight = gr_cfg.weight.unwrap_or(1.0);
        let fit_min = gr_cfg.fit_min.unwrap_or(0.0);
        let fit_max = gr_cfg.fit_max.unwrap_or(f64::INFINITY);
        // Q_max for inverse FT: use explicit config, or experimental S(Q) Qmax, or model Qmax
        let qmax = gr_cfg.qmax.unwrap_or_else(|| {
            if !experiments.is_empty() {
                *experiments[0].q.last().unwrap_or(&20.0)
            } else {
                *params.q_grid.last().unwrap_or(&20.0)
            }
        });
        let lorch = gr_cfg.lorch.unwrap_or(true);
        log_println!("  {} r points, r range: {:.2} - {:.2} A, fit range: [{:.2}, {:.2}], FT Qmax: {:.2}, Lorch: {}",
            r.len(), r[0], r[r.len() - 1], fit_min, fit_max, qmax, lorch);
        gr_datasets.push(ExperimentalGrData {
            r,
            gr,
            sigma,
            weight,
            fit_min,
            fit_max,
            qmax,
            lorch,
        });
    }

    if experiments.is_empty() && gr_datasets.is_empty() {
        log_eprintln!("No experimental data specified in config!");
        process::exit(1);
    }

    // --- Constraints ---
    let constraints = cfg.constraints();
    if !constraints.min_distances.is_empty() {
        log_println!("Minimum distance constraints:");
        for (pair, d) in &constraints.min_distances {
            log_println!("  {} > {:.2} A", pair, d);
        }
    }
    if !constraints.coordination.is_empty() {
        log_println!("Coordination constraints:");
        for cc in &constraints.coordination {
            log_println!(
                "  {} [{}, {}] within {:.2} A",
                cc.pair,
                cc.min,
                cc.max,
                cc.cutoff
            );
        }
    }

    // --- Build pair potentials ---
    let potential_set = if let Some(ref pot_cfg) = cfg.potential {
        match PotentialSet::from_config(pot_cfg, &config.species, params.rdf_cutoff, &config_dir) {
            Ok(ps) => {
                log_println!(
                    "\nPair potentials (weight = {:.6}, cutoff = {:.1} A):",
                    ps.weight,
                    ps.cutoff
                );
                for pot in &ps.potentials {
                    log_println!(
                        "  {}: {} bins, dr = {:.4} A",
                        pot.pair_label,
                        pot.n_bins,
                        pot.dr
                    );
                }
                Some(ps)
            }
            Err(e) => {
                log_eprintln!("Error building potentials: {}", e);
                process::exit(1);
            }
        }
    } else {
        None
    };

    // --- Save starting structure S(Q) and g(r) (skip on resume) ---
    if resume_state.is_none() {
        log_println!("\nComputing starting structure S(Q) and g(r)...");
        let rdf_dr = params.rdf_cutoff / params.rdf_nbins as f64;
        let histograms = rdf::compute_histograms(&config, params.rdf_nbins, params.rdf_cutoff);
        let partials_gr = rdf::normalise_histograms(&histograms, &config, params.rdf_nbins, rdf_dr);
        let r_grid: Vec<f64> = (0..params.rdf_nbins)
            .map(|i| (i as f64 + 0.5) * rdf_dr)
            .collect();
        let partial_sq =
            sq::compute_all_partial_sq(&r_grid, &partials_gr, rho0, &params.q_grid, params.lorch);
        // Write starting S(Q) for each data type, applying convention offset
        for exp in experiments.iter() {
            let (raw_sq, label) = match exp.kind {
                DataKind::Xray => (
                    xray::compute_xray_sq(&config, &partial_sq, &params.q_grid),
                    "xray",
                ),
                DataKind::Neutron => (
                    neutron::compute_sq(&config, &partial_sq, &params.q_grid),
                    "neutron",
                ),
            };
            let out_sq = exp.convention.transform_array(&raw_sq, &params.q_grid);
            let path = output_dir.join(format!("start_{}_sq.dat", label));
            io::write_sq(&path, &params.q_grid, &out_sq).unwrap();
            log_println!("  Saved starting {} S(Q) to {:?}", label, path);
        }

        // Save starting partial g(r)
        let gr_path = output_dir.join("start_gr.dat");
        {
            let n_types = config.species.len();
            let mut file = std::fs::File::create(&gr_path).unwrap();
            use std::io::Write;
            write!(file, "# r").unwrap();
            for a in 0..n_types {
                for b in a..n_types {
                    write!(file, " g_{}{}", config.species[a], config.species[b]).unwrap();
                }
            }
            writeln!(file).unwrap();
            for bin in 0..params.rdf_nbins {
                write!(file, "{:.6}", r_grid[bin]).unwrap();
                for a in 0..n_types {
                    for b in a..n_types {
                        let pair_idx = config.pair_index(a, b);
                        write!(file, " {:.6}", partials_gr[&pair_idx][bin]).unwrap();
                    }
                }
                writeln!(file).unwrap();
            }
        }
        log_println!("  Saved starting partial g(r) to {:?}", gr_path);

        // Save starting total g(r) via inverse FT
        if !gr_datasets.is_empty() {
            // Raw S(Q) in Faber-Ziman convention for FT
            let has_neutron_ft = experiments
                .iter()
                .any(|e| matches!(e.kind, DataKind::Neutron));
            let sq_for_ft = if has_neutron_ft {
                neutron::compute_sq(&config, &partial_sq, &params.q_grid)
            } else {
                xray::compute_xray_sq(&config, &partial_sq, &params.q_grid)
            };
            let gd0 = &gr_datasets[0];
            let qmax_gr = gd0.qmax;
            let use_lorch = gd0.lorch;
            let dq = if params.q_grid.len() > 1 {
                params.q_grid[1] - params.q_grid[0]
            } else {
                1.0
            };
            let total_gr: Vec<f64> = r_grid
                .iter()
                .map(|&ri| {
                    if ri < 1e-10 {
                        return 1.0;
                    }
                    let pref = dq / (2.0 * std::f64::consts::PI * std::f64::consts::PI * rho0 * ri);
                    let mut val = 1.0;
                    for (k, &qk) in params.q_grid.iter().enumerate() {
                        if qk > qmax_gr {
                            break;
                        }
                        let window = if use_lorch {
                            let arg = std::f64::consts::PI * qk / qmax_gr;
                            if arg > 1e-10 {
                                arg.sin() / arg
                            } else {
                                1.0
                            }
                        } else {
                            1.0
                        };
                        val += pref * qk * window * (sq_for_ft[k] - 1.0) * (qk * ri).sin();
                    }
                    val
                })
                .collect();
            let total_gr_path = output_dir.join("start_total_gr.dat");
            io::write_gr(&total_gr_path, &r_grid, &total_gr).unwrap();
            log_println!("  Saved starting total g(r) to {:?}", total_gr_path);
        }
    }

    // --- RMC refinement ---
    log_println!("\nStarting RMC refinement:");
    log_println!("  max_moves = {}", params.max_moves);
    log_println!("  max_step = {:.4} A", params.max_step);
    log_println!(
        "  RDF cutoff = {:.1} A, {} bins",
        params.rdf_cutoff,
        params.rdf_nbins
    );
    log_println!(
        "  Q grid: {} points up to {:.1} 1/A",
        params.q_grid.len(),
        params.q_grid.last().unwrap_or(&0.0)
    );
    if (params.anneal_start - params.anneal_end).abs() > 1e-10 {
        log_println!(
            "  Annealing: T = {:.2} -> {:.2}",
            params.anneal_start,
            params.anneal_end
        );
    }
    if params.convergence_threshold > 0.0 {
        log_println!(
            "  Convergence: threshold = {:.1e}, window = {} moves",
            params.convergence_threshold,
            params.convergence_window
        );
    }
    log_println!();

    let state = if let Some(ref epsr_cfg) = cfg.epsr {
        // --- EPSR outer loop ---
        let epsr_iters = epsr_cfg.iterations.unwrap_or(10);
        let epsr_feedback = epsr_cfg.feedback.unwrap_or(0.2);
        let epsr_smooth = epsr_cfg.smooth_sigma.unwrap_or(0.02);
        let epsr_kt = epsr_cfg.temperature.unwrap_or(0.025);
        let epsr_min_r = epsr_cfg.min_r.unwrap_or(1.0);
        let epsr_conv = epsr_cfg.convergence.unwrap_or(0.0);
        let epsr_conv_window = epsr_cfg.convergence_window.unwrap_or(3);
        let epsr_moves = epsr_cfg.moves_per_iteration.unwrap_or(params.max_moves);

        log_println!(
            "EPSR mode: {} iterations, feedback = {:.3}, kT = {:.4} eV, smooth_sigma = {:.3} A",
            epsr_iters,
            epsr_feedback,
            epsr_kt,
            epsr_smooth
        );
        if epsr_conv > 0.0 {
            log_println!(
                "  moves/iter = {}, min_r = {:.2} A, convergence = {:.1e} (relative, window = {})",
                epsr_moves,
                epsr_min_r,
                epsr_conv,
                epsr_conv_window
            );
        } else {
            log_println!(
                "  moves/iter = {}, min_r = {:.2} A, convergence = off",
                epsr_moves,
                epsr_min_r
            );
        }

        // Save reference potential set (before EP)
        let reference_potential = potential_set.clone();

        // Initialize EPSR state
        let ep_dr = 0.01; // EP grid spacing (coarser than potential table)
        let ep_cutoff = params.rdf_cutoff.min(
            potential_set
                .as_ref()
                .map_or(params.rdf_cutoff, |p| p.cutoff),
        );
        let mut epsr_state = EpsrState::new(&config.species, ep_cutoff, ep_dr);

        // Load previous EP tables if ep_restart is set
        if let Some(ref restart_dir) = epsr_cfg.ep_restart {
            let restart_path = resolve_path(&config_dir, restart_dir);
            log_println!("  Loading EP restart from {:?}", restart_path);
            match epsr_state.load_potentials(&restart_path) {
                Ok(n) => {
                    log_println!(
                        "  Loaded {}/{} EP tables from previous run",
                        n,
                        epsr_state.n_pairs
                    );
                }
                Err(e) => {
                    log_eprintln!("Error loading EP restart: {}", e);
                    process::exit(1);
                }
            }
        }

        // Find first S(Q) experiment for EPSR EP update (prefer X-ray, fall back to neutron)
        let epsr_exp = experiments
            .iter()
            .find(|e| matches!(e.kind, DataKind::Xray))
            .or_else(|| {
                experiments
                    .iter()
                    .find(|e| matches!(e.kind, DataKind::Neutron))
            });

        // Compute weights matching the experiment type
        let epsr_weights = match epsr_exp {
            Some(e) if matches!(e.kind, DataKind::Xray) => {
                log_println!("  EPSR EP update using X-ray data");
                EpsrState::compute_xray_weights(&config, &params.q_grid)
            }
            Some(e) if matches!(e.kind, DataKind::Neutron) => {
                log_println!("  EPSR EP update using neutron data");
                EpsrState::compute_neutron_weights(&config, params.q_grid.len())
            }
            _ => {
                log_println!("  Warning: no S(Q) data for EPSR EP update");
                vec![]
            }
        };

        // Interpolate experimental S(Q) onto simulation Q grid
        let exp_on_grid: Option<Vec<f64>> =
            epsr_exp.map(|exp| epsr::interpolate_exp_to_grid(&exp.q, &exp.sq, &params.q_grid));

        let mut last_state: Option<rmc::RmcState> = resume_state;
        let mut conv_streak: usize = 0;

        for iter in 0..epsr_iters {
            log_println!("\n{}", "=".repeat(60));
            log_println!("EPSR iteration {}/{}", iter + 1, epsr_iters);
            log_println!("{}", "=".repeat(60));

            // First iteration: use full [rmc] settings (convergence, best-restoration,
            // annealing) to get a well-converged starting structure.
            // Subsequent iterations: equilibrium runs for EP update — no convergence
            // early-stop, no best-restoration, fixed moves_per_iteration.
            let mut iter_params = if iter == 0 {
                log_println!("Using [rmc] settings for initial convergence");
                params.clone()
            } else {
                RmcParams {
                    max_moves: epsr_moves,
                    max_step: params.max_step,
                    checkpoint_every: params.checkpoint_every,
                    seed: params.seed,
                    rdf_cutoff: params.rdf_cutoff,
                    rdf_nbins: params.rdf_nbins,
                    q_grid: params.q_grid.clone(),
                    lorch: params.lorch,
                    print_every: params.print_every,
                    target_acceptance: params.target_acceptance,
                    adjust_step_every: params.adjust_step_every,
                    anneal_start: params.anneal_end, // no annealing
                    anneal_end: params.anneal_end,
                    anneal_steps: 0,
                    convergence_threshold: 0.0,
                    convergence_window: params.convergence_window,
                    restore_best: false,
                }
            };

            // Build combined potential = reference + EP
            let combined = epsr_state.build_combined_potential(
                reference_potential.as_ref(),
                &config.species,
                params.rdf_cutoff,
            );

            // Run MC with combined potential
            let checkpoint_dir = output_dir.clone();
            let checkpoint_fn: Option<Box<dyn Fn(&rmc::RmcState, &rsmith::atoms::Configuration)>> =
                Some(Box::new(move |state, cfg| {
                    let path = checkpoint_dir.join("checkpoint.dat");
                    if let Err(e) = io::write_checkpoint(&path, state, cfg) {
                        log_eprintln!("Warning: checkpoint failed: {}", e);
                    } else {
                        log_println!("  Checkpoint saved at move {}", state.move_count);
                    }
                }));

            // Re-seed to continue deterministically from previous iteration
            if let Some(ref prev) = last_state {
                iter_params.seed = prev.seed.wrapping_add(prev.move_count);
            }

            let state = rmc::run_rmc(
                &mut config,
                &experiments,
                &gr_datasets,
                &constraints,
                &iter_params,
                Some(&combined),
                checkpoint_fn,
                None, // fresh start each EPSR iter (positions carry over in config)
            );

            let acceptance = 100.0 * state.accepted as f64 / state.move_count.max(1) as f64;

            // Extract partials and compute EP update
            if let (Some(ref partial_sq), Some(ref total_sq), Some(ref exp_grid)) =
                (&state.partial_sq, &state.total_sq, &exp_on_grid)
            {
                let nq = params.q_grid.len();
                let delta_partials = EpsrState::compute_residual_partials(
                    partial_sq,
                    total_sq,
                    exp_grid,
                    &epsr_weights,
                    epsr_state.n_pairs,
                    nq,
                );

                let (max_delta, max_ep) = epsr_state.update(
                    &delta_partials,
                    &params.q_grid,
                    rho0,
                    epsr_feedback,
                    epsr_kt,
                    epsr_smooth,
                    epsr_min_r,
                );

                let rel_change = if max_ep > 1e-30 {
                    max_delta / max_ep
                } else {
                    f64::INFINITY
                };

                log_println!(
                    "\nEPSR iter {}: chi2 = {:.6}, max |ΔEP| = {:.6} eV, max |EP| = {:.6} eV, rel = {:.4}, acceptance = {:.1}%",
                    iter + 1,
                    state.chi2,
                    max_delta,
                    max_ep,
                    rel_change,
                    acceptance
                );

                // Write EP files and per-iteration structure
                if let Err(e) = epsr_state.write_potentials(&output_dir) {
                    log_eprintln!("Warning: failed to write EP files: {}", e);
                }
                let iter_xyz = output_dir.join(format!("refined_iter_{}.xyz", iter + 1));
                if let Err(e) = io::write_xyz(&iter_xyz, &config) {
                    log_eprintln!("Warning: failed to write iter structure: {}", e);
                }

                // Convergence check: relative change below threshold for N consecutive iterations
                if epsr_conv > 0.0 {
                    if rel_change < epsr_conv {
                        conv_streak += 1;
                        log_println!(
                            "  Convergence: rel change {:.4} < {:.1e} ({}/{})",
                            rel_change,
                            epsr_conv,
                            conv_streak,
                            epsr_conv_window
                        );
                        if conv_streak >= epsr_conv_window {
                            log_println!(
                                "EPSR converged: relative EP change < {:.1e} for {} consecutive iterations",
                                epsr_conv,
                                epsr_conv_window
                            );
                            last_state = Some(state);
                            break;
                        }
                    } else {
                        if conv_streak > 0 {
                            log_println!(
                                "  Convergence: rel change {:.4} >= {:.1e} — streak reset",
                                rel_change,
                                epsr_conv
                            );
                        }
                        conv_streak = 0;
                    }
                }
            } else {
                log_println!(
                    "\nEPSR iter {}: chi2 = {:.6}, acceptance = {:.1}% (no S(Q) data for EP update)",
                    iter + 1, state.chi2, acceptance
                );
            }

            last_state = Some(state);
        }

        log_println!("\nEPSR refinement complete.");
        last_state.unwrap()
    } else {
        // --- Standard RMC (no EPSR) ---
        let checkpoint_dir = output_dir.clone();
        let checkpoint_fn: Option<Box<dyn Fn(&rmc::RmcState, &rsmith::atoms::Configuration)>> =
            Some(Box::new(move |state, cfg| {
                let path = checkpoint_dir.join("checkpoint.dat");
                if let Err(e) = io::write_checkpoint(&path, state, cfg) {
                    log_eprintln!("Warning: checkpoint failed: {}", e);
                } else {
                    log_println!("  Checkpoint saved at move {}", state.move_count);
                }
            }));

        rmc::run_rmc(
            &mut config,
            &experiments,
            &gr_datasets,
            &constraints,
            &params,
            potential_set.as_ref(),
            checkpoint_fn,
            resume_state,
        )
    };

    // --- Save results ---
    let output_xyz = output_dir.join("refined.xyz");
    log_println!("\nSaving refined structure to {:?}", output_xyz);
    io::write_xyz(&output_xyz, &config).unwrap();

    if cfg.system.output_poscar.unwrap_or(false) {
        let output_poscar = output_dir.join("refined_POSCAR");
        log_println!("Saving refined structure to {:?}", output_poscar);
        io::write_poscar(&output_poscar, &config, "rsmith refined structure").unwrap();
    }

    // Compute and save final S(Q)
    let rdf_dr = params.rdf_cutoff / params.rdf_nbins as f64;
    let histograms = rdf::compute_histograms(&config, params.rdf_nbins, params.rdf_cutoff);
    let partials_gr = rdf::normalise_histograms(&histograms, &config, params.rdf_nbins, rdf_dr);
    let r_grid: Vec<f64> = (0..params.rdf_nbins)
        .map(|i| (i as f64 + 0.5) * rdf_dr)
        .collect();
    let partial_sq =
        sq::compute_all_partial_sq(&r_grid, &partials_gr, rho0, &params.q_grid, params.lorch);

    // Write refined S(Q) for each data type, applying convention offset
    for exp in experiments.iter() {
        let (raw_sq, label) = match exp.kind {
            DataKind::Xray => (
                xray::compute_xray_sq(&config, &partial_sq, &params.q_grid),
                "xray",
            ),
            DataKind::Neutron => (
                neutron::compute_sq(&config, &partial_sq, &params.q_grid),
                "neutron",
            ),
        };
        let out_sq = exp.convention.transform_array(&raw_sq, &params.q_grid);
        let path = output_dir.join(format!("refined_{}_sq.dat", label));
        log_println!("Saving refined {} S(Q) to {:?}", label, path);
        io::write_sq(&path, &params.q_grid, &out_sq).unwrap();
    }

    // Save refined partial g(r)
    let refined_gr_path = output_dir.join("refined_gr.dat");
    {
        let n_types = config.species.len();
        let mut file = std::fs::File::create(&refined_gr_path).unwrap();
        use std::io::Write;
        write!(file, "# r").unwrap();
        for a in 0..n_types {
            for b in a..n_types {
                write!(file, " g_{}{}", config.species[a], config.species[b]).unwrap();
            }
        }
        writeln!(file).unwrap();
        for bin in 0..params.rdf_nbins {
            write!(file, "{:.6}", r_grid[bin]).unwrap();
            for a in 0..n_types {
                for b in a..n_types {
                    let pair_idx = config.pair_index(a, b);
                    write!(file, " {:.6}", partials_gr[&pair_idx][bin]).unwrap();
                }
            }
            writeln!(file).unwrap();
        }
    }
    log_println!("Saving refined partial g(r) to {:?}", refined_gr_path);

    // Save refined total g(r) via inverse FT of S(Q), using same Lorch+Qmax as RMC
    if !gr_datasets.is_empty() {
        // Use primary experiment's weighting for the g(r) FT
        let has_neutron_gr = experiments
            .iter()
            .any(|e| matches!(e.kind, DataKind::Neutron));
        let sq_for_gr = if has_neutron_gr {
            neutron::compute_sq(&config, &partial_sq, &params.q_grid)
        } else {
            xray::compute_xray_sq(&config, &partial_sq, &params.q_grid)
        };
        let gd0 = &gr_datasets[0];
        let qmax_gr = gd0.qmax;
        let use_lorch = gd0.lorch;
        let dq = if params.q_grid.len() > 1 {
            params.q_grid[1] - params.q_grid[0]
        } else {
            1.0
        };
        let n_r_out = params.rdf_nbins;
        let dr_out = params.rdf_cutoff / n_r_out as f64;
        let r_out: Vec<f64> = (0..n_r_out).map(|i| (i as f64 + 0.5) * dr_out).collect();
        let total_gr: Vec<f64> = r_out
            .iter()
            .map(|&ri| {
                if ri < 1e-10 {
                    return 1.0;
                }
                let pref = dq / (2.0 * std::f64::consts::PI * std::f64::consts::PI * rho0 * ri);
                let mut val = 1.0;
                for (k, &qk) in params.q_grid.iter().enumerate() {
                    if qk > qmax_gr {
                        break;
                    }
                    let window = if use_lorch {
                        let arg = std::f64::consts::PI * qk / qmax_gr;
                        if arg > 1e-10 {
                            arg.sin() / arg
                        } else {
                            1.0
                        }
                    } else {
                        1.0
                    };
                    val += pref * qk * window * (sq_for_gr[k] - 1.0) * (qk * ri).sin();
                }
                val
            })
            .collect();
        let output_gr = output_dir.join("refined_total_gr.dat");
        log_println!("Saving refined total g(r) to {:?}", output_gr);
        io::write_gr(&output_gr, &r_out, &total_gr).unwrap();
    }

    log_println!("\nDone. Final chi2 = {:.6}", state.chi2);
    rsmith::logging::flush_log_file();
}

fn run_analysis(
    config: &rsmith::atoms::Configuration,
    cfg: &Config,
    output_dir: &Path,
    label: &str,
) {
    let pair_cutoffs = cfg.analysis_pairs();
    if pair_cutoffs.is_empty() {
        log_eprintln!("No pair cutoffs found for analysis. Add [analysis.cutoffs] or [[constraints.coordination]] to config.");
        process::exit(1);
    }

    let pairs = analyze::build_analysis_pairs(&pair_cutoffs);
    log_println!("\n--- {} structure ---", label);
    log_println!("Analysis pairs:");
    for p in &pairs {
        log_println!(
            "  {}-{}: cutoff = {:.2} A",
            p.species_a,
            p.species_b,
            p.cutoff
        );
    }

    log_println!("\nComputing coordination numbers...");
    let cn_results = analyze::compute_coordination_numbers(config, &pairs);

    let nbins = cfg
        .analysis
        .as_ref()
        .and_then(|a| a.angle_bins)
        .unwrap_or(180);
    let triplet_filter: Option<Vec<String>> =
        cfg.analysis.as_ref().and_then(|a| a.angle_triplets.clone());

    log_println!("Computing bond angle distributions ({} bins)...", nbins);
    let angle_results =
        analyze::compute_bond_angles(config, &pairs, nbins, triplet_filter.as_deref());

    analyze::print_analysis_summary(&cn_results, &angle_results);

    let cn_path = output_dir.join(format!("analysis_{}_cn.dat", label));
    analyze::write_cn_histograms(&cn_path, &cn_results).unwrap_or_else(|e| {
        log_eprintln!("Error writing CN histograms: {}", e);
    });
    log_println!("\nWrote CN histograms to {:?}", cn_path);

    if !angle_results.is_empty() {
        let angle_path = output_dir.join(format!("analysis_{}_angles.dat", label));
        analyze::write_angle_histograms(&angle_path, &angle_results).unwrap_or_else(|e| {
            log_eprintln!("Error writing angle histograms: {}", e);
        });
        log_println!("Wrote angle histograms to {:?}", angle_path);
    }
}

fn resolve_path(base: &Path, relative: &str) -> PathBuf {
    let p = Path::new(relative);
    if p.is_absolute() {
        p.to_path_buf()
    } else {
        base.join(p)
    }
}
