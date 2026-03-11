use std::path::{Path, PathBuf};
use std::process;

use reversesmith::analyze;
use reversesmith::config::Config;
use reversesmith::io;
use reversesmith::potential::PotentialSet;
use reversesmith::rdf;
use reversesmith::rmc::{self, DataKind, ExperimentalData, ExperimentalGrData, RmcParams};
use reversesmith::sq;
use reversesmith::xray;
use reversesmith::{log_println, log_eprintln};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: reversesmith <config.toml> [--compute-sq-only] [--analyze [structure.xyz]]");
        process::exit(1);
    }

    reversesmith::logging::init_log_file();

    let config_path = Path::new(&args[1]);
    let compute_sq_only = args.iter().any(|a| a == "--compute-sq-only");
    let analyze_mode = args.iter().any(|a| a == "--analyze");
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
        config.box_lengths[0], config.box_lengths[1], config.box_lengths[2]
    );
    log_println!(
        "  Volume = {:.1} A^3, rho0 = {:.6} atoms/A^3",
        config.volume(),
        config.number_density()
    );

    let params = cfg.rmc_params();
    let rho0 = config.number_density();

    // --- Analyze mode ---
    if analyze_mode {
        if analyze_structure.is_some() {
            // Explicit path: analyze just that one structure
            run_analysis(&config, &cfg, &config_dir, "analysis");
        } else {
            // No explicit path: analyze starting structure, then refined if it exists
            run_analysis(&config, &cfg, &config_dir, "starting");
            let refined_path = config_dir.join("refined.xyz");
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
                run_analysis(&refined, &cfg, &config_dir, "refined");
            }
        }
        reversesmith::logging::flush_log_file();
        return;
    }

    // --- Compute S(Q) mode ---
    if compute_sq_only {
        compute_sq_and_exit(&config, &params, rho0, &config_dir, &cfg);
        reversesmith::logging::flush_log_file();
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
        let sigma_val = xray_cfg.sigma.unwrap_or(0.01);
        let sigma = vec![sigma_val; sq.len()];
        let weight = xray_cfg.weight.unwrap_or(1.0);
        let fit_min = xray_cfg.fit_min.unwrap_or(0.0);
        let fit_max = xray_cfg.fit_max.unwrap_or(f64::INFINITY);
        log_println!("  {} Q points, Q range: {:.2} - {:.2}, fit range: [{:.2}, {:.2}]",
            q.len(), q[0], q[q.len() - 1], fit_min, fit_max);
        experiments.push(ExperimentalData {
            q,
            sq,
            sigma,
            weight,
            kind: DataKind::Xray,
            fit_min,
            fit_max,
        });
    }

    if let Some(ref neutron_cfg) = cfg.data.neutron_sq {
        let path = resolve_path(&config_dir, &neutron_cfg.file);
        log_println!("Loading neutron S(Q) from {:?} ...", path);
        let (q, sq) = io::read_sq_data(&path).unwrap_or_else(|e| {
            log_eprintln!("Error reading neutron S(Q): {}", e);
            process::exit(1);
        });
        let sigma_val = neutron_cfg.sigma.unwrap_or(0.01);
        let sigma = vec![sigma_val; sq.len()];
        let weight = neutron_cfg.weight.unwrap_or(1.0);
        let fit_min = neutron_cfg.fit_min.unwrap_or(0.0);
        let fit_max = neutron_cfg.fit_max.unwrap_or(f64::INFINITY);
        experiments.push(ExperimentalData {
            q,
            sq,
            sigma,
            weight,
            kind: DataKind::Neutron,
            fit_min,
            fit_max,
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
        let sigma_val = gr_cfg.sigma.unwrap_or(0.01);
        let sigma = vec![sigma_val; gr.len()];
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
                cc.pair, cc.min, cc.max, cc.cutoff
            );
        }
    }

    // --- Build pair potentials ---
    let potential_set = if let Some(ref pot_cfg) = cfg.potential {
        match PotentialSet::from_config(pot_cfg, &config.species, params.rdf_cutoff, &config_dir) {
            Ok(ps) => {
                log_println!("\nPair potentials (weight = {:.6}, cutoff = {:.1} A):", ps.weight, ps.cutoff);
                for pot in &ps.potentials {
                    log_println!("  {}: {} bins, dr = {:.4} A", pot.pair_label, pot.n_bins, pot.dr);
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

    // --- RMC refinement ---
    log_println!("\nStarting RMC refinement:");
    log_println!("  max_moves = {}", params.max_moves);
    log_println!("  max_step = {:.4} A", params.max_step);
    log_println!("  RDF cutoff = {:.1} A, {} bins", params.rdf_cutoff, params.rdf_nbins);
    log_println!("  Q grid: {} points up to {:.1} 1/A", params.q_grid.len(), params.q_grid.last().unwrap_or(&0.0));
    if (params.anneal_start - params.anneal_end).abs() > 1e-10 {
        log_println!("  Annealing: T = {:.2} -> {:.2}", params.anneal_start, params.anneal_end);
    }
    if params.convergence_threshold > 0.0 {
        log_println!("  Convergence: threshold = {:.1e}, window = {} moves", params.convergence_threshold, params.convergence_window);
    }
    log_println!();

    let checkpoint_dir = config_dir.clone();
    let checkpoint_fn: Option<Box<dyn Fn(&rmc::RmcState, &reversesmith::atoms::Configuration)>> =
        Some(Box::new(move |state, cfg| {
            let path = checkpoint_dir.join("checkpoint.dat");
            if let Err(e) = io::write_checkpoint(&path, state, cfg) {
                log_eprintln!("Warning: checkpoint failed: {}", e);
            } else {
                log_println!("  Checkpoint saved at move {}", state.move_count);
            }
        }));

    let state = rmc::run_rmc(
        &mut config,
        &experiments,
        &gr_datasets,
        &constraints,
        &params,
        potential_set.as_ref(),
        checkpoint_fn,
    );

    // --- Save results ---
    let output_xyz = config_dir.join("refined.xyz");
    log_println!("\nSaving refined structure to {:?}", output_xyz);
    io::write_xyz(&output_xyz, &config).unwrap();

    // Compute and save final S(Q)
    let rdf_dr = params.rdf_cutoff / params.rdf_nbins as f64;
    let histograms = rdf::compute_histograms(&config, params.rdf_nbins, params.rdf_cutoff);
    let partials_gr =
        rdf::normalise_histograms(&histograms, &config, params.rdf_nbins, rdf_dr);
    let r_grid: Vec<f64> = (0..params.rdf_nbins)
        .map(|i| (i as f64 + 0.5) * rdf_dr)
        .collect();
    let partial_sq =
        sq::compute_all_partial_sq(&r_grid, &partials_gr, rho0, &params.q_grid, params.lorch);
    let sx = xray::compute_xray_sq(&config, &partial_sq, &params.q_grid);

    let output_sq = config_dir.join("refined_sq.dat");
    log_println!("Saving refined S(Q) to {:?}", output_sq);
    io::write_sq(&output_sq, &params.q_grid, &sx).unwrap();

    // Save refined total X-ray g(r) via inverse FT of S_X(Q), using same Lorch+Qmax as RMC
    if !gr_datasets.is_empty() {
        let gd0 = &gr_datasets[0];
        let qmax_gr = gd0.qmax;
        let use_lorch = gd0.lorch;
        let dq = if params.q_grid.len() > 1 { params.q_grid[1] - params.q_grid[0] } else { 1.0 };
        let n_r_out = params.rdf_nbins;
        let dr_out = params.rdf_cutoff / n_r_out as f64;
        let r_out: Vec<f64> = (0..n_r_out).map(|i| (i as f64 + 0.5) * dr_out).collect();
        let total_gr: Vec<f64> = r_out.iter().map(|&ri| {
            if ri < 1e-10 { return 1.0; }
            let pref = dq / (2.0 * std::f64::consts::PI * std::f64::consts::PI * rho0 * ri);
            let mut val = 1.0;
            for (k, &qk) in params.q_grid.iter().enumerate() {
                if qk > qmax_gr { break; }
                let window = if use_lorch {
                    let arg = std::f64::consts::PI * qk / qmax_gr;
                    if arg > 1e-10 { arg.sin() / arg } else { 1.0 }
                } else {
                    1.0
                };
                val += pref * qk * window * (sx[k] - 1.0) * (qk * ri).sin();
            }
            val
        }).collect();
        let output_gr = config_dir.join("refined_total_gr.dat");
        log_println!("Saving refined total X-ray g(r) to {:?}", output_gr);
        io::write_gr(&output_gr, &r_out, &total_gr).unwrap();
    }

    log_println!("\nDone. Final chi2 = {:.6}", state.chi2);
    reversesmith::logging::flush_log_file();
}

/// Compute S(Q) from the initial structure and exit (no RMC).
fn compute_sq_and_exit(
    config: &reversesmith::atoms::Configuration,
    params: &RmcParams,
    rho0: f64,
    config_dir: &Path,
    cfg: &Config,
) {
    log_println!("\nComputing partial RDFs...");
    let rdfs = rdf::compute_partial_rdfs(config, params.rdf_nbins, params.rdf_cutoff);

    // Print pair labels
    let n_types = config.species.len();
    for a in 0..n_types {
        for b in a..n_types {
            let pair_idx = config.pair_index(a, b);
            let gr = &rdfs.partials[&pair_idx];
            let max_val = gr.iter().cloned().fold(0.0f64, f64::max);
            let max_pos = rdfs.r[gr
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                .unwrap()
                .0];
            log_println!(
                "  {}-{}: max g(r) = {:.2} at r = {:.2} A",
                config.species[a], config.species[b], max_val, max_pos
            );
        }
    }

    log_println!("\nComputing partial S(Q)...");
    let partial_sq = sq::compute_all_partial_sq(
        &rdfs.r,
        &rdfs.partials,
        rho0,
        &params.q_grid,
        params.lorch,
    );

    log_println!("Computing total X-ray S(Q)...");
    let sx = xray::compute_xray_sq(config, &partial_sq, &params.q_grid);

    let output_sq = config_dir.join("computed_sq.dat");
    log_println!("Saving computed S(Q) to {:?}", output_sq);
    io::write_sq(&output_sq, &params.q_grid, &sx).unwrap();

    // Also save partial g(r) for validation
    let output_gr = config_dir.join("computed_gr.dat");
    {
        let mut file = std::fs::File::create(&output_gr).unwrap();
        use std::io::Write;
        write!(file, "# r").unwrap();
        for a in 0..n_types {
            for b in a..n_types {
                write!(file, " g_{}{}", config.species[a], config.species[b]).unwrap();
            }
        }
        writeln!(file).unwrap();
        for bin in 0..rdfs.nbins {
            write!(file, "{:.6}", rdfs.r[bin]).unwrap();
            for a in 0..n_types {
                for b in a..n_types {
                    let pair_idx = config.pair_index(a, b);
                    write!(file, " {:.6}", rdfs.partials[&pair_idx][bin]).unwrap();
                }
            }
            writeln!(file).unwrap();
        }
    }
    log_println!("Saving partial g(r) to {:?}", output_gr);

    // Compute and save total X-ray g(r) via inverse FT (with Lorch+Qmax from g(r) config)
    {
        let qmax_gr = cfg.data.xray_gr.as_ref()
            .and_then(|gc| gc.qmax)
            .unwrap_or_else(|| *params.q_grid.last().unwrap_or(&20.0));
        let use_lorch = cfg.data.xray_gr.as_ref()
            .and_then(|gc| gc.lorch)
            .unwrap_or(true);
        let dq = if params.q_grid.len() > 1 { params.q_grid[1] - params.q_grid[0] } else { 1.0 };
        let n_r_out = params.rdf_nbins;
        let dr_out = params.rdf_cutoff / n_r_out as f64;
        let r_out: Vec<f64> = (0..n_r_out).map(|i| (i as f64 + 0.5) * dr_out).collect();
        let total_gr: Vec<f64> = r_out.iter().map(|&ri| {
            if ri < 1e-10 { return 1.0; }
            let pref = dq / (2.0 * std::f64::consts::PI * std::f64::consts::PI * rho0 * ri);
            let mut val = 1.0;
            for (k, &qk) in params.q_grid.iter().enumerate() {
                if qk > qmax_gr { break; }
                let window = if use_lorch {
                    let arg = std::f64::consts::PI * qk / qmax_gr;
                    if arg > 1e-10 { arg.sin() / arg } else { 1.0 }
                } else {
                    1.0
                };
                val += pref * qk * window * (sx[k] - 1.0) * (qk * ri).sin();
            }
            val
        }).collect();
        let output_total_gr = config_dir.join("computed_total_gr.dat");
        log_println!("Saving total X-ray g(r) to {:?}", output_total_gr);
        io::write_gr(&output_total_gr, &r_out, &total_gr).unwrap();
    }

    // Compare with experimental if available
    if let Some(ref xray_cfg) = cfg.data.xray_sq {
        let path = resolve_path(config_dir, &xray_cfg.file);
        if path.exists() {
            let (q_exp, sq_exp) = io::read_sq_data(&path).unwrap();
            let mut chi2 = 0.0;
            let mut count = 0;
            let sigma = xray_cfg.sigma.unwrap_or(0.01);
            for (i, &qe) in q_exp.iter().enumerate() {
                if qe >= params.q_grid[0] && qe <= *params.q_grid.last().unwrap() {
                    let sq_calc = interp(&params.q_grid, &sx, qe);
                    let diff = sq_calc - sq_exp[i];
                    chi2 += diff * diff / (sigma * sigma);
                    count += 1;
                }
            }
            log_println!("\nChi2 vs experimental X-ray S(Q): {:.2} ({} points, per-point: {:.6})", chi2, count, chi2 / count.max(1) as f64);
        }
    }
}

fn run_analysis(
    config: &reversesmith::atoms::Configuration,
    cfg: &Config,
    config_dir: &Path,
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
        log_println!("  {}-{}: cutoff = {:.2} A", p.species_a, p.species_b, p.cutoff);
    }

    log_println!("\nComputing coordination numbers...");
    let cn_results = analyze::compute_coordination_numbers(config, &pairs);

    let nbins = cfg
        .analysis
        .as_ref()
        .and_then(|a| a.angle_bins)
        .unwrap_or(180);
    let triplet_filter: Option<Vec<String>> = cfg
        .analysis
        .as_ref()
        .and_then(|a| a.angle_triplets.clone());

    log_println!("Computing bond angle distributions ({} bins)...", nbins);
    let angle_results = analyze::compute_bond_angles(
        config,
        &pairs,
        nbins,
        triplet_filter.as_deref(),
    );

    analyze::print_analysis_summary(&cn_results, &angle_results);

    let cn_path = config_dir.join(format!("analysis_{}_cn.dat", label));
    analyze::write_cn_histograms(&cn_path, &cn_results).unwrap_or_else(|e| {
        log_eprintln!("Error writing CN histograms: {}", e);
    });
    log_println!("\nWrote CN histograms to {:?}", cn_path);

    if !angle_results.is_empty() {
        let angle_path = config_dir.join(format!("analysis_{}_angles.dat", label));
        analyze::write_angle_histograms(&angle_path, &angle_results).unwrap_or_else(|e| {
            log_eprintln!("Error writing angle histograms: {}", e);
        });
        log_println!("Wrote angle histograms to {:?}", angle_path);
    }
}

fn interp(x: &[f64], y: &[f64], xi: f64) -> f64 {
    if xi <= x[0] {
        return y[0];
    }
    if xi >= x[x.len() - 1] {
        return y[y.len() - 1];
    }
    let mut lo = 0;
    let mut hi = x.len() - 1;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if x[mid] <= xi {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let t = (xi - x[lo]) / (x[hi] - x[lo]);
    y[lo] + t * (y[hi] - y[lo])
}

fn resolve_path(base: &Path, relative: &str) -> PathBuf {
    let p = Path::new(relative);
    if p.is_absolute() {
        p.to_path_buf()
    } else {
        base.join(p)
    }
}
