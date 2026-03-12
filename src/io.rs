use std::collections::HashMap;
use std::fs;
use std::io::{self, BufRead, Write};
use std::path::Path;

use crate::atoms::{Atom, Configuration};
use crate::rmc::RmcState;

/// Read a LAMMPS data file into a Configuration.
///
/// Expects format: "Atoms # charge"
///   atom_id type_id charge x y z [ix iy iz]
///
/// `type_map` maps LAMMPS type integers to species names, e.g. {1: "Ca", 2: "Si", 3: "O"}.
pub fn read_lammps_data(path: &Path, type_map: &HashMap<u32, String>) -> io::Result<Configuration> {
    let content = fs::read_to_string(path)?;
    let lines: Vec<&str> = content.lines().collect();

    let mut box_lo = [0.0f64; 3];
    let mut box_hi = [0.0f64; 3];
    let mut atoms = Vec::new();

    // Determine species ordering from type_map (sorted by LAMMPS type id)
    let mut type_keys: Vec<u32> = type_map.keys().copied().collect();
    type_keys.sort();
    let species: Vec<String> = type_keys.iter().map(|k| type_map[k].clone()).collect();

    // Build reverse map: species name -> type_id (0-based index into species vec)
    let species_to_id: HashMap<String, usize> = species
        .iter()
        .enumerate()
        .map(|(i, s)| (s.clone(), i))
        .collect();

    let mut in_atoms_section = false;
    for line in &lines {
        let line = line.trim();
        if line.is_empty() {
            if in_atoms_section && !atoms.is_empty() {
                break; // End of atoms section
            }
            continue;
        }

        // Parse box bounds
        if line.contains("xlo xhi") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            box_lo[0] = parts[0].parse().unwrap();
            box_hi[0] = parts[1].parse().unwrap();
        } else if line.contains("ylo yhi") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            box_lo[1] = parts[0].parse().unwrap();
            box_hi[1] = parts[1].parse().unwrap();
        } else if line.contains("zlo zhi") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            box_lo[2] = parts[0].parse().unwrap();
            box_hi[2] = parts[1].parse().unwrap();
        }

        // Detect atoms section
        if line.starts_with("Atoms") {
            in_atoms_section = true;
            continue;
        }

        if in_atoms_section {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 6 {
                continue;
            }
            // atom_id type charge x y z [ix iy iz]
            let lammps_type: u32 = match parts[1].parse() {
                Ok(v) => v,
                Err(_) => continue,
            };
            let x: f64 = parts[3].parse().unwrap();
            let y: f64 = parts[4].parse().unwrap();
            let z: f64 = parts[5].parse().unwrap();

            let species_name = type_map
                .get(&lammps_type)
                .unwrap_or_else(|| panic!("Unknown LAMMPS type {}", lammps_type));
            let type_id = species_to_id[species_name];

            atoms.push(Atom {
                position: [x, y, z],
                species: species_name.clone(),
                type_id,
            });
        }
    }

    let box_lengths = [
        box_hi[0] - box_lo[0],
        box_hi[1] - box_lo[1],
        box_hi[2] - box_lo[2],
    ];

    // Count composition
    let mut composition: HashMap<String, usize> = HashMap::new();
    for atom in &atoms {
        *composition.entry(atom.species.clone()).or_insert(0) += 1;
    }

    let mut config = Configuration {
        atoms,
        box_lengths,
        species,
        composition,
    };
    config.wrap_pbc();
    Ok(config)
}

/// Read XYZ file format.
pub fn read_xyz(path: &Path) -> io::Result<Configuration> {
    let content = fs::read_to_string(path)?;
    let lines: Vec<&str> = content.lines().collect();

    let n_atoms: usize = lines[0].trim().parse().map_err(|e| {
        io::Error::new(io::ErrorKind::InvalidData, format!("Bad atom count: {}", e))
    })?;

    // Try to parse box from comment line: "Lattice="Lx 0 0 0 Ly 0 0 0 Lz" ..." or "Lx Ly Lz"
    let comment = lines[1].trim();
    let box_lengths = parse_xyz_box(comment).unwrap_or([100.0, 100.0, 100.0]);

    let mut atoms = Vec::with_capacity(n_atoms);
    let mut species_set: Vec<String> = Vec::new();

    for i in 0..n_atoms {
        let parts: Vec<&str> = lines[2 + i].split_whitespace().collect();
        let name = parts[0].to_string();
        let x: f64 = parts[1].parse().unwrap();
        let y: f64 = parts[2].parse().unwrap();
        let z: f64 = parts[3].parse().unwrap();

        if !species_set.contains(&name) {
            species_set.push(name.clone());
        }
        let type_id = species_set.iter().position(|s| s == &name).unwrap();

        atoms.push(Atom {
            position: [x, y, z],
            species: name,
            type_id,
        });
    }

    let mut composition: HashMap<String, usize> = HashMap::new();
    for atom in &atoms {
        *composition.entry(atom.species.clone()).or_insert(0) += 1;
    }

    let mut config = Configuration {
        atoms,
        box_lengths,
        species: species_set,
        composition,
    };
    config.wrap_pbc();
    Ok(config)
}

fn parse_xyz_box(comment: &str) -> Option<[f64; 3]> {
    // Extended XYZ: Lattice="Lx 0 0 0 Ly 0 0 0 Lz"
    if let Some(start) = comment.find("Lattice=\"") {
        let rest = &comment[start + 9..];
        if let Some(end) = rest.find('"') {
            let vals: Vec<f64> = rest[..end]
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if vals.len() >= 9 {
                return Some([vals[0], vals[4], vals[8]]);
            }
        }
    }
    // Simple: "Lx Ly Lz"
    let parts: Vec<f64> = comment
        .split_whitespace()
        .filter_map(|s| s.parse().ok())
        .collect();
    if parts.len() >= 3 {
        return Some([parts[0], parts[1], parts[2]]);
    }
    None
}

/// Read two-column Q, S(Q) experimental data (whitespace-separated, # comments).
pub fn read_sq_data(path: &Path) -> io::Result<(Vec<f64>, Vec<f64>)> {
    let file = fs::File::open(path)?;
    let reader = io::BufReader::new(file);
    let mut q_vals = Vec::new();
    let mut sq_vals = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let line = line.trim().to_string();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            if let (Ok(q), Ok(s)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                q_vals.push(q);
                sq_vals.push(s);
            }
        }
    }
    Ok((q_vals, sq_vals))
}

/// Write configuration as XYZ file.
pub fn write_xyz(path: &Path, config: &Configuration) -> io::Result<()> {
    let mut file = fs::File::create(path)?;
    writeln!(file, "{}", config.atoms.len())?;
    writeln!(
        file,
        "Lattice=\"{:.6} 0 0 0 {:.6} 0 0 0 {:.6}\"",
        config.box_lengths[0], config.box_lengths[1], config.box_lengths[2]
    )?;
    for atom in &config.atoms {
        writeln!(
            file,
            "{} {:.8} {:.8} {:.8}",
            atom.species, atom.position[0], atom.position[1], atom.position[2]
        )?;
    }
    Ok(())
}

/// Write two-column Q, S(Q) data to file.
pub fn write_sq(path: &Path, q: &[f64], sq: &[f64]) -> io::Result<()> {
    let mut file = fs::File::create(path)?;
    writeln!(file, "# Q(1/A)  S(Q)")?;
    for (qi, si) in q.iter().zip(sq.iter()) {
        writeln!(file, "{:.6} {:.6}", qi, si)?;
    }
    Ok(())
}

/// Write two-column r, g(r) data to file.
pub fn write_gr(path: &Path, r: &[f64], gr: &[f64]) -> io::Result<()> {
    let mut file = fs::File::create(path)?;
    writeln!(file, "# r(A)  g(r)")?;
    for (ri, gi) in r.iter().zip(gr.iter()) {
        writeln!(file, "{:.6} {:.6}", ri, gi)?;
    }
    Ok(())
}

/// Serialize RMC state to a checkpoint file (simple text format).
pub fn write_checkpoint(path: &Path, state: &RmcState, config: &Configuration) -> io::Result<()> {
    let mut file = fs::File::create(path)?;
    writeln!(file, "# rsmith checkpoint")?;
    writeln!(file, "move_count {}", state.move_count)?;
    writeln!(file, "accepted {}", state.accepted)?;
    writeln!(file, "chi2 {:.10}", state.chi2)?;
    writeln!(file, "max_step {:.10}", state.max_step)?;
    writeln!(file, "seed {}", state.seed)?;
    writeln!(file, "natoms {}", config.atoms.len())?;
    writeln!(
        file,
        "box {:.10} {:.10} {:.10}",
        config.box_lengths[0], config.box_lengths[1], config.box_lengths[2]
    )?;
    for atom in &config.atoms {
        writeln!(
            file,
            "{} {} {:.10} {:.10} {:.10}",
            atom.species, atom.type_id, atom.position[0], atom.position[1], atom.position[2]
        )?;
    }
    Ok(())
}

/// Read checkpoint file, returning (RmcState, Configuration).
pub fn read_checkpoint(path: &Path, species: &[String]) -> io::Result<(RmcState, Configuration)> {
    let content = fs::read_to_string(path)?;
    let mut move_count = 0u64;
    let mut accepted = 0u64;
    let mut chi2 = 0.0f64;
    let mut max_step = 0.1f64;
    let mut seed = 42u64;
    let mut box_lengths = [0.0f64; 3];
    let mut atoms = Vec::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        match parts[0] {
            "move_count" => move_count = parts[1].parse().unwrap(),
            "accepted" => accepted = parts[1].parse().unwrap(),
            "chi2" => chi2 = parts[1].parse().unwrap(),
            "max_step" => max_step = parts[1].parse().unwrap(),
            "seed" => seed = parts[1].parse().unwrap(),
            "natoms" | "box" => {
                if parts[0] == "box" {
                    box_lengths[0] = parts[1].parse().unwrap();
                    box_lengths[1] = parts[2].parse().unwrap();
                    box_lengths[2] = parts[3].parse().unwrap();
                }
            }
            _ => {
                // Atom line: species type_id x y z
                if parts.len() >= 5 {
                    let name = parts[0].to_string();
                    let type_id: usize = parts[1].parse().unwrap();
                    let x: f64 = parts[2].parse().unwrap();
                    let y: f64 = parts[3].parse().unwrap();
                    let z: f64 = parts[4].parse().unwrap();
                    atoms.push(Atom {
                        position: [x, y, z],
                        species: name,
                        type_id,
                    });
                }
            }
        }
    }

    let mut composition: HashMap<String, usize> = HashMap::new();
    for atom in &atoms {
        *composition.entry(atom.species.clone()).or_insert(0) += 1;
    }

    let config = Configuration {
        atoms,
        box_lengths,
        species: species.to_vec(),
        composition,
    };

    let state = RmcState {
        move_count,
        accepted,
        chi2,
        max_step,
        seed,
        partial_sq: None,
        total_sq: None,
    };

    Ok((state, config))
}
