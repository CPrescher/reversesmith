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

/// Read a VASP POSCAR/CONTCAR file into a Configuration.
///
/// Supports VASP 5+ format with species names on line 6:
///   comment
///   scale_factor
///   a1x a1y a1z
///   a2x a2y a2z
///   a3x a3y a3z
///   Species1 Species2 ...
///   N1 N2 ...
///   Direct|Cartesian
///   x y z
///   ...
pub fn read_poscar(path: &Path) -> io::Result<Configuration> {
    let content = fs::read_to_string(path)?;
    let lines: Vec<&str> = content.lines().collect();

    if lines.len() < 8 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "POSCAR file too short",
        ));
    }

    // Line 1: comment (ignored)
    // Line 2: scale factor
    let scale: f64 = lines[1].trim().parse().map_err(|e| {
        io::Error::new(io::ErrorKind::InvalidData, format!("Bad scale factor: {}", e))
    })?;

    // Lines 3-5: lattice vectors (only orthorhombic supported)
    let mut lattice = [[0.0f64; 3]; 3];
    for i in 0..3 {
        let parts: Vec<f64> = lines[2 + i]
            .split_whitespace()
            .filter_map(|s| s.parse().ok())
            .collect();
        if parts.len() < 3 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Bad lattice vector on line {}", 3 + i),
            ));
        }
        lattice[i] = [parts[0] * scale, parts[1] * scale, parts[2] * scale];
    }

    let box_lengths = [lattice[0][0], lattice[1][1], lattice[2][2]];

    // Check for significant off-diagonal elements (non-orthorhombic)
    let off_diag = lattice[0][1].abs()
        + lattice[0][2].abs()
        + lattice[1][0].abs()
        + lattice[1][2].abs()
        + lattice[2][0].abs()
        + lattice[2][1].abs();
    if off_diag > 1e-6 * (box_lengths[0] + box_lengths[1] + box_lengths[2]) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Non-orthorhombic cells not supported (off-diagonal lattice elements too large)",
        ));
    }

    // Line 6: species names (VASP 5+)
    // Line 7: counts per species
    // Detect whether line 6 is species names or counts
    let line6_parts: Vec<&str> = lines[5].split_whitespace().collect();
    let line6_is_numeric = line6_parts.iter().all(|s| s.parse::<usize>().is_ok());

    let (species_names, counts, coord_start_line) = if line6_is_numeric {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "VASP 4 format (no species names on line 6) not supported; use VASP 5+ format",
        ));
    } else {
        let species: Vec<String> = line6_parts.iter().map(|s| s.to_string()).collect();
        let counts: Vec<usize> = lines[6]
            .split_whitespace()
            .map(|s| {
                s.parse().map_err(|e| {
                    io::Error::new(io::ErrorKind::InvalidData, format!("Bad atom count: {}", e))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        if species.len() != counts.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Species count ({}) != atom count entries ({})",
                    species.len(),
                    counts.len()
                ),
            ));
        }
        // Line 8: Selective dynamics (optional) or Direct/Cartesian
        let mut coord_line = 7;
        if lines[coord_line]
            .trim()
            .to_lowercase()
            .starts_with('s')
        {
            coord_line = 8; // skip "Selective dynamics"
        }
        (species, counts, coord_line)
    };

    let n_atoms: usize = counts.iter().sum();
    let is_direct = lines[coord_start_line]
        .trim()
        .to_lowercase()
        .starts_with('d');

    let mut atoms = Vec::with_capacity(n_atoms);
    let mut atom_line = coord_start_line + 1;

    for (type_id, (species_name, &count)) in species_names.iter().zip(counts.iter()).enumerate() {
        for _ in 0..count {
            if atom_line >= lines.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "POSCAR: fewer coordinate lines than expected",
                ));
            }
            let parts: Vec<f64> = lines[atom_line]
                .split_whitespace()
                .take(3)
                .filter_map(|s| s.parse().ok())
                .collect();
            if parts.len() < 3 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Bad coordinate on line {}", atom_line + 1),
                ));
            }

            let position = if is_direct {
                // Fractional → Cartesian (orthorhombic: simply scale)
                [
                    parts[0] * box_lengths[0],
                    parts[1] * box_lengths[1],
                    parts[2] * box_lengths[2],
                ]
            } else {
                [parts[0], parts[1], parts[2]]
            };

            atoms.push(Atom {
                position,
                species: species_name.clone(),
                type_id,
            });
            atom_line += 1;
        }
    }

    let mut composition: HashMap<String, usize> = HashMap::new();
    for atom in &atoms {
        *composition.entry(atom.species.clone()).or_insert(0) += 1;
    }

    let mut config = Configuration {
        atoms,
        box_lengths,
        species: species_names,
        composition,
    };
    config.wrap_pbc();
    Ok(config)
}

/// Write configuration as VASP POSCAR file (VASP 5+ format, Direct coordinates).
pub fn write_poscar(path: &Path, config: &Configuration, comment: &str) -> io::Result<()> {
    let mut file = fs::File::create(path)?;

    // Line 1: comment
    writeln!(file, "{}", comment)?;
    // Line 2: scale factor
    writeln!(file, "1.0")?;
    // Lines 3-5: lattice vectors (orthorhombic)
    writeln!(
        file,
        "  {:.10}  {:.10}  {:.10}",
        config.box_lengths[0], 0.0, 0.0
    )?;
    writeln!(
        file,
        "  {:.10}  {:.10}  {:.10}",
        0.0, config.box_lengths[1], 0.0
    )?;
    writeln!(
        file,
        "  {:.10}  {:.10}  {:.10}",
        0.0, 0.0, config.box_lengths[2]
    )?;

    // Line 6: species names (in species order)
    let species_line: Vec<&str> = config.species.iter().map(|s| s.as_str()).collect();
    writeln!(file, "  {}", species_line.join("  "))?;

    // Line 7: counts per species (in species order)
    let counts: Vec<String> = config
        .species
        .iter()
        .map(|s| config.composition.get(s).unwrap_or(&0).to_string())
        .collect();
    writeln!(file, "  {}", counts.join("  "))?;

    // Line 8: coordinate type
    writeln!(file, "Direct")?;

    // Atom coordinates in fractional, grouped by species
    for (type_id, _species_name) in config.species.iter().enumerate() {
        for atom in &config.atoms {
            if atom.type_id == type_id {
                writeln!(
                    file,
                    "  {:.10}  {:.10}  {:.10}",
                    atom.position[0] / config.box_lengths[0],
                    atom.position[1] / config.box_lengths[1],
                    atom.position[2] / config.box_lengths[2],
                )?;
            }
        }
    }

    Ok(())
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn poscar_roundtrip() {
        let dir = std::env::temp_dir().join("rsmith_poscar_test");
        std::fs::create_dir_all(&dir).unwrap();

        // Write a test POSCAR
        let poscar_in = dir.join("POSCAR");
        std::fs::write(
            &poscar_in,
            "CaSiO3 test\n\
             1.0\n\
               20.0  0.0  0.0\n\
               0.0  20.0  0.0\n\
               0.0  0.0  20.0\n\
             Ca  Si  O\n\
             1  1  3\n\
             Direct\n\
               0.1  0.2  0.3\n\
               0.4  0.5  0.6\n\
               0.15  0.25  0.35\n\
               0.55  0.65  0.75\n\
               0.85  0.95  0.05\n",
        )
        .unwrap();

        // Read
        let config = read_poscar(&poscar_in).unwrap();
        assert_eq!(config.atoms.len(), 5);
        assert_eq!(config.species, vec!["Ca", "Si", "O"]);
        assert_eq!(*config.composition.get("Ca").unwrap(), 1);
        assert_eq!(*config.composition.get("Si").unwrap(), 1);
        assert_eq!(*config.composition.get("O").unwrap(), 3);

        // Check first Ca atom: fractional (0.1, 0.2, 0.3) * 20 = (2, 4, 6)
        assert!((config.atoms[0].position[0] - 2.0).abs() < 1e-8);
        assert!((config.atoms[0].position[1] - 4.0).abs() < 1e-8);
        assert!((config.atoms[0].position[2] - 6.0).abs() < 1e-8);
        assert_eq!(config.atoms[0].species, "Ca");

        // Write back
        let poscar_out = dir.join("POSCAR_out");
        write_poscar(&poscar_out, &config, "roundtrip test").unwrap();

        // Read back and compare
        let config2 = read_poscar(&poscar_out).unwrap();
        assert_eq!(config2.atoms.len(), config.atoms.len());
        for (a, b) in config.atoms.iter().zip(config2.atoms.iter()) {
            assert_eq!(a.species, b.species);
            for d in 0..3 {
                assert!(
                    (a.position[d] - b.position[d]).abs() < 1e-6,
                    "Position mismatch: {:?} vs {:?}",
                    a.position,
                    b.position
                );
            }
        }

        // Cleanup
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn poscar_cartesian() {
        let dir = std::env::temp_dir().join("rsmith_poscar_cart_test");
        std::fs::create_dir_all(&dir).unwrap();

        let poscar = dir.join("POSCAR");
        std::fs::write(
            &poscar,
            "test cartesian\n\
             1.0\n\
               10.0  0.0  0.0\n\
               0.0  10.0  0.0\n\
               0.0  0.0  10.0\n\
             Si\n\
             2\n\
             Cartesian\n\
               1.0  2.0  3.0\n\
               5.0  6.0  7.0\n",
        )
        .unwrap();

        let config = read_poscar(&poscar).unwrap();
        assert_eq!(config.atoms.len(), 2);
        assert!((config.atoms[0].position[0] - 1.0).abs() < 1e-8);
        assert!((config.atoms[1].position[2] - 7.0).abs() < 1e-8);

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn poscar_scale_factor() {
        let dir = std::env::temp_dir().join("rsmith_poscar_scale_test");
        std::fs::create_dir_all(&dir).unwrap();

        let poscar = dir.join("POSCAR");
        std::fs::write(
            &poscar,
            "test scale\n\
             2.0\n\
               5.0  0.0  0.0\n\
               0.0  5.0  0.0\n\
               0.0  0.0  5.0\n\
             O\n\
             1\n\
             Direct\n\
               0.5  0.5  0.5\n",
        )
        .unwrap();

        let config = read_poscar(&poscar).unwrap();
        // Box should be 2.0 * 5.0 = 10.0
        assert!((config.box_lengths[0] - 10.0).abs() < 1e-8);
        // Position: 0.5 * 10.0 = 5.0
        assert!((config.atoms[0].position[0] - 5.0).abs() < 1e-8);

        std::fs::remove_dir_all(&dir).ok();
    }
}
