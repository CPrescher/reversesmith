//! Pair potential energy for hybrid RMC refinement.
//!
//! Supports Lennard-Jones 12-6, Buckingham, Pedone (Morse + r⁻¹²), Coulomb
//! DSF, and tabulated potentials.
//! All potentials are tabulated on a uniform grid (dr = 0.001 A) for fast lookup.
//! Multiple potential types can be combined additively for the same pair.

use std::f64::consts::PI;
use std::path::Path;

use crate::atoms::Configuration;
use crate::cells::{CellList, VerletNeighborList};
use crate::config::PotentialConfig;
use crate::energy::EnergyModel;
use crate::io;

/// Conversion factor: e^2 / (4πε₀) in eV·A units.
const COULOMB_CONST: f64 = 14.3997;

/// A single pair potential tabulated on a uniform grid.
///
/// All analytical potentials (Lennard-Jones, Buckingham, Pedone, Coulomb) are tabulated at
/// construction time for consistent, fast evaluation. The potential is shifted
/// so V(cutoff) = 0 to avoid energy discontinuities.
#[derive(Clone)]
pub struct PairPotential {
    pub pair_label: String,
    pub type_a: usize,
    pub type_b: usize,
    /// V(r) on uniform grid: `table[i] = V(i * dr)`. Shifted so V(cutoff) = 0.
    pub table: Vec<f64>,
    pub dr: f64,
    pub n_bins: usize,
    pub cutoff: f64,
}

/// Collection of pair potentials with O(1) type-pair lookup.
///
/// Built from a `[potential]` TOML config section. The `weight` field controls
/// the energy contribution to the RMC cost function:
/// `cost = chi2 + weight * E_total`.
#[derive(Clone)]
pub struct PotentialSet {
    pub potentials: Vec<PairPotential>,
    pub weight: f64,
    pub cutoff: f64,
    /// Lookup: potential_index[type_a * n_types + type_b] = index into potentials, or usize::MAX
    pub potential_index: Vec<usize>,
    pub n_types: usize,
}

impl PairPotential {
    /// Tabulate Lennard-Jones potential
    /// V(r) = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6].
    pub fn from_lennard_jones(
        pair_label: String,
        type_a: usize,
        type_b: usize,
        epsilon: f64,
        sigma: f64,
        cutoff: f64,
        dr: f64,
    ) -> Self {
        let n_bins = (cutoff / dr).ceil() as usize + 1;
        let mut table = vec![0.0; n_bins];

        let v_cut = lennard_jones_eval(cutoff, epsilon, sigma);
        for (i, value) in table.iter_mut().enumerate() {
            let r = (i as f64 * dr).max(0.5);
            *value = lennard_jones_eval(r, epsilon, sigma) - v_cut;
        }
        if let Some(last) = table.last_mut() {
            *last = 0.0;
        }

        PairPotential {
            pair_label,
            type_a,
            type_b,
            table,
            dr,
            n_bins,
            cutoff,
        }
    }

    /// Tabulate Buckingham potential V(r) = A*exp(-r/rho) - C/r^6.
    pub fn from_buckingham(
        pair_label: String,
        type_a: usize,
        type_b: usize,
        a_param: f64,
        rho: f64,
        c_param: f64,
        cutoff: f64,
        dr: f64,
    ) -> Self {
        let n_bins = (cutoff / dr).ceil() as usize + 1;
        let mut table = vec![0.0; n_bins];

        let v_cut = buckingham_eval(cutoff, a_param, rho, c_param);
        for i in 0..n_bins {
            let r = (i as f64 * dr).max(0.5); // cap at 0.5 A to avoid divergence
            table[i] = buckingham_eval(r, a_param, rho, c_param) - v_cut;
        }
        if n_bins > 0 {
            table[n_bins - 1] = 0.0;
        }

        PairPotential {
            pair_label,
            type_a,
            type_b,
            table,
            dr,
            n_bins,
            cutoff,
        }
    }

    /// Tabulate Pedone potential V(r) = D0*[1 - exp(-α(r-r0))]² - D0 + C0/r¹².
    pub fn from_pedone(
        pair_label: String,
        type_a: usize,
        type_b: usize,
        d0: f64,
        alpha: f64,
        r0: f64,
        c0: f64,
        cutoff: f64,
        dr: f64,
    ) -> Self {
        let n_bins = (cutoff / dr).ceil() as usize + 1;
        let mut table = vec![0.0; n_bins];

        let v_cut = pedone_eval(cutoff, d0, alpha, r0, c0);
        for i in 0..n_bins {
            let r = (i as f64 * dr).max(0.5);
            table[i] = pedone_eval(r, d0, alpha, r0, c0) - v_cut;
        }
        if n_bins > 0 {
            table[n_bins - 1] = 0.0;
        }

        PairPotential {
            pair_label,
            type_a,
            type_b,
            table,
            dr,
            n_bins,
            cutoff,
        }
    }

    /// Tabulate Coulomb DSF interaction for a specific pair: qi*qj * V_dsf(r).
    /// DSF (Wolf 1999 / Fennell & Gezelter 2006):
    ///   V(r) = K * qi*qj * [erfc(α*r)/r - erfc(α*rc)/rc
    ///          + (erfc(α*rc)/rc² + 2α/√π * exp(-α²rc²)/rc) * (r - rc)]
    pub fn from_coulomb_dsf(
        pair_label: String,
        type_a: usize,
        type_b: usize,
        qi: f64,
        qj: f64,
        alpha_dsf: f64,
        cutoff: f64,
        dr: f64,
    ) -> Self {
        let n_bins = (cutoff / dr).ceil() as usize + 1;
        let mut table = vec![0.0; n_bins];

        // DSF constants at cutoff
        let erfc_rc = erfc_approx(alpha_dsf * cutoff);
        let shift = erfc_rc / cutoff;
        let dshift = erfc_rc / (cutoff * cutoff)
            + 2.0 * alpha_dsf / PI.sqrt() * (-alpha_dsf * alpha_dsf * cutoff * cutoff).exp()
                / cutoff;

        let qq = COULOMB_CONST * qi * qj;

        for i in 0..n_bins {
            let r = (i as f64 * dr).max(0.3); // cap short range
            if r >= cutoff {
                table[i] = 0.0;
            } else {
                let erfc_r = erfc_approx(alpha_dsf * r);
                table[i] = qq * (erfc_r / r - shift + dshift * (r - cutoff));
            }
        }
        // The DSF form is self-shifting: V(rc) = 0 by construction.
        // But enforce exactly due to floating point:
        if n_bins > 0 {
            table[n_bins - 1] = 0.0;
        }

        PairPotential {
            pair_label,
            type_a,
            type_b,
            table,
            dr,
            n_bins,
            cutoff,
        }
    }

    /// Read a tabulated potential from a two-column file (r in A, V in eV),
    /// interpolate onto uniform grid, shift so V(cutoff) = 0.
    pub fn from_table(
        pair_label: String,
        type_a: usize,
        type_b: usize,
        path: &Path,
        cutoff: f64,
        dr: f64,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let (r_data, v_data) = io::read_sq_data(path)?;
        if r_data.len() < 2 {
            return Err(
                format!("Potential table {} has fewer than 2 points", path.display()).into(),
            );
        }

        let n_bins = (cutoff / dr).ceil() as usize + 1;
        let mut table = vec![0.0; n_bins];

        let v_cut = interp(&r_data, &v_data, cutoff);
        for i in 0..n_bins {
            let r = i as f64 * dr;
            table[i] = interp(&r_data, &v_data, r) - v_cut;
        }
        if n_bins > 0 {
            table[n_bins - 1] = 0.0;
        }

        Ok(PairPotential {
            pair_label,
            type_a,
            type_b,
            table,
            dr,
            n_bins,
            cutoff,
        })
    }

    /// Create a `PairPotential` from a pre-built table, with cutoff shift applied.
    pub fn from_vec(
        pair_label: String,
        type_a: usize,
        type_b: usize,
        mut table: Vec<f64>,
        cutoff: f64,
        dr: f64,
    ) -> Self {
        let n_bins = table.len();
        // Shift so V(cutoff) = 0
        let v_cut = if n_bins > 0 { table[n_bins - 1] } else { 0.0 };
        for v in &mut table {
            *v -= v_cut;
        }
        if n_bins > 0 {
            table[n_bins - 1] = 0.0;
        }
        PairPotential {
            pair_label,
            type_a,
            type_b,
            table,
            dr,
            n_bins,
            cutoff,
        }
    }

    /// Add another potential's table onto this one (element-wise).
    pub fn add_table(&mut self, other: &PairPotential) {
        let n = self.n_bins.min(other.n_bins);
        for i in 0..n {
            self.table[i] += other.table[i];
        }
    }

    /// Look up V(r) by linear interpolation on the table.
    #[inline]
    pub fn evaluate(&self, r: f64) -> f64 {
        if r >= self.cutoff {
            return 0.0;
        }
        let x = r / self.dr;
        let i = x as usize;
        if i + 1 >= self.n_bins {
            return 0.0;
        }
        let t = x - i as f64;
        self.table[i] + t * (self.table[i + 1] - self.table[i])
    }
}

impl PotentialSet {
    /// Build from TOML config.
    pub fn from_config(
        cfg: &PotentialConfig,
        species: &[String],
        rdf_cutoff: f64,
        base_dir: &Path,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let weight = cfg.weight.unwrap_or(0.001);
        let cutoff = cfg.cutoff.unwrap_or(rdf_cutoff);
        let dr = 0.001;
        let n_types = species.len();

        if cutoff > rdf_cutoff {
            return Err(format!(
                "Potential cutoff ({:.1} A) exceeds rdf_cutoff ({:.1} A). \
                 Increase [sq] rdf_cutoff or reduce [potential] cutoff.",
                cutoff, rdf_cutoff
            )
            .into());
        }

        let parse_pair = |pair: &str| -> Result<(usize, usize), Box<dyn std::error::Error>> {
            let parts: Vec<&str> = pair.split('-').collect();
            if parts.len() != 2 {
                return Err(format!("Invalid pair format '{}', expected 'A-B'", pair).into());
            }
            let a = species
                .iter()
                .position(|s| s == parts[0])
                .ok_or_else(|| format!("Unknown species '{}' in pair '{}'", parts[0], pair))?;
            let b = species
                .iter()
                .position(|s| s == parts[1])
                .ok_or_else(|| format!("Unknown species '{}' in pair '{}'", parts[1], pair))?;
            Ok((a, b))
        };

        let mut potentials: Vec<PairPotential> = Vec::new();
        let mut defined_pairs: Vec<(usize, usize)> = Vec::new();

        // Helper: find existing potential for a pair, or create a new slot
        let find_or_push = |potentials: &mut Vec<PairPotential>,
                            defined_pairs: &mut Vec<(usize, usize)>,
                            pot: PairPotential,
                            a: usize,
                            b: usize| {
            let existing = defined_pairs
                .iter()
                .position(|&(pa, pb)| (pa == a && pb == b) || (pa == b && pb == a));
            if let Some(idx) = existing {
                // Add to existing (accumulate short-range + Coulomb)
                potentials[idx].add_table(&pot);
            } else {
                potentials.push(pot);
                defined_pairs.push((a, b));
            }
        };

        // 1. Lennard-Jones potentials
        if let Some(ref terms) = cfg.lennard_jones {
            for term in terms {
                let (a, b) = parse_pair(&term.pair)?;
                let pot = PairPotential::from_lennard_jones(
                    term.pair.clone(),
                    a,
                    b,
                    term.epsilon,
                    term.sigma,
                    cutoff,
                    dr,
                );
                find_or_push(&mut potentials, &mut defined_pairs, pot, a, b);
            }
        }

        // 2. Buckingham potentials
        if let Some(ref bucks) = cfg.buckingham {
            for buck in bucks {
                let (a, b) = parse_pair(&buck.pair)?;
                let pot = PairPotential::from_buckingham(
                    buck.pair.clone(),
                    a,
                    b,
                    buck.a_param,
                    buck.rho,
                    buck.c_param,
                    cutoff,
                    dr,
                );
                find_or_push(&mut potentials, &mut defined_pairs, pot, a, b);
            }
        }

        // 3. Pedone potentials
        if let Some(ref peds) = cfg.pedone {
            for ped in peds {
                let (a, b) = parse_pair(&ped.pair)?;
                let pot = PairPotential::from_pedone(
                    ped.pair.clone(),
                    a,
                    b,
                    ped.d0,
                    ped.alpha,
                    ped.r0,
                    ped.c0,
                    cutoff,
                    dr,
                );
                find_or_push(&mut potentials, &mut defined_pairs, pot, a, b);
            }
        }

        // 4. Coulomb DSF (applies to all pairs with nonzero qi*qj)
        if let Some(ref coul) = cfg.coulomb {
            for a in 0..n_types {
                for b in a..n_types {
                    let qi = match coul.charges.get(&species[a]) {
                        Some(&q) => q,
                        None => continue,
                    };
                    let qj = match coul.charges.get(&species[b]) {
                        Some(&q) => q,
                        None => continue,
                    };
                    if (qi * qj).abs() < 1e-15 {
                        continue;
                    }

                    let label = format!("{}-{}", species[a], species[b]);
                    let pot = PairPotential::from_coulomb_dsf(
                        label, a, b, qi, qj, coul.alpha, cutoff, dr,
                    );
                    find_or_push(&mut potentials, &mut defined_pairs, pot, a, b);
                }
            }
        }

        // 5. Tabulated potentials (override everything on conflict)
        if let Some(ref tabs) = cfg.tabulated {
            for tab in tabs {
                let (a, b) = parse_pair(&tab.pair)?;
                let path = if Path::new(&tab.file).is_absolute() {
                    Path::new(&tab.file).to_path_buf()
                } else {
                    base_dir.join(&tab.file)
                };
                let pot = PairPotential::from_table(tab.pair.clone(), a, b, &path, cutoff, dr)?;

                // Tabulated replaces (not adds to) any existing entry
                let existing = defined_pairs
                    .iter()
                    .position(|&(pa, pb)| (pa == a && pb == b) || (pa == b && pb == a));
                if let Some(idx) = existing {
                    potentials[idx] = pot;
                } else {
                    potentials.push(pot);
                    defined_pairs.push((a, b));
                }
            }
        }

        // Build type-pair lookup table
        let mut potential_index = vec![usize::MAX; n_types * n_types];
        for (i, pot) in potentials.iter().enumerate() {
            potential_index[pot.type_a * n_types + pot.type_b] = i;
            potential_index[pot.type_b * n_types + pot.type_a] = i;
        }

        Ok(PotentialSet {
            potentials,
            weight,
            cutoff,
            potential_index,
            n_types,
        })
    }

    /// Insert or additively combine a `PairPotential` into the set.
    ///
    /// If a potential already exists for this pair, the tables are summed.
    /// Otherwise a new entry is created and the lookup table updated.
    pub fn add_potential(&mut self, pot: PairPotential) {
        let a = pot.type_a;
        let b = pot.type_b;
        let existing = self.potential_index[a * self.n_types + b];
        if existing != usize::MAX {
            self.potentials[existing].add_table(&pot);
        } else {
            let idx = self.potentials.len();
            self.potential_index[a * self.n_types + b] = idx;
            self.potential_index[b * self.n_types + a] = idx;
            self.potentials.push(pot);
        }
    }

    /// Compute the pair potential energy contribution of a single atom at a given position.
    /// Sums V(r_ij) for all j neighbors within cutoff using the cell list.
    #[inline]
    pub fn energy_of_atom(
        &self,
        config: &Configuration,
        atom_idx: usize,
        pos: &[f64; 3],
        cell_list: &CellList,
        pos_cell: usize,
    ) -> f64 {
        let ti = config.atoms[atom_idx].type_id;
        let cutoff2 = self.cutoff * self.cutoff;
        let n_types = self.n_types;

        let mut energy = 0.0f64;
        let neighbor_cells = cell_list.neighbor_cells(pos_cell);
        for &nc in &neighbor_cells {
            for j in cell_list.atoms_in_cell(nc) {
                if j == atom_idx {
                    continue;
                }

                let tj = config.atoms[j].type_id;
                let pot_idx = self.potential_index[ti * n_types + tj];
                if pot_idx == usize::MAX {
                    continue;
                }

                let r2 =
                    periodic_distance_squared(&config.atoms[j].position, pos, &config.box_lengths);

                if r2 < cutoff2 {
                    let r = r2.sqrt();
                    energy += self.potentials[pot_idx].evaluate(r);
                }
            }
        }
        energy
    }

    /// Compute energy delta for moving an atom from old_pos to new_pos.
    /// Single pass over neighbors when old and new cells share the same neighbor set.
    pub fn energy_delta_atom(
        &self,
        config: &Configuration,
        atom_idx: usize,
        old_pos: &[f64; 3],
        new_pos: &[f64; 3],
        cell_list: &CellList,
        old_cell: usize,
        new_cell: usize,
    ) -> f64 {
        // If cells differ, the neighbor sets may differ — fall back to two calls
        if old_cell != new_cell {
            let old_e = self.energy_of_atom(config, atom_idx, old_pos, cell_list, old_cell);
            let new_e = self.energy_of_atom(config, atom_idx, new_pos, cell_list, new_cell);
            return new_e - old_e;
        }

        let ti = config.atoms[atom_idx].type_id;
        let cutoff2 = self.cutoff * self.cutoff;
        let n_types = self.n_types;

        let mut delta = 0.0f64;
        let neighbor_cells = cell_list.neighbor_cells(old_cell);
        for &nc in &neighbor_cells {
            for j in cell_list.atoms_in_cell(nc) {
                if j == atom_idx {
                    continue;
                }

                let tj = config.atoms[j].type_id;
                let pot_idx = self.potential_index[ti * n_types + tj];
                if pot_idx == usize::MAX {
                    continue;
                }

                let position = &config.atoms[j].position;
                let r2_old = periodic_distance_squared(position, old_pos, &config.box_lengths);
                let r2_new = periodic_distance_squared(position, new_pos, &config.box_lengths);

                let pot = &self.potentials[pot_idx];
                let e_new = if r2_new < cutoff2 {
                    pot.evaluate(r2_new.sqrt())
                } else {
                    0.0
                };
                let e_old = if r2_old < cutoff2 {
                    pot.evaluate(r2_old.sqrt())
                } else {
                    0.0
                };
                delta += e_new - e_old;
            }
        }
        delta
    }

    /// Compute the energy delta from a pre-built cutoff-plus-skin neighbor
    /// slice. This is the hot path used by pure EPSR Monte Carlo.
    #[inline]
    pub fn energy_delta_atom_neighbors(
        &self,
        config: &Configuration,
        atom_idx: usize,
        old_pos: &[f64; 3],
        new_pos: &[f64; 3],
        neighbors: &[usize],
    ) -> f64 {
        let ti = config.atoms[atom_idx].type_id;
        let cutoff2 = self.cutoff * self.cutoff;
        let n_types = self.n_types;
        let mut delta = 0.0;

        for &j in neighbors {
            let tj = config.atoms[j].type_id;
            let pot_idx = self.potential_index[ti * n_types + tj];
            if pot_idx == usize::MAX {
                continue;
            }
            let position = &config.atoms[j].position;
            let r2_old = periodic_distance_squared(position, old_pos, &config.box_lengths);
            let r2_new = periodic_distance_squared(position, new_pos, &config.box_lengths);
            let potential = &self.potentials[pot_idx];
            let new_energy = if r2_new < cutoff2 {
                potential.evaluate(r2_new.sqrt())
            } else {
                0.0
            };
            let old_energy = if r2_old < cutoff2 {
                potential.evaluate(r2_old.sqrt())
            } else {
                0.0
            };
            delta += new_energy - old_energy;
        }
        delta
    }

    /// Total pair energy using a directed Verlet list. Each physical pair is
    /// counted once by retaining only neighbors with `j > i`.
    pub fn total_energy_neighbors(
        &self,
        config: &Configuration,
        neighbor_list: &VerletNeighborList,
    ) -> f64 {
        debug_assert!((neighbor_list.cutoff() - self.cutoff).abs() < 1.0e-12);
        let cutoff2 = self.cutoff * self.cutoff;
        let n_types = self.n_types;
        let mut energy = 0.0;

        for i in 0..config.atoms.len() {
            let type_i = config.atoms[i].type_id;
            let position_i = &config.atoms[i].position;
            for &j in neighbor_list.neighbors(i) {
                if j <= i {
                    continue;
                }
                let type_j = config.atoms[j].type_id;
                let potential_index = self.potential_index[type_i * n_types + type_j];
                if potential_index == usize::MAX {
                    continue;
                }
                let r2 = periodic_distance_squared(
                    position_i,
                    &config.atoms[j].position,
                    &config.box_lengths,
                );
                if r2 < cutoff2 {
                    energy += self.potentials[potential_index].evaluate(r2.sqrt());
                }
            }
        }
        energy
    }

    /// Compute total pair potential energy of the configuration.
    /// Each pair is counted once (i < j).
    pub fn total_energy(&self, config: &Configuration, cell_list: &CellList) -> f64 {
        let n_atoms = config.atoms.len();
        let cutoff2 = self.cutoff * self.cutoff;
        let n_types = self.n_types;

        let mut energy = 0.0f64;
        for i in 0..n_atoms {
            let ti = config.atoms[i].type_id;
            let pi = &config.atoms[i].position;
            let ci = cell_list.cell_of[i];
            let neighbor_cells = cell_list.neighbor_cells(ci);
            for &nc in &neighbor_cells {
                for j in cell_list.atoms_in_cell(nc) {
                    if j <= i {
                        continue;
                    }

                    let tj = config.atoms[j].type_id;
                    let pot_idx = self.potential_index[ti * n_types + tj];
                    if pot_idx == usize::MAX {
                        continue;
                    }

                    let r2 = periodic_distance_squared(
                        &config.atoms[j].position,
                        pi,
                        &config.box_lengths,
                    );

                    if r2 < cutoff2 {
                        let r = r2.sqrt();
                        energy += self.potentials[pot_idx].evaluate(r);
                    }
                }
            }
        }
        energy
    }
}

impl EnergyModel for PotentialSet {
    fn label(&self) -> &str {
        "pair potentials"
    }

    fn weight(&self) -> f64 {
        self.weight
    }

    fn cutoff(&self) -> f64 {
        self.cutoff
    }

    fn total_energy(&mut self, config: &Configuration, cell_list: &CellList) -> f64 {
        PotentialSet::total_energy(self, config, cell_list)
    }

    fn energy_delta_atom(
        &mut self,
        config: &Configuration,
        atom_idx: usize,
        old_pos: &[f64; 3],
        new_pos: &[f64; 3],
        cell_list: &CellList,
        old_cell: usize,
        new_cell: usize,
    ) -> f64 {
        PotentialSet::energy_delta_atom(
            self, config, atom_idx, old_pos, new_pos, cell_list, old_cell, new_cell,
        )
    }
}

// --- Helper functions ---

#[inline]
fn buckingham_eval(r: f64, a: f64, rho: f64, c: f64) -> f64 {
    a * (-r / rho).exp() - c / r.powi(6)
}

#[inline]
fn lennard_jones_eval(r: f64, epsilon: f64, sigma: f64) -> f64 {
    let sr6 = (sigma / r).powi(6);
    4.0 * epsilon * (sr6 * sr6 - sr6)
}

#[inline]
fn periodic_distance_squared(a: &[f64; 3], b: &[f64; 3], box_lengths: &[f64; 3]) -> f64 {
    let mut distance2 = 0.0;
    for axis in 0..3 {
        let mut delta = a[axis] - b[axis];
        let half = 0.5 * box_lengths[axis];
        if delta > half {
            delta -= box_lengths[axis];
        } else if delta < -half {
            delta += box_lengths[axis];
        }
        distance2 += delta * delta;
    }
    distance2
}

#[inline]
fn pedone_eval(r: f64, d0: f64, alpha: f64, r0: f64, c0: f64) -> f64 {
    let morse = 1.0 - (-alpha * (r - r0)).exp();
    d0 * morse * morse - d0 + c0 / r.powi(12)
}

/// Approximation of erfc(x) using Abramowitz & Stegun formula 7.1.26.
/// Max error ~1.5e-7, sufficient for potential tabulation.
#[inline]
fn erfc_approx(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.3275911 * x.abs());
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let result = poly * (-x * x).exp();
    if x >= 0.0 {
        result
    } else {
        2.0 - result
    }
}

/// Linear interpolation helper for tabulated data.
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

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::{lennard_jones_eval, periodic_distance_squared, PairPotential, PotentialSet};
    use crate::atoms::{Atom, Configuration};
    use crate::cells::VerletNeighborList;

    #[test]
    fn lennard_jones_has_the_expected_zero_and_minimum_before_cutoff_shift() {
        let epsilon = 0.125;
        let sigma = 2.4;
        let cutoff = 8.0;
        let potential = PairPotential::from_lennard_jones(
            "Si-O".to_string(),
            0,
            1,
            epsilon,
            sigma,
            cutoff,
            0.001,
        );
        let shift = lennard_jones_eval(cutoff, epsilon, sigma);

        assert!((potential.evaluate(sigma) + shift).abs() < 1.0e-12);

        let minimum = 2.0_f64.powf(1.0 / 6.0) * sigma;
        assert!((potential.evaluate(minimum) - (-epsilon - shift)).abs() < 2.0e-7);
        assert_eq!(potential.evaluate(cutoff), 0.0);
    }

    #[test]
    fn lennard_jones_matches_lammps_lj_cut_with_energy_shift() {
        // Oracle: LAMMPS 22 Jul 2025 Update 4, units metal,
        // pair_style lj/cut 8.0; pair_coeff 1 1 0.010 3.0;
        // pair_modify shift yes; two atoms separated by 3.0 A.
        let potential =
            PairPotential::from_lennard_jones("A-A".to_string(), 0, 0, 0.010, 3.0, 8.0, 0.001);
        let lammps_energy = 0.000_110_927_232_890_385_65;
        assert!((potential.evaluate(3.0) - lammps_energy).abs() < 1.0e-15);
    }

    fn two_species_fixture() -> (Configuration, PotentialSet) {
        let box_lengths = [13.0, 14.0, 15.0];
        let mut atoms = Vec::new();
        for z in 0..3 {
            for y in 0..4 {
                for x in 0..4 {
                    let type_id = (x + 2 * y + z) % 2;
                    atoms.push(Atom {
                        position: [
                            0.35 + x as f64 * 2.9,
                            0.45 + y as f64 * 3.1,
                            0.55 + z as f64 * 4.0,
                        ],
                        species: if type_id == 0 { "A" } else { "B" }.to_string(),
                        type_id,
                    });
                }
            }
        }
        let mut composition = HashMap::new();
        composition.insert(
            "A".to_string(),
            atoms.iter().filter(|atom| atom.type_id == 0).count(),
        );
        composition.insert(
            "B".to_string(),
            atoms.iter().filter(|atom| atom.type_id == 1).count(),
        );
        let config = Configuration {
            atoms,
            box_lengths,
            species: vec!["A".to_string(), "B".to_string()],
            composition,
        };

        let cutoff = 5.5;
        let potentials = vec![
            PairPotential::from_lennard_jones("A-A".to_string(), 0, 0, 0.010, 2.2, cutoff, 0.001),
            PairPotential::from_lennard_jones("A-B".to_string(), 0, 1, 0.013, 2.1, cutoff, 0.001),
            PairPotential::from_lennard_jones("B-B".to_string(), 1, 1, 0.008, 2.0, cutoff, 0.001),
        ];
        let potential_set = PotentialSet {
            potentials,
            weight: 1.0,
            cutoff,
            potential_index: vec![0, 1, 1, 2],
            n_types: 2,
        };
        (config, potential_set)
    }

    fn brute_atom_energy(
        potentials: &PotentialSet,
        config: &Configuration,
        atom_idx: usize,
        position: &[f64; 3],
    ) -> f64 {
        let type_i = config.atoms[atom_idx].type_id;
        let cutoff2 = potentials.cutoff * potentials.cutoff;
        let mut energy = 0.0;
        for (j, atom) in config.atoms.iter().enumerate() {
            if j == atom_idx {
                continue;
            }
            let potential_index =
                potentials.potential_index[type_i * potentials.n_types + atom.type_id];
            let r2 = periodic_distance_squared(position, &atom.position, &config.box_lengths);
            if r2 < cutoff2 {
                energy += potentials.potentials[potential_index].evaluate(r2.sqrt());
            }
        }
        energy
    }

    #[test]
    fn verlet_energy_deltas_match_brute_force_through_moves_and_rebuilds() {
        let (mut config, potentials) = two_species_fixture();
        let mut positions: Vec<[f64; 3]> = config.atoms.iter().map(|atom| atom.position).collect();
        let mut neighbors =
            VerletNeighborList::new(&positions, &config.box_lengths, potentials.cutoff, 1.0);

        let brute_total: f64 = (0..config.atoms.len())
            .map(|i| brute_atom_energy(&potentials, &config, i, &config.atoms[i].position))
            .sum::<f64>()
            * 0.5;
        let listed_total = potentials.total_energy_neighbors(&config, &neighbors);
        assert!((listed_total - brute_total).abs() < 1.0e-12 * (1.0 + brute_total.abs()));

        for trial in 0..240 {
            let atom_idx = trial * 17 % config.atoms.len();
            let old_pos = config.atoms[atom_idx].position;
            let displacement = [
                ((trial * 13 % 19) as f64 - 9.0) * 0.055,
                ((trial * 7 % 17) as f64 - 8.0) * 0.047,
                ((trial * 11 % 23) as f64 - 11.0) * 0.039,
            ];
            let mut new_pos = [0.0; 3];
            for axis in 0..3 {
                new_pos[axis] =
                    (old_pos[axis] + displacement[axis]).rem_euclid(config.box_lengths[axis]);
            }

            neighbors.prepare_trial(&positions, atom_idx, &new_pos);
            let listed_delta = potentials.energy_delta_atom_neighbors(
                &config,
                atom_idx,
                &old_pos,
                &new_pos,
                neighbors.neighbors(atom_idx),
            );
            let brute_delta = brute_atom_energy(&potentials, &config, atom_idx, &new_pos)
                - brute_atom_energy(&potentials, &config, atom_idx, &old_pos);
            assert!(
                (listed_delta - brute_delta).abs() < 1.0e-12 * (1.0 + brute_delta.abs()),
                "trial {trial}: listed={listed_delta}, brute={brute_delta}"
            );

            config.atoms[atom_idx].position = new_pos;
            positions[atom_idx] = new_pos;
            neighbors.accept_move(atom_idx, &new_pos);
        }
        assert!(neighbors.rebuilds() > 1, "test did not exercise rebuilds");
    }
}
