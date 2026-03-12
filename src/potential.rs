//! Pair potential energy for hybrid RMC refinement.
//!
//! Supports Buckingham, Pedone (Morse + r⁻¹²), Coulomb DSF, and tabulated potentials.
//! All potentials are tabulated on a uniform grid (dr = 0.001 A) for fast lookup.
//! Multiple potential types can be combined additively for the same pair.

use std::f64::consts::PI;
use std::path::Path;

use crate::atoms::Configuration;
use crate::cells::CellList;
use crate::config::PotentialConfig;
use crate::io;

/// Conversion factor: e^2 / (4πε₀) in eV·A units.
const COULOMB_CONST: f64 = 14.3997;

/// A single pair potential tabulated on a uniform grid.
///
/// All analytical potentials (Buckingham, Pedone, Coulomb) are tabulated at
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

        // 1. Buckingham potentials
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

        // 2. Pedone potentials
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

        // 3. Coulomb DSF (applies to all pairs with nonzero qi*qj)
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

        // 4. Tabulated potentials (override everything on conflict)
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
        let box_lengths = &config.box_lengths;
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

                let pj = &config.atoms[j].position;
                let mut r2 = 0.0f64;
                for d in 0..3 {
                    let mut delta = pj[d] - pos[d];
                    let l = box_lengths[d];
                    delta -= l * (delta / l).round();
                    r2 += delta * delta;
                }

                if r2 < cutoff2 {
                    let r = r2.sqrt();
                    energy += self.potentials[pot_idx].evaluate(r);
                }
            }
        }
        energy
    }

    /// Compute total pair potential energy of the configuration.
    /// Each pair is counted once (i < j).
    pub fn total_energy(&self, config: &Configuration, cell_list: &CellList) -> f64 {
        let n_atoms = config.atoms.len();
        let box_lengths = &config.box_lengths;
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

                    let pj = &config.atoms[j].position;
                    let mut r2 = 0.0f64;
                    for d in 0..3 {
                        let mut delta = pj[d] - pi[d];
                        let l = box_lengths[d];
                        delta -= l * (delta / l).round();
                        r2 += delta * delta;
                    }

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

// --- Helper functions ---

#[inline]
fn buckingham_eval(r: f64, a: f64, rho: f64, c: f64) -> f64 {
    a * (-r / rho).exp() - c / r.powi(6)
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
