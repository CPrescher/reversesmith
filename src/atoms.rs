use std::collections::HashMap;

/// Return the standard molar mass (g/mol) for a given element symbol.
/// IUPAC 2021 standard atomic weights.
pub fn molar_mass(element: &str) -> f64 {
    match element {
        "H" => 1.008,
        "He" => 4.0026,
        "Li" => 6.941,
        "Be" => 9.0122,
        "B" => 10.81,
        "C" => 12.011,
        "N" => 14.007,
        "O" => 15.999,
        "F" => 18.998,
        "Ne" => 20.180,
        "Na" => 22.990,
        "Mg" => 24.305,
        "Al" => 26.982,
        "Si" => 28.086,
        "P" => 30.974,
        "S" => 32.06,
        "Cl" => 35.45,
        "Ar" => 39.948,
        "K" => 39.098,
        "Ca" => 40.078,
        "Sc" => 44.956,
        "Ti" => 47.867,
        "V" => 50.942,
        "Cr" => 51.996,
        "Mn" => 54.938,
        "Fe" => 55.845,
        "Co" => 58.933,
        "Ni" => 58.693,
        "Cu" => 63.546,
        "Zn" => 65.38,
        "Ga" => 69.723,
        "Ge" => 72.63,
        "As" => 74.922,
        "Se" => 78.971,
        "Br" => 79.904,
        "Kr" => 83.798,
        "Rb" => 85.468,
        "Sr" => 87.62,
        "Y" => 88.906,
        "Zr" => 91.224,
        "Nb" => 92.906,
        "Mo" => 95.95,
        "Tc" => 98.0,
        "Ru" => 101.07,
        "Rh" => 102.906,
        "Pd" => 106.42,
        "Ag" => 107.868,
        "Cd" => 112.414,
        "In" => 114.818,
        "Sn" => 118.710,
        "Sb" => 121.760,
        "Te" => 127.60,
        "I" => 126.904,
        "Xe" => 131.293,
        "Cs" => 132.905,
        "Ba" => 137.327,
        "La" => 138.905,
        "Ce" => 140.116,
        "Pr" => 140.908,
        "Nd" => 144.242,
        "Pm" => 145.0,
        "Sm" => 150.36,
        "Eu" => 151.964,
        "Gd" => 157.25,
        "Tb" => 158.925,
        "Dy" => 162.500,
        "Ho" => 164.930,
        "Er" => 167.259,
        "Tm" => 168.934,
        "Yb" => 173.045,
        "Lu" => 174.967,
        "Hf" => 178.49,
        "Ta" => 180.948,
        "W" => 183.84,
        "Re" => 186.207,
        "Os" => 190.23,
        "Ir" => 192.217,
        "Pt" => 195.084,
        "Au" => 196.967,
        "Hg" => 200.592,
        "Tl" => 204.38,
        "Pb" => 207.2,
        "Bi" => 208.980,
        "Po" => 209.0,
        "At" => 210.0,
        "Rn" => 222.0,
        "Fr" => 223.0,
        "Ra" => 226.0,
        "Ac" => 227.0,
        "Th" => 232.038,
        "Pa" => 231.036,
        "U" => 238.029,
        "Np" => 237.0,
        "Pu" => 244.0,
        "Am" => 243.0,
        "Cm" => 247.0,
        "Bk" => 247.0,
        "Cf" => 251.0,
        _ => panic!(
            "molar_mass: unknown element '{}'. Supported: H-Cf (Z=1-98).",
            element
        ),
    }
}

#[derive(Debug, Clone)]
pub struct Atom {
    pub position: [f64; 3],
    pub species: String,
    pub type_id: usize,
}

#[derive(Debug, Clone)]
pub struct Configuration {
    pub atoms: Vec<Atom>,
    pub box_lengths: [f64; 3],
    pub species: Vec<String>,
    pub composition: HashMap<String, usize>,
}

impl Configuration {
    /// Wrap atom positions into the primary box [0, L).
    pub fn wrap_pbc(&mut self) {
        for atom in &mut self.atoms {
            for d in 0..3 {
                let l = self.box_lengths[d];
                atom.position[d] -= l * (atom.position[d] / l).floor();
            }
        }
    }

    /// Minimum-image distance between atoms i and j.
    pub fn distance(&self, i: usize, j: usize) -> f64 {
        self.distance_vec(i, j)
            .iter()
            .map(|d| d * d)
            .sum::<f64>()
            .sqrt()
    }

    /// Minimum-image displacement vector from atom i to atom j.
    pub fn distance_vec(&self, i: usize, j: usize) -> [f64; 3] {
        let pi = &self.atoms[i].position;
        let pj = &self.atoms[j].position;
        let mut dr = [0.0f64; 3];
        for d in 0..3 {
            let mut delta = pj[d] - pi[d];
            let l = self.box_lengths[d];
            delta -= l * (delta / l).round();
            dr[d] = delta;
        }
        dr
    }

    /// Minimum-image distance between a position and atom j.
    pub fn distance_from_pos(&self, pos: &[f64; 3], j: usize) -> f64 {
        let pj = &self.atoms[j].position;
        let mut r2 = 0.0f64;
        for d in 0..3 {
            let mut delta = pj[d] - pos[d];
            let l = self.box_lengths[d];
            delta -= l * (delta / l).round();
            r2 += delta * delta;
        }
        r2.sqrt()
    }

    /// Move atom i by displacement, wrapping into box.
    pub fn move_atom(&mut self, i: usize, displacement: [f64; 3]) {
        for d in 0..3 {
            self.atoms[i].position[d] += displacement[d];
            let l = self.box_lengths[d];
            self.atoms[i].position[d] -= l * (self.atoms[i].position[d] / l).floor();
        }
    }

    /// Box volume.
    pub fn volume(&self) -> f64 {
        self.box_lengths[0] * self.box_lengths[1] * self.box_lengths[2]
    }

    /// Total number density.
    pub fn number_density(&self) -> f64 {
        self.atoms.len() as f64 / self.volume()
    }

    /// Mass density in g/cm^3.
    pub fn mass_density(&self) -> f64 {
        let na: f64 = 6.02214076e23;
        // total mass in g/mol
        let total_mass: f64 = self
            .composition
            .iter()
            .map(|(el, &n)| n as f64 * molar_mass(el))
            .sum();
        // volume in A^3 -> cm^3: 1 A^3 = 1e-24 cm^3
        (total_mass / self.volume()) * (1e24 / na)
    }

    /// Number of distinct type pairs (upper triangle including diagonal).
    pub fn num_type_pairs(&self) -> usize {
        let n = self.species.len();
        n * (n + 1) / 2
    }

    /// Map (type_a, type_b) with a <= b to a linear index.
    pub fn pair_index(&self, a: usize, b: usize) -> usize {
        let (a, b) = if a <= b { (a, b) } else { (b, a) };
        let n = self.species.len();
        a * n - a * (a + 1) / 2 + b
    }

    /// Get concentration of species by type_id.
    pub fn concentration(&self, type_id: usize) -> f64 {
        let name = &self.species[type_id];
        *self.composition.get(name).unwrap_or(&0) as f64 / self.atoms.len() as f64
    }

    /// Count of atoms with given type_id.
    pub fn count_type(&self, type_id: usize) -> usize {
        let name = &self.species[type_id];
        *self.composition.get(name).unwrap_or(&0)
    }
}
