use std::collections::HashMap;

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
        self.distance_vec(i, j).iter().map(|d| d * d).sum::<f64>().sqrt()
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
