/// Cell list for O(1)-per-neighbor spatial lookups with periodic boundary conditions.
pub struct CellList {
    /// Number of cells per dimension.
    pub nc: [usize; 3],
    /// Cell size per dimension (>= cutoff).
    pub cell_size: [f64; 3],
    /// Box lengths.
    pub box_lengths: [f64; 3],
    /// `head[cell_idx]` = first atom in cell, or usize::MAX if empty.
    pub head: Vec<usize>,
    /// `next[atom_idx]` = next atom in same cell, or usize::MAX if end.
    pub next: Vec<usize>,
    /// `cell_of[atom_idx]` = cell index for this atom.
    pub cell_of: Vec<usize>,
    /// Number of cells total.
    pub n_cells: usize,
}

/// Contiguous per-atom Verlet neighbor list for repeated short-range trials.
///
/// Pairs are included out to `cutoff + skin`. The list remains valid while
/// the sum of the two largest atom displacements from the reference
/// configuration does not exceed `skin`.
pub struct VerletNeighborList {
    cutoff: f64,
    skin: f64,
    list_radius: f64,
    box_lengths: [f64; 3],
    reference_positions: Vec<[f64; 3]>,
    offsets: Vec<usize>,
    neighbors: Vec<usize>,
    max_displacement: f64,
    rebuilds: u64,
}

impl CellList {
    /// Build a cell list for the given positions and box.
    pub fn new(positions: &[[f64; 3]], box_lengths: &[f64; 3], cutoff: f64) -> Self {
        let n = positions.len();
        // Determine number of cells per dimension (at least 3 to avoid self-interaction issues)
        let nc = [
            (box_lengths[0] / cutoff).floor().max(3.0) as usize,
            (box_lengths[1] / cutoff).floor().max(3.0) as usize,
            (box_lengths[2] / cutoff).floor().max(3.0) as usize,
        ];
        let cell_size = [
            box_lengths[0] / nc[0] as f64,
            box_lengths[1] / nc[1] as f64,
            box_lengths[2] / nc[2] as f64,
        ];
        let n_cells = nc[0] * nc[1] * nc[2];

        let mut head = vec![usize::MAX; n_cells];
        let mut next = vec![usize::MAX; n];
        let mut cell_of = vec![0usize; n];

        for i in 0..n {
            let ci = cell_index(&positions[i], &cell_size, &nc, box_lengths);
            cell_of[i] = ci;
            next[i] = head[ci];
            head[ci] = i;
        }

        CellList {
            nc,
            cell_size,
            box_lengths: *box_lengths,
            head,
            next,
            cell_of,
            n_cells,
        }
    }

    /// Move atom `idx` from old_pos to new_pos, updating the cell list.
    pub fn move_atom(&mut self, idx: usize, new_pos: &[f64; 3]) {
        let new_ci = cell_index(new_pos, &self.cell_size, &self.nc, &self.box_lengths);
        let old_ci = self.cell_of[idx];
        if new_ci == old_ci {
            return;
        }

        // Remove from old cell
        if self.head[old_ci] == idx {
            self.head[old_ci] = self.next[idx];
        } else {
            let mut prev = self.head[old_ci];
            while prev != usize::MAX {
                if self.next[prev] == idx {
                    self.next[prev] = self.next[idx];
                    break;
                }
                prev = self.next[prev];
            }
        }

        // Insert into new cell
        self.next[idx] = self.head[new_ci];
        self.head[new_ci] = idx;
        self.cell_of[idx] = new_ci;
    }

    /// Get all 27 neighbor cell indices for a given cell (no allocation).
    #[inline]
    pub fn neighbor_cells(&self, cell_idx: usize) -> [usize; 27] {
        let cz = cell_idx / (self.nc[0] * self.nc[1]);
        let rem = cell_idx % (self.nc[0] * self.nc[1]);
        let cy = rem / self.nc[0];
        let cx = rem % self.nc[0];

        let mut cells = [0usize; 27];
        let mut idx = 0;
        for dz in [-1i32, 0, 1] {
            for dy in [-1i32, 0, 1] {
                for dx in [-1i32, 0, 1] {
                    let nx = ((cx as i32 + dx + self.nc[0] as i32) % self.nc[0] as i32) as usize;
                    let ny = ((cy as i32 + dy + self.nc[1] as i32) % self.nc[1] as i32) as usize;
                    let nz = ((cz as i32 + dz + self.nc[2] as i32) % self.nc[2] as i32) as usize;
                    cells[idx] = nz * self.nc[0] * self.nc[1] + ny * self.nc[0] + nx;
                    idx += 1;
                }
            }
        }
        cells
    }

    /// Iterate over all atoms in a cell.
    #[inline]
    pub fn atoms_in_cell(&self, cell_idx: usize) -> CellIter<'_> {
        CellIter {
            cell_list: self,
            current: self.head[cell_idx],
        }
    }
}

impl VerletNeighborList {
    /// Build a directed neighbor list. Every atom owns a contiguous slice of
    /// candidate indices, and every physical pair therefore appears twice.
    pub fn new(positions: &[[f64; 3]], box_lengths: &[f64; 3], cutoff: f64, skin: f64) -> Self {
        assert!(cutoff.is_finite() && cutoff > 0.0);
        assert!(skin.is_finite() && skin > 0.0);
        let mut list = Self {
            cutoff,
            skin,
            list_radius: cutoff + skin,
            box_lengths: *box_lengths,
            reference_positions: Vec::new(),
            offsets: Vec::new(),
            neighbors: Vec::new(),
            max_displacement: 0.0,
            rebuilds: 0,
        };
        list.rebuild(positions);
        list
    }

    /// Candidate atoms for `atom_idx`, including the Verlet skin.
    #[inline]
    pub fn neighbors(&self, atom_idx: usize) -> &[usize] {
        &self.neighbors[self.offsets[atom_idx]..self.offsets[atom_idx + 1]]
    }

    /// Prepare the list for a proposed move. A rebuild is performed before
    /// evaluating the trial whenever the current skin can no longer guarantee
    /// that all cutoff pairs are present.
    pub fn prepare_trial(&mut self, positions: &[[f64; 3]], atom_idx: usize, new_pos: &[f64; 3]) {
        let trial_displacement = periodic_distance(
            new_pos,
            &self.reference_positions[atom_idx],
            &self.box_lengths,
        );
        if trial_displacement + self.max_displacement <= self.skin {
            return;
        }

        // If one proposed step itself exceeds the current skin, enlarge the
        // skin before rebuilding. This is unusual for RMC but keeps the
        // validity rule exact for arbitrary configured step sizes.
        let step = periodic_distance(new_pos, &positions[atom_idx], &self.box_lengths);
        if step >= self.skin {
            self.skin = step * 1.05 + f64::EPSILON;
            self.list_radius = self.cutoff + self.skin;
        }
        self.rebuild(positions);
    }

    /// Record an accepted coordinate for conservative rebuild tracking.
    #[inline]
    pub fn accept_move(&mut self, atom_idx: usize, new_pos: &[f64; 3]) {
        let displacement = periodic_distance(
            new_pos,
            &self.reference_positions[atom_idx],
            &self.box_lengths,
        );
        self.max_displacement = self.max_displacement.max(displacement);
    }

    pub fn cutoff(&self) -> f64 {
        self.cutoff
    }

    pub fn skin(&self) -> f64 {
        self.skin
    }

    pub fn list_radius(&self) -> f64 {
        self.list_radius
    }

    pub fn rebuilds(&self) -> u64 {
        self.rebuilds
    }

    pub fn total_neighbors(&self) -> usize {
        self.neighbors.len()
    }

    pub fn max_neighbors(&self) -> usize {
        self.offsets
            .windows(2)
            .map(|window| window[1] - window[0])
            .max()
            .unwrap_or(0)
    }

    fn rebuild(&mut self, positions: &[[f64; 3]]) {
        self.reference_positions.clear();
        self.reference_positions.extend_from_slice(positions);
        self.offsets.clear();
        self.offsets.reserve(positions.len() + 1);
        self.neighbors.clear();
        self.offsets.push(0);

        // Three cells per list radius is a useful compromise: it keeps bucket
        // occupancy low while making construction much cheaper than O(N^2).
        let builder_width = self.list_radius / 3.0;
        let cells = CellList::new(positions, &self.box_lengths, builder_width);
        let cell_radius = [
            (self.list_radius / cells.cell_size[0]).ceil() as i32,
            (self.list_radius / cells.cell_size[1]).ceil() as i32,
            (self.list_radius / cells.cell_size[2]).ceil() as i32,
        ];
        let list_radius2 = self.list_radius * self.list_radius;
        let mut seen_cells = vec![0u32; cells.n_cells];
        let mut stamp = 0u32;

        for (atom_idx, pos) in positions.iter().enumerate() {
            stamp = stamp.wrapping_add(1);
            if stamp == 0 {
                seen_cells.fill(0);
                stamp = 1;
            }
            let center = cell_coordinates(cells.cell_of[atom_idx], &cells.nc);

            for dz in -cell_radius[2]..=cell_radius[2] {
                let min_z = minimum_cell_separation(dz, cells.cell_size[2]);
                for dy in -cell_radius[1]..=cell_radius[1] {
                    let min_y = minimum_cell_separation(dy, cells.cell_size[1]);
                    for dx in -cell_radius[0]..=cell_radius[0] {
                        let min_x = minimum_cell_separation(dx, cells.cell_size[0]);
                        if min_x * min_x + min_y * min_y + min_z * min_z > list_radius2 {
                            continue;
                        }
                        let coordinates = [
                            periodic_cell(center[0], dx, cells.nc[0]),
                            periodic_cell(center[1], dy, cells.nc[1]),
                            periodic_cell(center[2], dz, cells.nc[2]),
                        ];
                        let cell = linear_cell(coordinates, &cells.nc);
                        if seen_cells[cell] == stamp {
                            continue;
                        }
                        seen_cells[cell] = stamp;

                        for candidate in cells.atoms_in_cell(cell) {
                            if candidate == atom_idx {
                                continue;
                            }
                            if periodic_distance_squared(
                                pos,
                                &positions[candidate],
                                &self.box_lengths,
                            ) < list_radius2
                            {
                                self.neighbors.push(candidate);
                            }
                        }
                    }
                }
            }
            self.offsets.push(self.neighbors.len());
        }
        self.max_displacement = 0.0;
        self.rebuilds += 1;
    }
}

pub struct CellIter<'a> {
    cell_list: &'a CellList,
    current: usize,
}

impl<'a> Iterator for CellIter<'a> {
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<usize> {
        if self.current == usize::MAX {
            return None;
        }
        let idx = self.current;
        self.current = self.cell_list.next[idx];
        Some(idx)
    }
}

#[inline]
fn cell_index(
    pos: &[f64; 3],
    _cell_size: &[f64; 3],
    nc: &[usize; 3],
    box_lengths: &[f64; 3],
) -> usize {
    let cx = ((pos[0] / box_lengths[0]).fract() * nc[0] as f64).floor() as usize;
    let cy = ((pos[1] / box_lengths[1]).fract() * nc[1] as f64).floor() as usize;
    let cz = ((pos[2] / box_lengths[2]).fract() * nc[2] as f64).floor() as usize;
    // Clamp to valid range (handles edge cases from floating point)
    let cx = cx.min(nc[0] - 1);
    let cy = cy.min(nc[1] - 1);
    let cz = cz.min(nc[2] - 1);
    cz * nc[0] * nc[1] + cy * nc[0] + cx
}

#[inline]
fn cell_coordinates(cell_idx: usize, nc: &[usize; 3]) -> [usize; 3] {
    let xy = nc[0] * nc[1];
    let z = cell_idx / xy;
    let remainder = cell_idx % xy;
    [remainder % nc[0], remainder / nc[0], z]
}

#[inline]
fn linear_cell(coordinates: [usize; 3], nc: &[usize; 3]) -> usize {
    coordinates[2] * nc[0] * nc[1] + coordinates[1] * nc[0] + coordinates[0]
}

#[inline]
fn periodic_cell(cell: usize, offset: i32, count: usize) -> usize {
    (cell as i32 + offset).rem_euclid(count as i32) as usize
}

#[inline]
fn minimum_cell_separation(offset: i32, cell_size: f64) -> f64 {
    (offset.abs() - 1).max(0) as f64 * cell_size
}

#[inline]
fn periodic_distance(a: &[f64; 3], b: &[f64; 3], box_lengths: &[f64; 3]) -> f64 {
    periodic_distance_squared(a, b, box_lengths).sqrt()
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

#[cfg(test)]
mod tests {
    use super::VerletNeighborList;

    fn distance2(a: &[f64; 3], b: &[f64; 3], box_lengths: &[f64; 3]) -> f64 {
        super::periodic_distance_squared(a, b, box_lengths)
    }

    #[test]
    fn verlet_list_contains_every_pair_inside_list_radius_without_duplicates() {
        let box_lengths = [17.0, 19.0, 23.0];
        let positions: Vec<[f64; 3]> = (0..91)
            .map(|i| {
                [
                    (i * 37 % 169) as f64 / 169.0 * box_lengths[0],
                    (i * 53 % 181) as f64 / 181.0 * box_lengths[1],
                    (i * 71 % 211) as f64 / 211.0 * box_lengths[2],
                ]
            })
            .collect();
        let list = VerletNeighborList::new(&positions, &box_lengths, 4.0, 1.2);
        let radius2 = list.list_radius().powi(2);

        for i in 0..positions.len() {
            let mut observed = list.neighbors(i).to_vec();
            observed.sort_unstable();
            let original_len = observed.len();
            observed.dedup();
            assert_eq!(
                observed.len(),
                original_len,
                "duplicate neighbor for atom {i}"
            );
            assert!(!observed.contains(&i));

            let expected: Vec<usize> = (0..positions.len())
                .filter(|&j| {
                    j != i && distance2(&positions[i], &positions[j], &box_lengths) < radius2
                })
                .collect();
            assert_eq!(observed, expected, "neighbor mismatch for atom {i}");
        }
    }

    #[test]
    fn verlet_list_rebuilds_before_a_pair_can_enter_the_cutoff() {
        let box_lengths = [20.0; 3];
        let mut positions = vec![[1.0, 1.0, 1.0], [7.1, 1.0, 1.0], [15.0, 15.0, 15.0]];
        let mut list = VerletNeighborList::new(&positions, &box_lengths, 5.0, 1.0);
        assert!(!list.neighbors(0).contains(&1));

        let first_move = [1.6, 1.0, 1.0];
        list.prepare_trial(&positions, 0, &first_move);
        positions[0] = first_move;
        list.accept_move(0, &first_move);

        let second_move = [6.6, 1.0, 1.0];
        list.prepare_trial(&positions, 1, &second_move);
        assert!(list.rebuilds() >= 2);
        assert!(list.neighbors(1).contains(&0));
    }
}
