/// Cell list for O(1)-per-neighbor spatial lookups with periodic boundary conditions.

pub struct CellList {
    /// Number of cells per dimension.
    pub nc: [usize; 3],
    /// Cell size per dimension (>= cutoff).
    pub cell_size: [f64; 3],
    /// Box lengths.
    pub box_lengths: [f64; 3],
    /// head[cell_idx] = first atom in cell, or usize::MAX if empty.
    pub head: Vec<usize>,
    /// next[atom_idx] = next atom in same cell, or usize::MAX if end.
    pub next: Vec<usize>,
    /// cell_of[atom_idx] = cell index for this atom.
    pub cell_of: Vec<usize>,
    /// Precomputed neighbor cell offsets (including self).
    pub neighbor_offsets: Vec<usize>,
    /// Number of cells total.
    pub n_cells: usize,
}

impl CellList {
    /// Build a cell list for the given positions and box.
    pub fn new(
        positions: &[[f64; 3]],
        box_lengths: &[f64; 3],
        cutoff: f64,
    ) -> Self {
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

        // Precompute neighbor cell offsets (27 neighbors in 3D)
        let mut neighbor_offsets = Vec::with_capacity(27);
        for dz in [-1i32, 0, 1] {
            for dy in [-1i32, 0, 1] {
                for dx in [-1i32, 0, 1] {
                    // Compute offset using modular arithmetic trick
                    let ox = ((dx + nc[0] as i32) % nc[0] as i32) as usize;
                    let oy = ((dy + nc[1] as i32) % nc[1] as i32) as usize;
                    let oz = ((dz + nc[2] as i32) % nc[2] as i32) as usize;
                    neighbor_offsets.push(oz * nc[0] * nc[1] + oy * nc[0] + ox);
                }
            }
        }

        CellList {
            nc,
            cell_size,
            box_lengths: *box_lengths,
            head,
            next,
            cell_of,
            neighbor_offsets,
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
        CellIter { cell_list: self, current: self.head[cell_idx] }
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
fn cell_index(pos: &[f64; 3], _cell_size: &[f64; 3], nc: &[usize; 3], box_lengths: &[f64; 3]) -> usize {
    let cx = ((pos[0] / box_lengths[0]).fract() * nc[0] as f64).floor() as usize;
    let cy = ((pos[1] / box_lengths[1]).fract() * nc[1] as f64).floor() as usize;
    let cz = ((pos[2] / box_lengths[2]).fract() * nc[2] as f64).floor() as usize;
    // Clamp to valid range (handles edge cases from floating point)
    let cx = cx.min(nc[0] - 1);
    let cy = cy.min(nc[1] - 1);
    let cz = cz.min(nc[2] - 1);
    cz * nc[0] * nc[1] + cy * nc[0] + cx
}
