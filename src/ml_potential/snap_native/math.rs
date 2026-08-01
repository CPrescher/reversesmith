//! Quantum-number indexing and coupling coefficients used by SNAP.
//!
//! SNAP stores doubled angular momenta as integers, so half-integer quantum
//! numbers never need floating-point comparisons. The coefficient convention
//! here is the Condon--Shortley convention used in the published SNAP
//! bispectrum definition.

use rustfft::num_complex::Complex64;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct BispectrumIndex {
    pub two_j1: usize,
    pub two_j2: usize,
    pub two_j: usize,
}

/// Return independent bispectrum components in the order used by LAMMPS SNAP
/// coefficient files.
pub(super) fn bispectrum_indices(two_j_max: usize) -> Vec<BispectrumIndex> {
    let mut indices = Vec::new();

    for two_j1 in 0..=two_j_max {
        for two_j2 in 0..=two_j1 {
            let upper = two_j_max.min(two_j1 + two_j2);
            for two_j in ((two_j1 - two_j2)..=upper).step_by(2) {
                if two_j >= two_j1 {
                    indices.push(BispectrumIndex {
                        two_j1,
                        two_j2,
                        two_j,
                    });
                }
            }
        }
    }

    indices
}

#[derive(Clone, Debug)]
pub(super) struct ClebschGordanBlock {
    two_j1: usize,
    two_j2: usize,
    two_j: usize,
    coefficients: Vec<f64>,
}

impl ClebschGordanBlock {
    pub fn new(index: BispectrumIndex) -> Self {
        let mut coefficients = vec![0.0; (index.two_j1 + 1) * (index.two_j2 + 1)];

        for m1_index in 0..=index.two_j1 {
            let two_m1 = two_m(m1_index, index.two_j1);
            for m2_index in 0..=index.two_j2 {
                let two_m2 = two_m(m2_index, index.two_j2);
                let two_m = two_m1 + two_m2;
                coefficients[m1_index * (index.two_j2 + 1) + m2_index] = clebsch_gordan(
                    index.two_j1 as i32,
                    two_m1,
                    index.two_j2 as i32,
                    two_m2,
                    index.two_j as i32,
                    two_m,
                );
            }
        }

        Self {
            two_j1: index.two_j1,
            two_j2: index.two_j2,
            two_j: index.two_j,
            coefficients,
        }
    }

    pub fn coefficient(&self, m1_index: usize, m2_index: usize) -> f64 {
        debug_assert!(m1_index <= self.two_j1);
        debug_assert!(m2_index <= self.two_j2);
        self.coefficients[m1_index * (self.two_j2 + 1) + m2_index]
    }

    pub fn coupled_m_index(&self, m1_index: usize, m2_index: usize) -> Option<usize> {
        let two_m = two_m(m1_index, self.two_j1) + two_m(m2_index, self.two_j2);
        magnetic_index(two_m, self.two_j)
    }
}

pub(super) fn coupling_blocks(two_j_max: usize) -> Vec<ClebschGordanBlock> {
    bispectrum_indices(two_j_max)
        .into_iter()
        .map(ClebschGordanBlock::new)
        .collect()
}

#[derive(Clone, Debug)]
pub(super) struct SnapBasis {
    two_j_max: usize,
    coupling_blocks: Vec<ClebschGordanBlock>,
}

impl SnapBasis {
    pub fn new(two_j_max: usize) -> Result<Self, String> {
        let coupling_blocks = coupling_blocks(two_j_max);
        for block in &coupling_blocks {
            if !block.is_normalized(1.0e-12) {
                return Err(format!(
                    "failed to construct normalized SNAP coupling table ({}, {}, {})",
                    block.two_j1, block.two_j2, block.two_j
                ));
            }
        }
        Ok(Self {
            two_j_max,
            coupling_blocks,
        })
    }

    pub fn descriptor_count(&self) -> usize {
        self.coupling_blocks.len()
    }

    pub fn empty_density(&self) -> Vec<AngularMatrix> {
        (0..=self.two_j_max).map(AngularMatrix::identity).collect()
    }

    pub fn zero_density(&self) -> Vec<AngularMatrix> {
        (0..=self.two_j_max).map(AngularMatrix::zeros).collect()
    }

    pub fn add_neighbor(
        &self,
        density: &mut [AngularMatrix],
        displacement: [f64; 3],
        theta: f64,
        weight: f64,
    ) {
        debug_assert_eq!(density.len(), self.two_j_max + 1);
        let matrices = wigner_u_matrices(self.two_j_max, displacement, theta);
        for (total, neighbor) in density.iter_mut().zip(matrices) {
            total.add_scaled(&neighbor, weight);
        }
    }

    pub fn bispectrum(
        &self,
        density: &[AngularMatrix],
        normalize: bool,
        subtract_isolated_atom: bool,
    ) -> Vec<f64> {
        self.mixed_bispectrum(density, density, density, normalize, subtract_isolated_atom)
    }

    pub fn mixed_bispectrum(
        &self,
        density1: &[AngularMatrix],
        density2: &[AngularMatrix],
        density3: &[AngularMatrix],
        normalize: bool,
        subtract_isolated_atom: bool,
    ) -> Vec<f64> {
        debug_assert_eq!(density1.len(), self.two_j_max + 1);
        debug_assert_eq!(density2.len(), self.two_j_max + 1);
        debug_assert_eq!(density3.len(), self.two_j_max + 1);
        self.coupling_blocks
            .iter()
            .map(|block| {
                let u1 = &density1[block.two_j1];
                let u2 = &density2[block.two_j2];
                let u = &density3[block.two_j];
                let mut coupled = AngularMatrix::zeros(block.two_j);

                for m1_index in 0..=block.two_j1 {
                    for m2_index in 0..=block.two_j2 {
                        let Some(m_index) = block.coupled_m_index(m1_index, m2_index) else {
                            continue;
                        };
                        let m_coupling = block.coefficient(m1_index, m2_index);
                        for mp1_index in 0..=block.two_j1 {
                            for mp2_index in 0..=block.two_j2 {
                                let Some(mp_index) = block.coupled_m_index(mp1_index, mp2_index)
                                else {
                                    continue;
                                };
                                let coupling = m_coupling * block.coefficient(mp1_index, mp2_index);
                                coupled[(m_index, mp_index)] += coupling
                                    * u1[(m1_index, mp1_index)]
                                    * u2[(m2_index, mp2_index)];
                            }
                        }
                    }
                }

                let normalization = if normalize {
                    (block.two_j + 1) as f64
                } else {
                    1.0
                };
                let mut value = u
                    .values
                    .iter()
                    .zip(&coupled.values)
                    .map(|(u_value, coupled_value)| (u_value.conj() * coupled_value).re)
                    .sum::<f64>()
                    / normalization;
                if subtract_isolated_atom {
                    value -= if normalize {
                        1.0
                    } else {
                        (block.two_j + 1) as f64
                    };
                }
                value
            })
            .collect()
    }
}

#[derive(Clone, Debug)]
pub(super) struct AngularMatrix {
    two_j: usize,
    values: Vec<Complex64>,
}

impl AngularMatrix {
    fn zeros(two_j: usize) -> Self {
        Self {
            two_j,
            values: vec![Complex64::new(0.0, 0.0); (two_j + 1).pow(2)],
        }
    }

    fn identity(two_j: usize) -> Self {
        let mut matrix = Self::zeros(two_j);
        for index in 0..=two_j {
            matrix[(index, index)] = Complex64::new(1.0, 0.0);
        }
        matrix
    }

    fn add_scaled(&mut self, other: &Self, scale: f64) {
        debug_assert_eq!(self.two_j, other.two_j);
        for (target, value) in self.values.iter_mut().zip(&other.values) {
            *target += scale * value;
        }
    }
}

impl std::ops::Index<(usize, usize)> for AngularMatrix {
    type Output = Complex64;

    fn index(&self, (m_index, mp_index): (usize, usize)) -> &Self::Output {
        &self.values[mp_index * (self.two_j + 1) + m_index]
    }
}

impl std::ops::IndexMut<(usize, usize)> for AngularMatrix {
    fn index_mut(&mut self, (m_index, mp_index): (usize, usize)) -> &mut Self::Output {
        &mut self.values[mp_index * (self.two_j + 1) + m_index]
    }
}

fn wigner_u_matrices(two_j_max: usize, displacement: [f64; 3], theta: f64) -> Vec<AngularMatrix> {
    let (a, b) = cayley_klein_parameters(displacement, theta);

    let mut matrices = Vec::with_capacity(two_j_max + 1);
    matrices.push(AngularMatrix::identity(0));

    for two_j in 1..=two_j_max {
        let previous = &matrices[two_j - 1];
        let mut current = AngularMatrix::zeros(two_j);

        // Recouple the previous representation with the spin-1/2
        // representation. Only the left half is independent.
        for mp_index in 0..=two_j / 2 {
            for m_index in 0..two_j {
                let previous_value = previous[(m_index, mp_index)];
                let denominator = (two_j - mp_index) as f64;
                let a_factor = ((two_j - m_index) as f64 / denominator).sqrt();
                let b_factor = ((m_index + 1) as f64 / denominator).sqrt();
                current[(m_index, mp_index)] += a_factor * a.conj() * previous_value;
                current[(m_index + 1, mp_index)] -= b_factor * b.conj() * previous_value;
            }
        }

        // Recover the right half from inversion symmetry. For even two_j,
        // the middle column was already calculated by the recurrence.
        for mp_index in 0..two_j.div_ceil(2) {
            for m_index in 0..=two_j {
                let sign = if (m_index + mp_index) % 2 == 0 {
                    1.0
                } else {
                    -1.0
                };
                current[(two_j - m_index, two_j - mp_index)] =
                    sign * current[(m_index, mp_index)].conj();
            }
        }

        matrices.push(current);
    }

    matrices
}

fn cayley_klein_parameters(displacement: [f64; 3], theta: f64) -> (Complex64, Complex64) {
    let [x, y, z] = displacement;
    let radius = (x * x + y * y + z * z).sqrt();
    debug_assert!(radius > 0.0);

    // Cayley--Klein parameters for the unit quaternion obtained by mapping
    // the three-dimensional neighbor position onto the three-sphere. This is
    // algebraically identical to z0=r/tan(theta), followed by normalization
    // of (z0, -z, y, -x). The sign is therefore required when theta is
    // negative; using cos(theta) alone would not match that mapping.
    let sine = theta.sin();
    let direction_scale = sine.abs() / radius;
    let sine_sign = if sine < 0.0 { -1.0 } else { 1.0 };
    let a = Complex64::new(sine_sign * theta.cos(), -direction_scale * z);
    let b = Complex64::new(direction_scale * y, -direction_scale * x);
    (a, b)
}

impl ClebschGordanBlock {
    fn is_normalized(&self, tolerance: f64) -> bool {
        (0..=self.two_j).all(|m_index| {
            let norm = (0..=self.two_j1)
                .flat_map(|m1_index| (0..=self.two_j2).map(move |m2_index| (m1_index, m2_index)))
                .filter(|(m1_index, m2_index)| {
                    self.coupled_m_index(*m1_index, *m2_index) == Some(m_index)
                })
                .map(|(m1_index, m2_index)| self.coefficient(m1_index, m2_index).powi(2))
                .sum::<f64>();
            (norm - 1.0).abs() <= tolerance
        })
    }
}

fn two_m(index: usize, two_j: usize) -> i32 {
    2 * index as i32 - two_j as i32
}

fn magnetic_index(two_m: i32, two_j: usize) -> Option<usize> {
    let shifted = two_m + two_j as i32;
    if shifted < 0 || shifted % 2 != 0 {
        return None;
    }
    let index = (shifted / 2) as usize;
    (index <= two_j).then_some(index)
}

fn valid_magnetic_quantum_number(two_j: i32, two_m: i32) -> bool {
    two_j >= 0 && two_m.abs() <= two_j && (two_j - two_m) % 2 == 0
}

fn half_integer(value: i32) -> Option<usize> {
    (value >= 0 && value % 2 == 0).then_some((value / 2) as usize)
}

fn factorial(value: usize) -> f64 {
    (1..=value).fold(1.0, |product, factor| product * factor as f64)
}

/// Clebsch--Gordan coefficient using doubled integer quantum numbers.
fn clebsch_gordan(
    two_j1: i32,
    two_m1: i32,
    two_j2: i32,
    two_m2: i32,
    two_j: i32,
    two_m: i32,
) -> f64 {
    if !valid_magnetic_quantum_number(two_j1, two_m1)
        || !valid_magnetic_quantum_number(two_j2, two_m2)
        || !valid_magnetic_quantum_number(two_j, two_m)
        || two_m1 + two_m2 != two_m
        || two_j < (two_j1 - two_j2).abs()
        || two_j > two_j1 + two_j2
        || (two_j1 + two_j2 + two_j) % 2 != 0
    {
        return 0.0;
    }

    let Some(triangle_a) = half_integer(two_j + two_j1 - two_j2) else {
        return 0.0;
    };
    let Some(triangle_b) = half_integer(two_j - two_j1 + two_j2) else {
        return 0.0;
    };
    let Some(triangle_c) = half_integer(two_j1 + two_j2 - two_j) else {
        return 0.0;
    };
    let triangle_denominator = ((two_j1 + two_j2 + two_j) / 2 + 1) as usize;

    let triangle_factor = ((two_j + 1) as f64
        * factorial(triangle_a)
        * factorial(triangle_b)
        * factorial(triangle_c)
        / factorial(triangle_denominator))
    .sqrt();

    let magnetic_arguments = [
        (two_j + two_m) / 2,
        (two_j - two_m) / 2,
        (two_j1 - two_m1) / 2,
        (two_j1 + two_m1) / 2,
        (two_j2 - two_m2) / 2,
        (two_j2 + two_m2) / 2,
    ];
    let magnetic_factor = magnetic_arguments
        .into_iter()
        .map(|argument| factorial(argument as usize))
        .product::<f64>()
        .sqrt();

    let sum_limit = triangle_c
        .max(((two_j1 - two_m1) / 2) as usize)
        .max(((two_j2 + two_m2) / 2) as usize);
    let mut sum = 0.0;
    for k in 0..=sum_limit {
        let denominator_arguments = [
            k as i32,
            triangle_c as i32 - k as i32,
            (two_j1 - two_m1) / 2 - k as i32,
            (two_j2 + two_m2) / 2 - k as i32,
            (two_j - two_j2 + two_m1) / 2 + k as i32,
            (two_j - two_j1 - two_m2) / 2 + k as i32,
        ];
        if denominator_arguments.iter().any(|argument| *argument < 0) {
            continue;
        }
        let denominator = denominator_arguments
            .into_iter()
            .map(|argument| factorial(argument as usize))
            .product::<f64>();
        let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
        sum += sign / denominator;
    }

    triangle_factor * magnetic_factor * sum
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 1.0e-13,
            "expected {expected:.16e}, got {actual:.16e}"
        );
    }

    #[test]
    fn bispectrum_index_order_matches_snap_coefficient_files() {
        assert_eq!(
            bispectrum_indices(2),
            [
                BispectrumIndex {
                    two_j1: 0,
                    two_j2: 0,
                    two_j: 0
                },
                BispectrumIndex {
                    two_j1: 1,
                    two_j2: 0,
                    two_j: 1
                },
                BispectrumIndex {
                    two_j1: 1,
                    two_j2: 1,
                    two_j: 2
                },
                BispectrumIndex {
                    two_j1: 2,
                    two_j2: 0,
                    two_j: 2
                },
                BispectrumIndex {
                    two_j1: 2,
                    two_j2: 2,
                    two_j: 2
                }
            ]
        );
    }

    #[test]
    fn descriptor_counts_match_the_lammps_pyramidal_sequence() {
        let expected = [1, 2, 5, 8, 14, 20, 30, 40, 55];
        for (two_j_max, expected_count) in expected.into_iter().enumerate() {
            assert_eq!(bispectrum_indices(two_j_max).len(), expected_count);
        }
    }

    #[test]
    fn clebsch_gordan_known_half_integer_coefficients() {
        let inverse_sqrt_two = 1.0 / 2.0_f64.sqrt();
        assert_close(clebsch_gordan(1, 1, 1, -1, 0, 0), inverse_sqrt_two);
        assert_close(clebsch_gordan(1, -1, 1, 1, 0, 0), -inverse_sqrt_two);
        assert_close(clebsch_gordan(1, 1, 1, -1, 2, 0), inverse_sqrt_two);
        assert_close(clebsch_gordan(1, -1, 1, 1, 2, 0), inverse_sqrt_two);
    }

    #[test]
    fn clebsch_gordan_blocks_are_normalized_for_each_coupled_m() {
        for block in coupling_blocks(8) {
            for m_index in 0..=block.two_j {
                let mut norm = 0.0;
                for m1_index in 0..=block.two_j1 {
                    for m2_index in 0..=block.two_j2 {
                        if block.coupled_m_index(m1_index, m2_index) == Some(m_index) {
                            norm += block.coefficient(m1_index, m2_index).powi(2);
                        }
                    }
                }
                assert_close(norm, 1.0);
            }
        }
    }

    #[test]
    fn wigner_matrices_are_unitary() {
        let matrices = wigner_u_matrices(8, [1.2, -0.7, 2.1], 1.137);
        for matrix in matrices {
            for row in 0..=matrix.two_j {
                for other_row in 0..=matrix.two_j {
                    let inner_product = (0..=matrix.two_j)
                        .map(|column| matrix[(row, column)] * matrix[(other_row, column)].conj())
                        .sum::<Complex64>();
                    let expected = if row == other_row { 1.0 } else { 0.0 };
                    assert_close(inner_product.re, expected);
                    assert_close(inner_product.im, 0.0);
                }
            }
        }
    }

    #[test]
    fn isolated_atom_bispectrum_is_removed_by_bzero() {
        let basis = SnapBasis::new(8).unwrap();
        let density = basis.empty_density();
        let unnormalized = basis.bispectrum(&density, false, false);
        let normalized = basis.bispectrum(&density, true, false);
        for (block, (unnormalized, normalized)) in basis
            .coupling_blocks
            .iter()
            .zip(unnormalized.into_iter().zip(normalized))
        {
            assert_close(unnormalized, (block.two_j + 1) as f64);
            assert_close(normalized, 1.0);
        }
        assert!(basis
            .bispectrum(&density, false, true)
            .iter()
            .all(|value| value.abs() < 1.0e-13));
        assert!(basis
            .bispectrum(&density, true, true)
            .iter()
            .all(|value| value.abs() < 1.0e-13));
    }

    #[test]
    fn signed_cayley_mapping_matches_the_z0_definition() {
        let displacement = [1.2, -0.7, 2.1];
        let theta = -0.731_f64;
        let radius = displacement
            .iter()
            .map(|value| value * value)
            .sum::<f64>()
            .sqrt();
        let z0 = radius / theta.tan();
        let inverse_norm = 1.0 / (radius * radius + z0 * z0).sqrt();
        let (a, b) = cayley_klein_parameters(displacement, theta);

        assert_close(a.re, z0 * inverse_norm);
        assert_close(a.im, -displacement[2] * inverse_norm);
        assert_close(b.re, displacement[1] * inverse_norm);
        assert_close(b.im, -displacement[0] * inverse_norm);
    }

    #[test]
    fn invalid_couplings_are_zero() {
        assert_eq!(clebsch_gordan(1, 0, 1, 1, 1, 1), 0.0);
        assert_eq!(clebsch_gordan(1, 1, 1, 1, 0, 2), 0.0);
        assert_eq!(clebsch_gordan(1, 1, 1, 1, 4, 2), 0.0);
        assert_eq!(clebsch_gordan(1, 1, 1, -1, 2, 2), 0.0);
    }
}
