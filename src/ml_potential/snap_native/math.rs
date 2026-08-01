//! Quantum-number indexing and coupling coefficients used by SNAP.
//!
//! SNAP stores doubled angular momenta as integers, so half-integer quantum
//! numbers never need floating-point comparisons. The coefficient convention
//! here is the Condon--Shortley convention used in the published SNAP
//! bispectrum definition.

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
        Ok(Self { coupling_blocks })
    }

    pub fn descriptor_count(&self) -> usize {
        self.coupling_blocks.len()
    }
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
    fn invalid_couplings_are_zero() {
        assert_eq!(clebsch_gordan(1, 0, 1, 1, 1, 1), 0.0);
        assert_eq!(clebsch_gordan(1, 1, 1, 1, 0, 2), 0.0);
        assert_eq!(clebsch_gordan(1, 1, 1, 1, 4, 2), 0.0);
        assert_eq!(clebsch_gordan(1, 1, 1, -1, 2, 2), 0.0);
    }
}
