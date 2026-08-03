//! Native energy-only PACE/ACE evaluator.
//!
//! This module reads the compact C-tilde YAML (`.yace`) files emitted by
//! pacemaker and consumed by LAMMPS `pair_style pace`.  It implements the
//! product-basis energy expression directly in Rust; LAMMPS is only a test
//! oracle and is not a runtime dependency.

mod runtime;

pub use runtime::PaceNativeModel;

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

const MIN_DISTANCE: f64 = 1.0e-10;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PaceNeighbor {
    /// Minimum-image displacement from the central atom, in Angstrom.
    pub displacement: [f64; 3],
    /// Zero-based rsmith atom-type index.
    pub type_index: usize,
}

#[derive(Debug, Clone)]
pub struct PaceModel {
    pub elements: Vec<String>,
    pub source_path: PathBuf,
    e0: Vec<f64>,
    embeddings: Vec<Embedding>,
    bonds: Vec<Vec<Bond>>,
    functions: Vec<Vec<BasisFunction>>,
    type_to_element: Vec<usize>,
    maximum_cutoff: f64,
    spline_step: f64,
}

#[derive(Debug, Clone)]
struct Embedding {
    ndensity: usize,
    parameters: Vec<f64>,
    kind: EmbeddingKind,
    rho_core_cutoff: f64,
    drho_core_cutoff: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum EmbeddingKind {
    FinnisSinclair,
    FinnisSinclairShiftedScaled,
}

#[derive(Debug, Clone)]
struct Bond {
    nradmax: usize,
    lmax: usize,
    nradbasemax: usize,
    radial_kind: RadialKind,
    lambda: f64,
    coefficients: Vec<f64>,
    prehc: f64,
    lambdahc: f64,
    cutoff: f64,
    cutoff_width: f64,
    inner_cutoff: f64,
    inner_width: f64,
    inner_kind: InnerCutoffKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RadialKind {
    PowerScaled,
    ExponentialCosine,
    LinearScaled,
    SimplifiedBessel,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum InnerCutoffKind {
    Density,
    Distance,
}

#[derive(Debug, Clone)]
struct BasisFunction {
    mu0: usize,
    rank: usize,
    ndensity: usize,
    num_ms_combs: usize,
    mus: Vec<usize>,
    ns: Vec<usize>,
    ls: Vec<usize>,
    ms_combs: Vec<isize>,
    ctildes: Vec<f64>,
}

#[derive(Debug, Clone, Copy, Default)]
struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    const ONE: Self = Self { re: 1.0, im: 0.0 };

    fn conjugate(self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }
}

impl std::ops::AddAssign for Complex {
    fn add_assign(&mut self, rhs: Self) {
        self.re += rhs.re;
        self.im += rhs.im;
    }
}

impl std::ops::Mul for Complex {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self {
            re: self.re * rhs.re - self.im * rhs.im,
            im: self.re * rhs.im + self.im * rhs.re,
        }
    }
}

impl std::ops::Mul<f64> for Complex {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self {
        Self {
            re: self.re * rhs,
            im: self.im * rhs,
        }
    }
}

impl PaceModel {
    pub fn load(path: &Path, type_elements: &[String]) -> Result<Self, String> {
        let text = fs::read_to_string(path)
            .map_err(|error| format!("failed to read PACE model {}: {error}", path.display()))?;
        if text
            .lines()
            .any(|line| line.trim_start().starts_with("nelements="))
        {
            return Err(
                "legacy plain-text PACE .ace files are not supported by pace_native; convert the model to C-tilde YAML (.yace)"
                    .to_string(),
            );
        }
        let mut model = parse_yace(&text, path)?;
        model.type_to_element = resolve_type_elements(&model.elements, type_elements)?;
        Ok(model)
    }

    pub fn maximum_cutoff(&self) -> f64 {
        self.maximum_cutoff
    }

    pub fn element_index_for_type(&self, type_index: usize) -> Result<usize, String> {
        self.type_to_element
            .get(type_index)
            .copied()
            .ok_or_else(|| {
                format!(
                    "PACE atom type index {type_index} is outside the configured species mapping"
                )
            })
    }

    pub fn cutoff_for_types(
        &self,
        central_type: usize,
        neighbor_type: usize,
    ) -> Result<f64, String> {
        let central = self.element_index_for_type(central_type)?;
        let neighbor = self.element_index_for_type(neighbor_type)?;
        Ok(self.bonds[central][neighbor].cutoff)
    }

    pub fn atomic_energy(
        &self,
        central_type: usize,
        neighbors: &[PaceNeighbor],
    ) -> Result<f64, String> {
        let central = self.element_index_for_type(central_type)?;
        let element_count = self.elements.len();
        let max_rank1 = self.functions[central]
            .iter()
            .filter(|function| function.rank == 1)
            .flat_map(|function| function.ns.iter().copied())
            .max()
            .unwrap_or(0);
        let max_n = self.functions[central]
            .iter()
            .filter(|function| function.rank > 1)
            .flat_map(|function| function.ns.iter().copied())
            .max()
            .unwrap_or(0);
        let max_l = self.functions[central]
            .iter()
            .filter(|function| function.rank > 1)
            .flat_map(|function| function.ls.iter().copied())
            .max()
            .unwrap_or(0);

        let mut rank1 = vec![vec![0.0; max_rank1]; element_count];
        let lm_width = 2 * max_l + 1;
        let mut atomic_base =
            vec![vec![vec![Complex::default(); (max_l + 1) * lm_width]; max_n]; element_count];
        let mut rho_core = 0.0;

        for neighbor in neighbors {
            let neighbor_element = self.element_index_for_type(neighbor.type_index)?;
            let bond = &self.bonds[central][neighbor_element];
            if neighbor
                .displacement
                .iter()
                .any(|coordinate| !coordinate.is_finite())
            {
                return Err("PACE neighbor displacement must be finite".to_string());
            }
            let distance_squared = neighbor
                .displacement
                .iter()
                .map(|coordinate| coordinate * coordinate)
                .sum::<f64>();
            if !distance_squared.is_finite() {
                return Err("PACE neighbor distance must be finite".to_string());
            }
            if distance_squared < MIN_DISTANCE * MIN_DISTANCE {
                return Err(format!(
                    "PACE neighbor distance must be at least {MIN_DISTANCE} Angstrom"
                ));
            }
            let distance = distance_squared.sqrt();
            if distance >= bond.cutoff {
                continue;
            }

            let radial = bond.spline_radial(distance, self.spline_step)?;
            for n in 0..max_rank1.min(radial.len()) {
                rank1[neighbor_element][n] += radial[n];
            }

            if max_n > 0 {
                let unit = neighbor.displacement.map(|component| component / distance);
                let harmonics = spherical_harmonics(unit, max_l);
                for n in 0..max_n.min(bond.nradmax) {
                    for l in 0..=max_l.min(bond.lmax) {
                        let radial_nl = (0..bond.nradbasemax)
                            .map(|k| bond.coefficient(n, l, k) * radial[k])
                            .sum::<f64>();
                        for m in -(l as isize)..=(l as isize) {
                            let index = lm_index(l, m, max_l);
                            atomic_base[neighbor_element][n][index] += harmonics[index] * radial_nl;
                        }
                    }
                }
            }
            rho_core += bond.spline_core(distance, self.spline_step)?;
        }

        let embedding = &self.embeddings[central];
        let mut densities = vec![0.0; embedding.ndensity];
        for function in &self.functions[central] {
            if function.rank == 1 {
                let value = rank1[function.mus[0]][function.ns[0] - 1];
                for (density, coefficient) in densities
                    .iter_mut()
                    .zip(function.ctildes.iter().take(function.ndensity))
                {
                    *density += coefficient * value;
                }
                continue;
            }

            for combination in 0..function.num_ms_combs {
                let mut product = Complex::ONE;
                for factor in 0..function.rank {
                    let m = function.ms_combs[combination * function.rank + factor];
                    let index = lm_index(function.ls[factor], m, max_l);
                    product =
                        product * atomic_base[function.mus[factor]][function.ns[factor] - 1][index];
                }
                for (density_index, density) in densities.iter_mut().enumerate() {
                    *density += product.re
                        * function.ctildes[combination * function.ndensity + density_index];
                }
            }
        }

        let ace_energy = embedding.energy(&densities);
        let core_switch = descending_polynomial_switch(
            rho_core,
            embedding.rho_core_cutoff,
            embedding.drho_core_cutoff,
        );
        Ok(ace_energy * core_switch + rho_core + self.e0[central])
    }
}

impl Embedding {
    fn energy(&self, densities: &[f64]) -> f64 {
        densities
            .iter()
            .zip(self.parameters.chunks_exact(2))
            .map(|(density, parameters)| {
                let prefactor = parameters[0];
                let exponent = parameters[1];
                let value = match self.kind {
                    EmbeddingKind::FinnisSinclair => finnis_sinclair(*density, exponent),
                    EmbeddingKind::FinnisSinclairShiftedScaled => {
                        finnis_sinclair_shifted_scaled(*density, exponent)
                    }
                };
                prefactor * value
            })
            .sum()
    }
}

impl Bond {
    fn coefficient(&self, n: usize, l: usize, k: usize) -> f64 {
        self.coefficients[(n * (self.lmax + 1) + l) * self.nradbasemax + k]
    }

    fn spline_radial(&self, distance: f64, step: f64) -> Result<Vec<f64>, String> {
        cubic_lookup(distance, self.cutoff, step, |r| self.raw_radial(r))
    }

    fn spline_core(&self, distance: f64, step: f64) -> Result<f64, String> {
        Ok(cubic_lookup(distance, self.cutoff, step, |r| {
            let (value, derivative) = self.raw_core(r);
            (vec![value], vec![derivative])
        })?[0])
    }

    fn raw_radial(&self, distance: f64) -> (Vec<f64>, Vec<f64>) {
        if distance <= self.inner_cutoff - self.inner_width || distance >= self.cutoff {
            return (vec![0.0; self.nradbasemax], vec![0.0; self.nradbasemax]);
        }

        let (mut values, mut derivatives) = match self.radial_kind {
            RadialKind::PowerScaled => {
                cheb_pow(self.nradbasemax, self.lambda, self.cutoff, distance)
            }
            RadialKind::ExponentialCosine => cheb_exp_cos(
                self.nradbasemax,
                self.lambda,
                self.cutoff,
                self.cutoff_width,
                distance,
            ),
            RadialKind::LinearScaled => cheb_linear(self.nradbasemax, self.cutoff, distance),
            RadialKind::SimplifiedBessel => {
                simplified_bessel(self.nradbasemax, self.cutoff, distance)
            }
        };

        if self.inner_kind == InnerCutoffKind::Distance {
            let (descending, derivative) = descending_polynomial_switch_with_derivative(
                distance,
                self.inner_cutoff,
                self.inner_width,
            );
            let ascending = 1.0 - descending;
            for (value, radial_derivative) in values.iter_mut().zip(derivatives.iter_mut()) {
                *radial_derivative = *radial_derivative * ascending - *value * derivative;
                *value *= ascending;
            }
        }
        (values, derivatives)
    }

    fn raw_core(&self, distance: f64) -> (f64, f64) {
        if self.prehc == 0.0 || distance >= self.cutoff {
            return (0.0, 0.0);
        }
        let prefactor = self.prehc.abs();
        let lambda = self.lambdahc.abs();
        let radius_squared = distance * distance;
        let scaled = lambda * radius_squared;
        if scaled >= 50.0 {
            return (0.0, 0.0);
        }
        let exponential = (-scaled).exp();
        let mut value = prefactor * exponential / distance;
        let mut derivative = -prefactor * exponential * (2.0 * scaled + 1.0) / radius_squared;
        let envelope = 0.5 * (1.0 + (std::f64::consts::PI * distance / self.cutoff).cos());
        let envelope_derivative =
            -0.5 * (std::f64::consts::PI * distance / self.cutoff).sin() * std::f64::consts::PI
                / self.cutoff;
        derivative = value * envelope_derivative + derivative * envelope;
        value *= envelope;

        if self.inner_kind == InnerCutoffKind::Distance {
            let (switch, switch_derivative) = descending_polynomial_switch_with_derivative(
                distance,
                self.inner_cutoff,
                self.inner_width,
            );
            derivative = derivative * switch + value * switch_derivative;
            value *= switch;
        }
        (value, derivative)
    }
}

fn cubic_lookup<F>(distance: f64, cutoff: f64, step: f64, evaluate: F) -> Result<Vec<f64>, String>
where
    F: Fn(f64) -> (Vec<f64>, Vec<f64>),
{
    let bins = (cutoff / step) as usize;
    if bins == 0 {
        return Err("PACE deltaSplineBins is larger than a bond cutoff".to_string());
    }
    let scaled = distance * bins as f64 / cutoff;
    let lower = scaled.floor() as usize;
    if lower == 0 {
        return Err("PACE encountered a distance in the first radial spline bin".to_string());
    }
    if lower >= bins {
        return Ok(evaluate(cutoff).0.into_iter().map(|_| 0.0).collect());
    }
    let spacing = cutoff / bins as f64;
    let (lower_values, lower_derivatives) = evaluate(lower as f64 * spacing);
    let (upper_values, upper_derivatives) = evaluate((lower + 1) as f64 * spacing);
    let fraction = scaled - lower as f64;
    let fraction2 = fraction * fraction;
    let fraction3 = fraction2 * fraction;
    Ok(lower_values
        .iter()
        .zip(&upper_values)
        .zip(lower_derivatives.iter().zip(&upper_derivatives))
        .map(
            |((lower_value, upper_value), (lower_derivative, upper_derivative))| {
                let lower_slope = lower_derivative * spacing;
                let upper_slope = upper_derivative * spacing;
                let c2 = 3.0 * (upper_value - lower_value) - upper_slope - 2.0 * lower_slope;
                let c3 = -2.0 * (upper_value - lower_value) + upper_slope + lower_slope;
                lower_value + lower_slope * fraction + c2 * fraction2 + c3 * fraction3
            },
        )
        .collect())
}

fn chebyshev(count: usize, x: f64) -> (Vec<f64>, Vec<f64>) {
    let mut first = vec![0.0; count + 1];
    let mut second = vec![0.0; count + 1];
    let mut derivative = vec![0.0; count + 1];
    first[0] = 1.0;
    second[0] = 1.0;
    if count > 0 {
        first[1] = x;
        second[1] = 2.0 * x;
    }
    for order in 1..count {
        first[order + 1] = 2.0 * x * first[order] - first[order - 1];
        second[order + 1] = 2.0 * x * second[order] - second[order - 1];
    }
    for order in 1..=count {
        derivative[order] = order as f64 * second[order - 1];
    }
    (first, derivative)
}

fn cheb_pow(count: usize, lambda: f64, cutoff: f64, distance: f64) -> (Vec<f64>, Vec<f64>) {
    let remaining = 1.0 - distance / cutoff;
    let power = remaining.powf(lambda);
    let power_derivative = -lambda / cutoff * remaining.powf(lambda - 1.0);
    let x = 2.0 * (1.0 - power) - 1.0;
    let x_derivative = -2.0 * power_derivative;
    let (polynomials, derivatives) = chebyshev(count, x);
    let mut values = Vec::with_capacity(count);
    let mut radial_derivatives = Vec::with_capacity(count);
    for order in 1..=count {
        values.push(0.5 - 0.5 * polynomials[order]);
        radial_derivatives.push(-0.5 * derivatives[order] * x_derivative);
    }
    (values, radial_derivatives)
}

fn cheb_exp_cos(
    count: usize,
    lambda: f64,
    cutoff: f64,
    cutoff_width: f64,
    distance: f64,
) -> (Vec<f64>, Vec<f64>) {
    let exponential = (-lambda * distance / cutoff).exp();
    let endpoint = (-lambda).exp();
    let x = 1.0 - 2.0 * (exponential - endpoint) / (1.0 - endpoint);
    let x_derivative = 2.0 * lambda / cutoff * exponential / (1.0 - endpoint);
    let (polynomials, polynomial_derivatives) = chebyshev(count.saturating_sub(1), x);
    let mut values = vec![1.0; count];
    let mut derivatives = vec![0.0; count];
    for index in 1..count {
        values[index] = 0.5 - 0.5 * polynomials[index];
        derivatives[index] = -0.5 * polynomial_derivatives[index] * x_derivative;
    }

    let angle = std::f64::consts::PI * distance / cutoff;
    let envelope = 0.5 * (1.0 + angle.cos());
    let envelope_derivative = -0.5 * angle.sin() * std::f64::consts::PI / cutoff;
    for (value, derivative) in values.iter_mut().zip(derivatives.iter_mut()) {
        *derivative = *value * envelope_derivative + *derivative * envelope;
        *value *= envelope;
    }

    let transition_start = cutoff - cutoff_width;
    if cutoff_width > 0.0 && distance > transition_start {
        let angle = std::f64::consts::PI * (distance - transition_start) / cutoff_width;
        let switch = 0.5 * (1.0 + angle.cos());
        let switch_derivative = -0.5 * angle.sin() * std::f64::consts::PI / cutoff_width;
        for (value, derivative) in values.iter_mut().zip(derivatives.iter_mut()) {
            *derivative = *value * switch_derivative + *derivative * switch;
            *value *= switch;
        }
    }
    (values, derivatives)
}

fn cheb_linear(count: usize, cutoff: f64, distance: f64) -> (Vec<f64>, Vec<f64>) {
    let x = 1.0 - distance / cutoff;
    let (polynomials, polynomial_derivatives) = chebyshev(count, x);
    let mut values = Vec::with_capacity(count);
    let mut derivatives = Vec::with_capacity(count);
    for order in 1..=count {
        values.push(0.5 - 0.5 * polynomials[order]);
        derivatives.push(0.5 * polynomial_derivatives[order] / cutoff);
    }
    (values, derivatives)
}

fn sinc(value: f64) -> f64 {
    if value.abs() < 1.0e-8 {
        1.0 - value * value / 6.0
    } else {
        value.sin() / value
    }
}

fn sinc_derivative(value: f64) -> f64 {
    if value.abs() < 1.0e-6 {
        -value / 3.0 + value.powi(3) / 30.0
    } else {
        (value.cos() * value - value.sin()) / (value * value)
    }
}

fn simplified_bessel_aux(distance: f64, cutoff: f64, order: usize) -> (f64, f64) {
    let first = (order + 1) as f64;
    let second = (order + 2) as f64;
    let sign = if order % 2 == 0 { 1.0 } else { -1.0 };
    let prefactor =
        sign * 2.0_f64.sqrt() * std::f64::consts::PI / cutoff.powf(1.5) * first * second
            / (first * first + second * second).sqrt();
    let first_scale = first * std::f64::consts::PI / cutoff;
    let second_scale = second * std::f64::consts::PI / cutoff;
    let value = prefactor * (sinc(distance * first_scale) + sinc(distance * second_scale));
    let derivative = prefactor
        * (sinc_derivative(distance * first_scale) * first_scale
            + sinc_derivative(distance * second_scale) * second_scale);
    (value, derivative)
}

fn simplified_bessel(count: usize, cutoff: f64, distance: f64) -> (Vec<f64>, Vec<f64>) {
    if distance >= cutoff {
        return (vec![0.0; count], vec![0.0; count]);
    }
    let mut values = Vec::with_capacity(count);
    let mut derivatives = Vec::with_capacity(count);
    if count == 0 {
        return (values, derivatives);
    }
    let (first, first_derivative) = simplified_bessel_aux(distance, cutoff, 0);
    values.push(first);
    derivatives.push(first_derivative);
    let mut previous_d = 1.0;
    for order in 1..count {
        let n = order as f64;
        let e = n.powi(2) * (n + 2.0).powi(2) / (4.0 * (n + 1.0).powi(4) + 1.0);
        let d = 1.0 - e / previous_d;
        let mixing = (e / previous_d).sqrt();
        let normalization = d.sqrt();
        let (auxiliary, auxiliary_derivative) = simplified_bessel_aux(distance, cutoff, order);
        values.push((auxiliary + mixing * values[order - 1]) / normalization);
        derivatives.push((auxiliary_derivative + mixing * derivatives[order - 1]) / normalization);
        previous_d = d;
    }
    (values, derivatives)
}

fn descending_polynomial_switch(value: f64, upper: f64, width: f64) -> f64 {
    descending_polynomial_switch_with_derivative(value, upper, width).0
}

fn descending_polynomial_switch_with_derivative(value: f64, upper: f64, width: f64) -> (f64, f64) {
    if width <= 0.0 {
        return ((value < upper) as u8 as f64, 0.0);
    }
    if value <= upper - width {
        (1.0, 0.0)
    } else if value >= upper {
        (0.0, 0.0)
    } else {
        let x = 1.0 - 2.0 * (1.0 + (value - upper) / width);
        let switch = 0.5 + 3.75 * (x / 4.0 - x.powi(3) / 6.0 + x.powi(5) / 20.0);
        let derivative = -7.5 / width * (0.25 - x * x / 2.0 + x.powi(4) / 4.0);
        (switch, derivative)
    }
}

fn finnis_sinclair(value: f64, exponent: f64) -> f64 {
    const SCALE: f64 = 1.0e6;
    if value.abs() <= 1.0e-10 {
        return SCALE.powf(1.0 - exponent) * value;
    }
    let magnitude = value.abs();
    let scaled_cube = (SCALE * magnitude).powi(3);
    let blend = if scaled_cube > 30.0 {
        0.0
    } else {
        (-scaled_cube).exp()
    };
    value.signum()
        * ((1.0 - blend) * magnitude.powf(exponent)
            + SCALE.powf(1.0 - exponent) * blend * magnitude)
}

fn finnis_sinclair_shifted_scaled(value: f64, exponent: f64) -> f64 {
    if (exponent - 1.0).abs() < 1.0e-10 {
        return value;
    }
    let magnitude = value.abs();
    let inverse = 1.0 / exponent;
    let exponential = (-magnitude).exp();
    let x_offset = inverse.powf(inverse / (1.0 - inverse)) * exponential;
    let y_offset = inverse.powf(1.0 / (1.0 - inverse)) * exponential;
    value.signum() * ((x_offset + magnitude).powf(exponent) - y_offset)
}

fn spherical_harmonics(unit: [f64; 3], lmax: usize) -> Vec<Complex> {
    let width = 2 * lmax + 1;
    let mut result = vec![Complex::default(); (lmax + 1) * width];
    let z = unit[2].clamp(-1.0, 1.0);
    let azimuth = unit[1].atan2(unit[0]);
    for m in 0..=lmax {
        let mut p_mm = 1.0;
        if m > 0 {
            let transverse = (1.0 - z * z).max(0.0).sqrt();
            for factor in 1..=m {
                p_mm *= -((2 * factor - 1) as f64) * transverse;
            }
        }
        let mut values = vec![0.0; lmax + 1];
        values[m] = p_mm;
        if m < lmax {
            values[m + 1] = z * (2 * m + 1) as f64 * p_mm;
        }
        for l in (m + 2)..=lmax {
            values[l] = ((2 * l - 1) as f64 * z * values[l - 1]
                - (l + m - 1) as f64 * values[l - 2])
                / (l - m) as f64;
        }
        for l in m..=lmax {
            let normalization = ((2 * l + 1) as f64 * factorial(l - m) / factorial(l + m)).sqrt();
            let phase = m as f64 * azimuth;
            let positive = Complex {
                re: normalization * values[l] * phase.cos(),
                im: normalization * values[l] * phase.sin(),
            };
            result[lm_index(l, m as isize, lmax)] = positive;
            if m > 0 {
                let sign = if m % 2 == 0 { 1.0 } else { -1.0 };
                result[lm_index(l, -(m as isize), lmax)] = positive.conjugate() * sign;
            }
        }
    }
    result
}

fn factorial(value: usize) -> f64 {
    (1..=value).fold(1.0, |product, factor| product * factor as f64)
}

fn lm_index(l: usize, m: isize, lmax: usize) -> usize {
    l * (2 * lmax + 1) + (m + lmax as isize) as usize
}

fn parse_yace(text: &str, path: &Path) -> Result<PaceModel, String> {
    let records = flow_records(text);
    let elements_record = records
        .iter()
        .find(|record| record.trim_start().starts_with("elements:"))
        .ok_or_else(|| "PACE YAML is missing elements".to_string())?;
    let elements = parse_string_array(value_after_colon(elements_record)?)?;
    if elements.is_empty() {
        return Err("PACE model must contain at least one element".to_string());
    }
    let e0_record = records
        .iter()
        .find(|record| record.trim_start().starts_with("E0:"))
        .ok_or_else(|| "PACE YAML is missing E0".to_string())?;
    let e0 = parse_f64_array(value_after_colon(e0_record)?)?;
    if e0.len() != elements.len() {
        return Err("PACE E0 length must match elements".to_string());
    }
    let spline_record = records
        .iter()
        .find(|record| record.trim_start().starts_with("deltaSplineBins:"))
        .ok_or_else(|| "PACE YAML is missing deltaSplineBins".to_string())?;
    let spline_step = parse_f64(value_after_colon(spline_record)?, "deltaSplineBins")?;
    if !spline_step.is_finite() || spline_step <= 0.0 {
        return Err("PACE deltaSplineBins must be finite and positive".to_string());
    }

    let mut embeddings: Vec<Option<Embedding>> = vec![None; elements.len()];
    let mut bonds: Vec<Vec<Option<Bond>>> = vec![vec![None; elements.len()]; elements.len()];
    let mut functions = vec![Vec::new(); elements.len()];
    let mut section = "";
    let mut function_element = None;
    for record in &records {
        let trimmed = record.trim();
        match trimmed {
            "embeddings:" => {
                section = "embeddings";
                continue;
            }
            "bonds:" => {
                section = "bonds";
                continue;
            }
            "functions:" => {
                section = "functions";
                continue;
            }
            _ => {}
        }
        match section {
            "embeddings" if trimmed.contains('{') => {
                let (key, body) = split_key_body(trimmed)?;
                let index = parse_usize(key.trim(), "embedding element index")?;
                if index >= elements.len() {
                    return Err(format!(
                        "PACE embedding element index {index} is out of range"
                    ));
                }
                embeddings[index] = Some(parse_embedding(body)?);
            }
            "bonds" if trimmed.contains('{') => {
                let (key, body) = split_key_body(trimmed)?;
                let indices = parse_usize_array(key)?;
                if indices.len() != 2 || indices.iter().any(|index| *index >= elements.len()) {
                    return Err(format!("invalid PACE bond key {key}"));
                }
                bonds[indices[0]][indices[1]] = Some(parse_bond(body)?);
            }
            "functions" if trimmed.ends_with(':') && !trimmed.starts_with('-') => {
                let index = parse_usize(trimmed.trim_end_matches(':'), "function element index")?;
                if index >= elements.len() {
                    return Err(format!(
                        "PACE function element index {index} is out of range"
                    ));
                }
                function_element = Some(index);
            }
            "functions" if trimmed.starts_with('-') => {
                let element = function_element
                    .ok_or_else(|| "PACE function appears before its element block".to_string())?;
                let function = parse_function(trimmed.trim_start_matches('-').trim())?;
                if function.mu0 != element {
                    return Err(format!(
                        "PACE function mu0={} is stored in element block {element}",
                        function.mu0
                    ));
                }
                functions[element].push(function);
            }
            _ => {}
        }
    }

    let embeddings = embeddings
        .into_iter()
        .enumerate()
        .map(|(index, embedding)| {
            embedding.ok_or_else(|| format!("PACE model is missing embedding for element {index}"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let bonds = bonds
        .into_iter()
        .enumerate()
        .map(|(central, row)| {
            row.into_iter()
                .enumerate()
                .map(|(neighbor, bond)| {
                    bond.ok_or_else(|| {
                        format!("PACE model is missing directed bond [{central}, {neighbor}]")
                    })
                })
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()?;
    for (index, element_functions) in functions.iter().enumerate() {
        if element_functions.is_empty() {
            return Err(format!("PACE model has no functions for element {index}"));
        }
        for function in element_functions {
            if function.ndensity != embeddings[index].ndensity {
                return Err(format!(
                    "PACE function ndensity does not match embedding for element {index}"
                ));
            }
            if function.ctildes.iter().any(|value| !value.is_finite()) {
                return Err(format!(
                    "PACE function for element {index} has non-finite ctildes"
                ));
            }
            for factor in 0..function.rank {
                let neighbor = function.mus[factor];
                if neighbor >= elements.len() {
                    return Err(format!(
                        "PACE function for element {index} references species {neighbor}"
                    ));
                }
                let bond = &bonds[index][neighbor];
                let radial_limit = if function.rank == 1 {
                    bond.nradbasemax
                } else {
                    bond.nradmax
                };
                if function.ns[factor] > radial_limit || function.ls[factor] > bond.lmax {
                    return Err(format!(
                        "PACE function for element {index} exceeds bond [{index}, {neighbor}] basis dimensions"
                    ));
                }
                if function.rank == 1 && function.ls[factor] != 0 {
                    return Err("PACE rank-one functions must use l=0".to_string());
                }
            }
        }
    }
    if e0.iter().any(|value| !value.is_finite()) {
        return Err("PACE E0 values must be finite".to_string());
    }
    for (index, embedding) in embeddings.iter().enumerate() {
        if embedding.parameters.iter().any(|value| !value.is_finite())
            || !embedding.rho_core_cutoff.is_finite()
            || !embedding.drho_core_cutoff.is_finite()
            || embedding.drho_core_cutoff < 0.0
        {
            return Err(format!(
                "PACE embedding {index} contains invalid numeric values"
            ));
        }
    }
    for (central, row) in bonds.iter().enumerate() {
        for (neighbor, bond) in row.iter().enumerate() {
            if bond.nradmax == 0
                || bond.nradbasemax == 0
                || !bond.lambda.is_finite()
                || bond.lambda <= 0.0
                || !bond.cutoff_width.is_finite()
                || bond.cutoff_width < 0.0
                || !bond.inner_cutoff.is_finite()
                || !bond.inner_width.is_finite()
                || bond.inner_width < 0.0
                || !bond.prehc.is_finite()
                || !bond.lambdahc.is_finite()
                || bond.coefficients.iter().any(|value| !value.is_finite())
            {
                return Err(format!(
                    "PACE bond [{central}, {neighbor}] contains invalid dimensions or numeric values"
                ));
            }
        }
    }
    let maximum_cutoff = bonds
        .iter()
        .flatten()
        .map(|bond| bond.cutoff)
        .fold(0.0, f64::max);
    Ok(PaceModel {
        elements,
        source_path: path.to_path_buf(),
        e0,
        embeddings,
        bonds,
        functions,
        type_to_element: Vec::new(),
        maximum_cutoff,
        spline_step,
    })
}

fn parse_embedding(body: &str) -> Result<Embedding, String> {
    let fields = parse_flow_map(body)?;
    let ndensity = field_usize(&fields, "ndensity")?;
    let parameters = field_f64_array(&fields, "FS_parameters")?;
    if parameters.len() != 2 * ndensity {
        return Err("PACE FS_parameters must contain two values per density".to_string());
    }
    let kind = match field(&fields, "npoti")? {
        "FinnisSinclair" => EmbeddingKind::FinnisSinclair,
        "FinnisSinclairShiftedScaled" => EmbeddingKind::FinnisSinclairShiftedScaled,
        value => return Err(format!("unsupported PACE embedding {value}")),
    };
    Ok(Embedding {
        ndensity,
        parameters,
        kind,
        rho_core_cutoff: field_f64(&fields, "rho_core_cutoff")?,
        drho_core_cutoff: field_f64(&fields, "drho_core_cutoff")?,
    })
}

fn parse_bond(body: &str) -> Result<Bond, String> {
    let fields = parse_flow_map(body)?;
    let nradmax = field_usize(&fields, "nradmax")?;
    let lmax = field_usize(&fields, "lmax")?;
    let nradbasemax = field_usize(&fields, "nradbasemax")?;
    let radial_kind = match field(&fields, "radbasename")? {
        "ChebPow" => RadialKind::PowerScaled,
        "ChebExpCos" => RadialKind::ExponentialCosine,
        "ChebLinear" => RadialKind::LinearScaled,
        "SBessel" => RadialKind::SimplifiedBessel,
        value => return Err(format!("unsupported PACE radial basis {value}")),
    };
    let radial_parameters = field_f64_array(&fields, "radparameters")?;
    let lambda = *radial_parameters
        .first()
        .ok_or_else(|| "PACE radparameters must not be empty".to_string())?;
    let coefficients = field_f64_array(&fields, "radcoefficients")?;
    let expected = nradmax * (lmax + 1) * nradbasemax;
    if coefficients.len() != expected {
        return Err(format!(
            "PACE radcoefficients has {} values, expected {expected}",
            coefficients.len()
        ));
    }
    let inner_kind = match field(&fields, "inner_cutoff_type")? {
        "density" => InnerCutoffKind::Density,
        "distance" => InnerCutoffKind::Distance,
        "zbl" => return Err("PACE zbl inner cutoff is not supported by pace_native".to_string()),
        value => return Err(format!("unsupported PACE inner_cutoff_type {value}")),
    };
    let bond = Bond {
        nradmax,
        lmax,
        nradbasemax,
        radial_kind,
        lambda,
        coefficients,
        prehc: field_f64(&fields, "prehc")?,
        lambdahc: field_f64(&fields, "lambdahc")?,
        cutoff: field_f64(&fields, "rcut")?,
        cutoff_width: field_f64(&fields, "dcut")?,
        inner_cutoff: field_f64(&fields, "rcut_in")?,
        inner_width: field_f64(&fields, "dcut_in")?,
        inner_kind,
    };
    if bond.cutoff <= 0.0 || !bond.cutoff.is_finite() {
        return Err("PACE bond cutoff must be finite and positive".to_string());
    }
    Ok(bond)
}

fn parse_function(body: &str) -> Result<BasisFunction, String> {
    let fields = parse_flow_map(body)?;
    let mu0 = field_usize(&fields, "mu0")?;
    let rank = field_usize(&fields, "rank")?;
    let ndensity = field_usize(&fields, "ndensity")?;
    let num_ms_combs = field_usize(&fields, "num_ms_combs")?;
    let mus = field_usize_array(&fields, "mus")?;
    let ns = field_usize_array(&fields, "ns")?;
    let ls = field_usize_array(&fields, "ls")?;
    let ms_combs = field_isize_array(&fields, "ms_combs")?;
    let ctildes = field_f64_array(&fields, "ctildes")?;
    if rank == 0 || mus.len() != rank || ns.len() != rank || ls.len() != rank {
        return Err("PACE function rank does not match mus/ns/ls".to_string());
    }
    if ns.contains(&0) {
        return Err("PACE radial indices are one-based and must be positive".to_string());
    }
    if ms_combs.len() != rank * num_ms_combs {
        return Err("PACE ms_combs length does not match rank*num_ms_combs".to_string());
    }
    if ctildes.len() != ndensity * num_ms_combs {
        return Err("PACE ctildes length does not match ndensity*num_ms_combs".to_string());
    }
    for combination in ms_combs.chunks_exact(rank) {
        for ((m, l), _factor) in combination.iter().zip(&ls).zip(0..) {
            if m.unsigned_abs() > *l {
                return Err("PACE magnetic index exceeds angular index".to_string());
            }
        }
    }
    Ok(BasisFunction {
        mu0,
        rank,
        ndensity,
        num_ms_combs,
        mus,
        ns,
        ls,
        ms_combs,
        ctildes,
    })
}

fn flow_records(text: &str) -> Vec<String> {
    let mut records = Vec::new();
    let mut current = String::new();
    let mut depth = 0isize;
    for line in text.lines() {
        let line = line.split('#').next().unwrap_or("").trim_end();
        if line.trim().is_empty() {
            continue;
        }
        if current.is_empty() {
            current.push_str(line);
        } else {
            current.push(' ');
            current.push_str(line.trim());
        }
        depth += line
            .chars()
            .map(|character| match character {
                '[' | '{' => 1,
                ']' | '}' => -1,
                _ => 0,
            })
            .sum::<isize>();
        if depth == 0 {
            records.push(std::mem::take(&mut current));
        }
    }
    if !current.is_empty() {
        records.push(current);
    }
    records
}

fn parse_flow_map(body: &str) -> Result<HashMap<String, String>, String> {
    let body = body.trim();
    let body = body
        .strip_prefix('{')
        .and_then(|value| value.strip_suffix('}'))
        .ok_or_else(|| format!("expected flow mapping, got {body}"))?;
    let mut fields = HashMap::new();
    for entry in split_top_level(body, ',') {
        let (key, value) = entry
            .split_once(':')
            .ok_or_else(|| format!("invalid PACE mapping entry {entry}"))?;
        fields.insert(key.trim().to_string(), value.trim().to_string());
    }
    Ok(fields)
}

fn split_top_level(text: &str, separator: char) -> Vec<&str> {
    let mut result = Vec::new();
    let mut depth = 0isize;
    let mut start = 0;
    for (index, character) in text.char_indices() {
        match character {
            '[' | '{' => depth += 1,
            ']' | '}' => depth -= 1,
            _ => {}
        }
        if character == separator && depth == 0 {
            result.push(text[start..index].trim());
            start = index + character.len_utf8();
        }
    }
    result.push(text[start..].trim());
    result
}

fn split_key_body(text: &str) -> Result<(&str, &str), String> {
    let brace = text
        .find('{')
        .ok_or_else(|| format!("expected mapping body in {text}"))?;
    let prefix = text[..brace].trim();
    let key = prefix
        .strip_suffix(':')
        .ok_or_else(|| format!("expected mapping key in {text}"))?;
    Ok((key.trim(), &text[brace..]))
}

fn value_after_colon(text: &str) -> Result<&str, String> {
    text.split_once(':')
        .map(|(_, value)| value.trim())
        .ok_or_else(|| format!("expected key/value record {text}"))
}

fn field<'a>(fields: &'a HashMap<String, String>, name: &str) -> Result<&'a str, String> {
    fields
        .get(name)
        .map(String::as_str)
        .ok_or_else(|| format!("PACE mapping is missing {name}"))
}

fn field_f64(fields: &HashMap<String, String>, name: &str) -> Result<f64, String> {
    parse_f64(field(fields, name)?, name)
}

fn field_usize(fields: &HashMap<String, String>, name: &str) -> Result<usize, String> {
    parse_usize(field(fields, name)?, name)
}

fn field_f64_array(fields: &HashMap<String, String>, name: &str) -> Result<Vec<f64>, String> {
    parse_f64_array(field(fields, name)?)
}

fn field_usize_array(fields: &HashMap<String, String>, name: &str) -> Result<Vec<usize>, String> {
    parse_usize_array(field(fields, name)?)
}

fn field_isize_array(fields: &HashMap<String, String>, name: &str) -> Result<Vec<isize>, String> {
    numeric_tokens(field(fields, name)?)
        .into_iter()
        .map(|token| {
            token
                .parse::<isize>()
                .map_err(|error| format!("invalid integer {token} in {name}: {error}"))
        })
        .collect()
}

fn parse_string_array(text: &str) -> Result<Vec<String>, String> {
    let body = text
        .trim()
        .strip_prefix('[')
        .and_then(|value| value.strip_suffix(']'))
        .ok_or_else(|| format!("expected string array, got {text}"))?;
    Ok(split_top_level(body, ',')
        .into_iter()
        .map(|value| value.trim().trim_matches(['\'', '"']).to_string())
        .collect())
}

fn parse_f64_array(text: &str) -> Result<Vec<f64>, String> {
    numeric_tokens(text)
        .into_iter()
        .map(|token| parse_f64(token, "numeric array"))
        .collect()
}

fn parse_usize_array(text: &str) -> Result<Vec<usize>, String> {
    numeric_tokens(text)
        .into_iter()
        .map(|token| parse_usize(token, "integer array"))
        .collect()
}

fn numeric_tokens(text: &str) -> Vec<&str> {
    text.split(|character: char| {
        character == '[' || character == ']' || character == ',' || character.is_whitespace()
    })
    .filter(|token| !token.is_empty())
    .collect()
}

fn parse_f64(text: &str, field_name: &str) -> Result<f64, String> {
    text.trim()
        .parse::<f64>()
        .map_err(|error| format!("invalid PACE {field_name} value {text}: {error}"))
}

fn parse_usize(text: &str, field_name: &str) -> Result<usize, String> {
    text.trim()
        .parse::<usize>()
        .map_err(|error| format!("invalid PACE {field_name} value {text}: {error}"))
}

fn resolve_type_elements(
    model_elements: &[String],
    type_elements: &[String],
) -> Result<Vec<usize>, String> {
    type_elements
        .iter()
        .map(|element| {
            model_elements
                .iter()
                .position(|candidate| candidate == element)
                .ok_or_else(|| {
                    format!(
                        "rsmith species {element} is not present in PACE elements [{}]",
                        model_elements.join(", ")
                    )
                })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scaled_spherical_harmonics_have_y00_equal_one() {
        let harmonics = spherical_harmonics([1.0, 0.0, 0.0], 2);
        assert!((harmonics[lm_index(0, 0, 2)].re - 1.0).abs() < 1.0e-15);
        assert!((harmonics[lm_index(1, 0, 2)].re).abs() < 1.0e-15);
        assert!((harmonics[lm_index(1, 1, 2)].re + (1.5f64).sqrt()).abs() < 1.0e-15);
    }

    #[test]
    fn chebpow_matches_closed_form_first_basis() {
        for distance in [0.5, 1.0, 2.0, 4.25] {
            let (values, _) = cheb_pow(3, 2.0, 5.0, distance);
            let expected = (1.0 - distance / 5.0_f64).powi(2);
            assert!((values[0] - expected).abs() < 1.0e-14);
        }
    }

    #[test]
    fn polynomial_switch_has_flat_endpoints() {
        assert_eq!(descending_polynomial_switch(1.0, 2.0, 0.5), 1.0);
        assert_eq!(descending_polynomial_switch(2.0, 2.0, 0.5), 0.0);
        let (middle, _) = descending_polynomial_switch_with_derivative(1.75, 2.0, 0.5);
        assert!((middle - 0.5).abs() < 1.0e-14);
    }

    #[test]
    fn simplified_bessel_derivatives_match_finite_differences() {
        let distance = 2.137;
        let step = 1.0e-6;
        let (values, derivatives) = simplified_bessel(15, 5.0, distance);
        let (lower, _) = simplified_bessel(15, 5.0, distance - step);
        let (upper, _) = simplified_bessel(15, 5.0, distance + step);
        for index in 0..values.len() {
            let numerical = (upper[index] - lower[index]) / (2.0 * step);
            assert!(
                (derivatives[index] - numerical).abs() < 2.0e-9,
                "basis {index}: analytic={}, numerical={numerical}",
                derivatives[index]
            );
        }
    }
}
