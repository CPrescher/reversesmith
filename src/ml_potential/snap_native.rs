//! Native SNAP model-file support.
//!
//! The evaluator is intentionally independent of LAMMPS.  LAMMPS is used only
//! in integration tests as a numerical reference for FitSNAP-compatible files.

mod math;

use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};

/// Parameters controlling the SNAP bispectrum calculation.
#[derive(Debug, Clone, PartialEq)]
pub struct SnapParameters {
    pub rcutfac: f64,
    pub twojmax: usize,
    pub rfac0: f64,
    pub rmin0: f64,
    pub switchflag: bool,
    pub bzeroflag: bool,
    pub quadraticflag: bool,
    pub chemflag: bool,
    pub bnormflag: bool,
    pub wselfallflag: bool,
    pub switchinnerflag: bool,
    pub sinner: Vec<f64>,
    pub dinner: Vec<f64>,
}

impl Default for SnapParameters {
    fn default() -> Self {
        Self {
            rcutfac: 0.0,
            twojmax: 0,
            rfac0: 0.99363,
            rmin0: 0.0,
            switchflag: true,
            bzeroflag: true,
            quadraticflag: false,
            chemflag: false,
            bnormflag: false,
            wselfallflag: false,
            switchinnerflag: false,
            sinner: Vec::new(),
            dinner: Vec::new(),
        }
    }
}

/// Coefficients and element metadata from a `.snapcoeff` file.
#[derive(Debug, Clone, PartialEq)]
pub struct SnapCoefficients {
    pub ncoeff: usize,
    pub elements: Vec<SnapElement>,
    pub units: Option<String>,
}

/// One element block from a SNAP coefficient file.
#[derive(Debug, Clone, PartialEq)]
pub struct SnapElement {
    pub name: String,
    pub radius: f64,
    pub weight: f64,
    pub coefficients: Vec<f64>,
}

/// One neighbor of a central atom for native SNAP evaluation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SnapNeighbor {
    /// Minimum-image displacement from the central atom, in Angstrom.
    pub displacement: [f64; 3],
    /// Zero-based rsmith atom-type index, resolved through `type_to_element`.
    pub type_index: usize,
}

/// Fully resolved pair of FitSNAP output files.
#[derive(Debug, Clone)]
pub struct SnapModelFiles {
    pub parameters: SnapParameters,
    pub coefficients: SnapCoefficients,
    pub type_to_element: Vec<usize>,
    pub coefficient_path: PathBuf,
    pub parameter_path: PathBuf,
    basis: math::SnapBasis,
}

impl SnapModelFiles {
    pub fn load(
        coefficient_path: &Path,
        parameter_path: &Path,
        type_elements: &[String],
    ) -> Result<Self, String> {
        let coefficient_text = fs::read_to_string(coefficient_path).map_err(|error| {
            format!(
                "failed to read SNAP coefficient file {}: {error}",
                coefficient_path.display()
            )
        })?;
        let parameter_text = fs::read_to_string(parameter_path).map_err(|error| {
            format!(
                "failed to read SNAP parameter file {}: {error}",
                parameter_path.display()
            )
        })?;

        let coefficients = parse_coefficients(&coefficient_text)?;
        let parameters = parse_parameters(&parameter_text)?;
        let type_to_element = resolve_type_elements(&coefficients, type_elements)?;
        let basis = math::SnapBasis::new(parameters.twojmax)?;

        validate_model(&parameters, &coefficients, basis.descriptor_count())?;

        Ok(Self {
            parameters,
            coefficients,
            type_to_element,
            coefficient_path: coefficient_path.to_path_buf(),
            parameter_path: parameter_path.to_path_buf(),
            basis,
        })
    }

    pub fn descriptor_count(&self) -> usize {
        self.basis.descriptor_count()
    }

    /// Evaluate the standard, non-chemical SNAP descriptors of one atom.
    pub fn atomic_descriptors(
        &self,
        central_type_index: usize,
        neighbors: &[SnapNeighbor],
    ) -> Result<Vec<f64>, String> {
        if self.parameters.chemflag {
            return Err("native chemical SNAP descriptors are not implemented yet".to_string());
        }
        let central_element_index = self.element_index_for_type(central_type_index)?;
        let central_element = &self.coefficients.elements[central_element_index];
        let mut density = self.basis.empty_density();

        for neighbor in neighbors {
            let neighbor_element_index = self.element_index_for_type(neighbor.type_index)?;
            let neighbor_element = &self.coefficients.elements[neighbor_element_index];
            if neighbor
                .displacement
                .iter()
                .any(|coordinate| !coordinate.is_finite())
            {
                return Err("SNAP neighbor displacement must be finite".to_string());
            }
            let radius_squared = neighbor
                .displacement
                .iter()
                .map(|coordinate| coordinate * coordinate)
                .sum::<f64>();
            if !radius_squared.is_finite() {
                return Err("SNAP neighbor distance must be finite".to_string());
            }
            if radius_squared == 0.0 {
                return Err("SNAP neighbor displacement must be nonzero".to_string());
            }
            let radius = radius_squared.sqrt();
            let cutoff =
                self.parameters.rcutfac * (central_element.radius + neighbor_element.radius);
            if radius > cutoff {
                continue;
            }

            let (inner_midpoint, inner_half_width) = if self.parameters.switchinnerflag {
                (
                    0.5 * (self.parameters.sinner[central_element_index]
                        + self.parameters.sinner[neighbor_element_index]),
                    0.5 * (self.parameters.dinner[central_element_index]
                        + self.parameters.dinner[neighbor_element_index]),
                )
            } else {
                (0.0, 0.0)
            };
            let switching_weight = radial_switch(
                radius,
                cutoff,
                self.parameters.rmin0,
                self.parameters.switchflag,
                self.parameters.switchinnerflag,
                inner_midpoint,
                inner_half_width,
            );
            if switching_weight == 0.0 {
                continue;
            }
            let theta =
                (radius - self.parameters.rmin0) * self.parameters.rfac0 * std::f64::consts::PI
                    / (cutoff - self.parameters.rmin0);
            self.basis.add_neighbor(
                &mut density,
                neighbor.displacement,
                theta,
                switching_weight * neighbor_element.weight,
            );
        }

        Ok(self.basis.bispectrum(
            &density,
            self.parameters.bnormflag,
            self.parameters.bzeroflag,
        ))
    }

    /// Evaluate one atom's linear SNAP energy in eV.
    pub fn atomic_energy(
        &self,
        central_type_index: usize,
        neighbors: &[SnapNeighbor],
    ) -> Result<f64, String> {
        if self.parameters.quadraticflag {
            return Err("native quadratic SNAP energies are not implemented yet".to_string());
        }
        let element_index = self.element_index_for_type(central_type_index)?;
        let descriptors = self.atomic_descriptors(central_type_index, neighbors)?;
        let coefficients = &self.coefficients.elements[element_index].coefficients;
        Ok(coefficients[0]
            + coefficients[1..]
                .iter()
                .zip(descriptors)
                .map(|(coefficient, descriptor)| coefficient * descriptor)
                .sum::<f64>())
    }

    fn element_index_for_type(&self, type_index: usize) -> Result<usize, String> {
        self.type_to_element
            .get(type_index)
            .copied()
            .ok_or_else(|| {
                format!(
                "SNAP atom type index {type_index} is outside the configured mapping of {} types",
                self.type_to_element.len()
            )
            })
    }
}

#[allow(clippy::too_many_arguments)]
fn radial_switch(
    radius: f64,
    cutoff: f64,
    minimum_radius: f64,
    outer_enabled: bool,
    inner_enabled: bool,
    inner_midpoint: f64,
    inner_half_width: f64,
) -> f64 {
    let mut value = if !outer_enabled || radius <= minimum_radius {
        1.0
    } else if radius > cutoff {
        0.0
    } else {
        0.5 * ((std::f64::consts::PI * (radius - minimum_radius) / (cutoff - minimum_radius)).cos()
            + 1.0)
    };

    if inner_enabled && radius < inner_midpoint + inner_half_width {
        if radius > inner_midpoint - inner_half_width {
            value *= 0.5
                * (1.0
                    - (std::f64::consts::FRAC_PI_2
                        + std::f64::consts::FRAC_PI_2 * (radius - inner_midpoint)
                            / inner_half_width)
                        .cos());
        } else {
            value = 0.0;
        }
    }

    value
}

pub fn parse_parameters(input: &str) -> Result<SnapParameters, String> {
    let mut parameters = SnapParameters::default();
    let mut seen = HashSet::new();

    for (line_index, raw_line) in input.lines().enumerate() {
        let line_number = line_index + 1;
        let line = strip_comment(raw_line);
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        let key = fields[0];
        if !seen.insert(key.to_string()) {
            return Err(format!(
                "duplicate SNAP parameter '{key}' on line {line_number}"
            ));
        }

        match key {
            "rcutfac" => parameters.rcutfac = parse_one(fields.as_slice(), line_number)?,
            "twojmax" => parameters.twojmax = parse_one(fields.as_slice(), line_number)?,
            "rfac0" => parameters.rfac0 = parse_one(fields.as_slice(), line_number)?,
            "rmin0" => parameters.rmin0 = parse_one(fields.as_slice(), line_number)?,
            "switchflag" => parameters.switchflag = parse_flag(fields.as_slice(), line_number)?,
            "bzeroflag" => parameters.bzeroflag = parse_flag(fields.as_slice(), line_number)?,
            "quadraticflag" => {
                parameters.quadraticflag = parse_flag(fields.as_slice(), line_number)?
            }
            "chemflag" => parameters.chemflag = parse_flag(fields.as_slice(), line_number)?,
            "bnormflag" => parameters.bnormflag = parse_flag(fields.as_slice(), line_number)?,
            "wselfallflag" => parameters.wselfallflag = parse_flag(fields.as_slice(), line_number)?,
            "switchinnerflag" => {
                parameters.switchinnerflag = parse_flag(fields.as_slice(), line_number)?
            }
            "sinner" => parameters.sinner = parse_many(fields.as_slice(), line_number)?,
            "dinner" => parameters.dinner = parse_many(fields.as_slice(), line_number)?,
            // These tune LAMMPS execution and do not change the mathematical model.
            "chunksize" | "parallelthresh" => {
                let _: usize = parse_one(fields.as_slice(), line_number)?;
            }
            _ => {
                return Err(format!(
                    "unsupported SNAP parameter '{key}' on line {line_number}"
                ))
            }
        }
    }

    if !seen.contains("rcutfac") {
        return Err("SNAP parameter file is missing required 'rcutfac'".to_string());
    }
    if !seen.contains("twojmax") {
        return Err("SNAP parameter file is missing required 'twojmax'".to_string());
    }
    if parameters.rcutfac <= 0.0 || !parameters.rcutfac.is_finite() {
        return Err("SNAP rcutfac must be finite and greater than zero".to_string());
    }
    if !parameters.rfac0.is_finite() || parameters.rfac0 <= 0.0 || parameters.rfac0 >= 1.0 {
        return Err("SNAP rfac0 must be finite, greater than 0, and less than 1".to_string());
    }
    if !parameters.rmin0.is_finite() || parameters.rmin0 < 0.0 {
        return Err("SNAP rmin0 must be finite and greater than or equal to 0".to_string());
    }
    if parameters.switchinnerflag {
        if parameters.sinner.is_empty() || parameters.dinner.is_empty() {
            return Err("SNAP inner switching requires both sinner and dinner".to_string());
        }
        if parameters.sinner.len() != parameters.dinner.len() {
            return Err(
                "SNAP sinner and dinner must contain the same number of values".to_string(),
            );
        }
        if parameters
            .sinner
            .iter()
            .any(|value| *value <= 0.0 || !value.is_finite())
        {
            return Err("SNAP sinner values must be finite and greater than zero".to_string());
        }
        if parameters
            .dinner
            .iter()
            .any(|value| *value <= 0.0 || !value.is_finite())
        {
            return Err("SNAP dinner values must be finite and greater than zero".to_string());
        }
    }

    Ok(parameters)
}

pub fn parse_coefficients(input: &str) -> Result<SnapCoefficients, String> {
    let units = input
        .lines()
        .find_map(|line| metadata_value(line, "UNITS:"));
    let lines: Vec<(usize, &str)> = input
        .lines()
        .enumerate()
        .filter_map(|(index, line)| {
            let line = strip_comment(line);
            (!line.is_empty()).then_some((index + 1, line))
        })
        .collect();
    let Some(&(header_line, header)) = lines.first() else {
        return Err("SNAP coefficient file is empty".to_string());
    };
    let header_fields: Vec<&str> = header.split_whitespace().collect();
    if header_fields.len() != 2 {
        return Err(format!(
            "SNAP coefficient header on line {header_line} must contain nelem and ncoeff"
        ));
    }
    let nelem = parse_value::<usize>(header_fields[0], header_line, "nelem")?;
    let ncoeff = parse_value::<usize>(header_fields[1], header_line, "ncoeff")?;
    if nelem == 0 || ncoeff == 0 {
        return Err("SNAP nelem and ncoeff must be greater than zero".to_string());
    }

    let mut cursor = 1;
    let mut elements = Vec::with_capacity(nelem);
    let mut names = HashSet::new();
    for element_index in 0..nelem {
        let Some(&(line_number, element_line)) = lines.get(cursor) else {
            return Err(format!(
                "SNAP coefficient file ended before element block {}",
                element_index + 1
            ));
        };
        cursor += 1;
        let fields: Vec<&str> = element_line.split_whitespace().collect();
        if fields.len() != 3 {
            return Err(format!(
                "SNAP element header on line {line_number} must contain name, radius, and weight"
            ));
        }
        let name = fields[0].to_string();
        if !names.insert(name.clone()) {
            return Err(format!(
                "duplicate SNAP element '{name}' on line {line_number}"
            ));
        }
        let radius = parse_value::<f64>(fields[1], line_number, "radius")?;
        let weight = parse_value::<f64>(fields[2], line_number, "weight")?;
        if radius <= 0.0 || !radius.is_finite() {
            return Err(format!(
                "SNAP element radius on line {line_number} must be finite and greater than zero"
            ));
        }
        if !weight.is_finite() {
            return Err(format!(
                "SNAP element weight on line {line_number} must be finite"
            ));
        }

        let mut coefficients = Vec::with_capacity(ncoeff);
        for coefficient_index in 0..ncoeff {
            let Some(&(coefficient_line, value_line)) = lines.get(cursor) else {
                return Err(format!(
                    "SNAP element '{name}' ended after {coefficient_index} of {ncoeff} coefficients"
                ));
            };
            cursor += 1;
            let fields: Vec<&str> = value_line.split_whitespace().collect();
            if fields.len() != 1 {
                return Err(format!(
                    "SNAP coefficient on line {coefficient_line} must contain one value"
                ));
            }
            let coefficient = parse_value::<f64>(fields[0], coefficient_line, "coefficient")?;
            if !coefficient.is_finite() {
                return Err(format!(
                    "SNAP coefficient on line {coefficient_line} must be finite"
                ));
            }
            coefficients.push(coefficient);
        }
        elements.push(SnapElement {
            name,
            radius,
            weight,
            coefficients,
        });
    }
    if let Some((line_number, _)) = lines.get(cursor) {
        return Err(format!(
            "unexpected data after SNAP element blocks on line {line_number}"
        ));
    }

    Ok(SnapCoefficients {
        ncoeff,
        elements,
        units,
    })
}

fn resolve_type_elements(
    coefficients: &SnapCoefficients,
    type_elements: &[String],
) -> Result<Vec<usize>, String> {
    if type_elements.is_empty() {
        return Err("SNAP element mapping must not be empty".to_string());
    }
    let indices: HashMap<&str, usize> = coefficients
        .elements
        .iter()
        .enumerate()
        .map(|(index, element)| (element.name.as_str(), index))
        .collect();
    type_elements
        .iter()
        .map(|name| {
            indices.get(name.as_str()).copied().ok_or_else(|| {
                format!("SNAP element mapping refers to '{name}', which is absent from the coefficient file")
            })
        })
        .collect()
}

fn validate_model(
    parameters: &SnapParameters,
    coefficients: &SnapCoefficients,
    descriptor_count: usize,
) -> Result<(), String> {
    if let Some(units) = &coefficients.units {
        if units != "metal" {
            return Err(format!(
                "SNAP model uses '{units}' units; rsmith requires LAMMPS 'metal' units"
            ));
        }
    }
    if parameters.switchinnerflag && parameters.sinner.len() != coefficients.elements.len() {
        return Err(format!(
            "SNAP inner switching supplies {} values for {} coefficient elements",
            parameters.sinner.len(),
            coefficients.elements.len()
        ));
    }
    let minimum_cutoff = coefficients
        .elements
        .iter()
        .flat_map(|central| {
            coefficients
                .elements
                .iter()
                .map(move |neighbor| parameters.rcutfac * (central.radius + neighbor.radius))
        })
        .fold(f64::INFINITY, f64::min);
    if parameters.rmin0 >= minimum_cutoff {
        return Err(format!(
            "SNAP rmin0={} must be smaller than every pair cutoff; minimum cutoff is {minimum_cutoff}",
            parameters.rmin0
        ));
    }
    let chemical_channels = if parameters.chemflag {
        coefficients
            .elements
            .len()
            .checked_pow(3)
            .ok_or_else(|| "SNAP chemical descriptor count overflows usize".to_string())?
    } else {
        1
    };
    let linear_count = descriptor_count
        .checked_mul(chemical_channels)
        .ok_or_else(|| "SNAP linear descriptor count overflows usize".to_string())?;
    let quadratic_count = if parameters.quadraticflag {
        linear_count
            .checked_add(1)
            .and_then(|next| linear_count.checked_mul(next))
            .map(|count| count / 2)
            .ok_or_else(|| "SNAP quadratic descriptor count overflows usize".to_string())?
    } else {
        0
    };
    let expected_ncoeff = 1usize
        .checked_add(linear_count)
        .and_then(|count| count.checked_add(quadratic_count))
        .ok_or_else(|| "SNAP coefficient count overflows usize".to_string())?;
    if coefficients.ncoeff != expected_ncoeff {
        return Err(format!(
            "SNAP coefficient file has {} coefficients per element, but twojmax={} and the enabled model flags require {expected_ncoeff}",
            coefficients.ncoeff, parameters.twojmax
        ));
    }
    Ok(())
}

fn parse_flag(fields: &[&str], line_number: usize) -> Result<bool, String> {
    let value: usize = parse_one(fields, line_number)?;
    match value {
        0 => Ok(false),
        1 => Ok(true),
        _ => Err(format!("SNAP flag on line {line_number} must be 0 or 1")),
    }
}

fn parse_one<T: std::str::FromStr>(fields: &[&str], line_number: usize) -> Result<T, String> {
    if fields.len() != 2 {
        return Err(format!(
            "SNAP parameter '{}' on line {line_number} requires exactly one value",
            fields[0]
        ));
    }
    parse_value(fields[1], line_number, fields[0])
}

fn parse_many<T: std::str::FromStr>(fields: &[&str], line_number: usize) -> Result<Vec<T>, String> {
    if fields.len() < 2 {
        return Err(format!(
            "SNAP parameter '{}' on line {line_number} requires at least one value",
            fields[0]
        ));
    }
    fields[1..]
        .iter()
        .map(|value| parse_value(value, line_number, fields[0]))
        .collect()
}

fn parse_value<T: std::str::FromStr>(
    value: &str,
    line_number: usize,
    label: &str,
) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("invalid {label} value '{value}' on line {line_number}"))
}

fn strip_comment(line: &str) -> &str {
    line.split('#').next().unwrap_or_default().trim()
}

fn metadata_value(line: &str, key: &str) -> Option<String> {
    let (_, tail) = line.split_once(key)?;
    tail.split_whitespace().next().map(ToOwned::to_owned)
}

#[cfg(test)]
mod tests {
    use super::{parse_coefficients, parse_parameters, validate_model};

    #[test]
    fn parameter_defaults_match_lammps_snap_defaults() {
        let parameters = parse_parameters("rcutfac 4.9\ntwojmax 8\n").unwrap();
        assert_eq!(parameters.rfac0, 0.99363);
        assert!(parameters.switchflag);
        assert!(parameters.bzeroflag);
        assert!(!parameters.quadraticflag);
    }

    #[test]
    fn coefficient_parser_reads_metadata_and_multiple_elements() {
        let coefficients = parse_coefficients(
            "# DATE: 2026-08-01 UNITS: metal\n\
             2 2\n\
             Si 0.5 1.0\n\
             -1.0\n0.25\n\
             O 0.4 0.75\n\
             -2.0\n0.5\n",
        )
        .unwrap();
        assert_eq!(coefficients.units.as_deref(), Some("metal"));
        assert_eq!(coefficients.elements.len(), 2);
        assert_eq!(coefficients.elements[1].name, "O");
        assert_eq!(coefficients.elements[1].coefficients, vec![-2.0, 0.5]);
    }

    #[test]
    fn unknown_parameter_is_rejected() {
        let error = parse_parameters("rcutfac 4.9\ntwojmax 8\nmystery 1\n").unwrap_err();
        assert!(error.contains("unsupported SNAP parameter 'mystery'"));
    }

    #[test]
    fn coefficient_count_is_enforced_per_element() {
        let error = parse_coefficients("1 2\nSi 0.5 1.0\n-1.0\n").unwrap_err();
        assert!(error.contains("ended after 1 of 2 coefficients"));
    }

    #[test]
    fn coefficient_count_must_match_the_descriptor_basis() {
        let parameters = parse_parameters("rcutfac 4.0\ntwojmax 2\n").unwrap();
        let coefficients =
            parse_coefficients("1 5\nSi 0.5 1.0\n-1.0\n0.1\n0.2\n0.3\n0.4\n").unwrap();

        let error = validate_model(&parameters, &coefficients, 5).unwrap_err();
        assert!(error.contains("has 5 coefficients per element"), "{error}");
        assert!(error.contains("require 6"), "{error}");
    }

    #[test]
    fn non_finite_coefficients_are_rejected() {
        let error = parse_coefficients("1 1\nSi 0.5 1.0\nNaN\n").unwrap_err();
        assert!(error.contains("coefficient on line 3 must be finite"));
    }

    #[test]
    fn radial_mapping_parameters_are_strictly_validated() {
        for rfac0 in ["0", "1", "NaN", "inf"] {
            let input = format!("rcutfac 4.9\ntwojmax 8\nrfac0 {rfac0}\n");
            let error = parse_parameters(&input).unwrap_err();
            assert!(error.contains("rfac0 must be finite"), "{error}");
        }

        for rmin0 in ["-0.1", "NaN", "inf"] {
            let input = format!("rcutfac 4.9\ntwojmax 8\nrmin0 {rmin0}\n");
            let error = parse_parameters(&input).unwrap_err();
            assert!(error.contains("rmin0 must be finite"), "{error}");
        }
    }

    #[test]
    fn inner_switch_values_must_be_finite_and_positive() {
        for sinner in ["0", "-0.1", "NaN", "inf"] {
            let input =
                format!("rcutfac 4.9\ntwojmax 8\nswitchinnerflag 1\nsinner {sinner}\ndinner 0.2\n");
            let error = parse_parameters(&input).unwrap_err();
            assert!(error.contains("sinner values must be finite"), "{error}");
        }
    }
}
