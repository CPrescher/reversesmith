# Supported Elements

## X-ray form factors

Cromer-Mann parameterization of atomic form factors f(Q):

| Element | Supported |
|---------|-----------|
| Ca | Yes |
| Si | Yes |
| O  | Yes |
| Na | Yes |
| Al | Yes |
| Mg | Yes |

## Neutron scattering lengths

Coherent neutron scattering lengths b (fm):

| Element | b (fm) |
|---------|--------|
| Ca | 4.70 |
| Si | 4.149 |
| O  | 5.803 |
| Na | 3.63 |
| Al | 3.449 |
| Mg | 5.375 |

## Adding new elements

Additional elements can be added in `src/xray.rs` by providing:
1. Cromer-Mann coefficients (a1-a4, b1-b4, c) for X-ray form factors
2. Coherent scattering length for neutron weighting

The Cromer-Mann coefficients are tabulated in the International Tables for Crystallography, Volume C, or in Waasmaier & Kirfel (1995), *Acta Cryst.* A51, 416-431.
