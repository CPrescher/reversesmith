# Supported Elements

All elements from hydrogen (Z=1) through californium (Z=98) are supported for X-ray form factors, molar masses, and (where applicable) neutron scattering.

## X-ray form factors

Cromer-Mann 9-parameter fits from the International Tables for Crystallography, Vol C (2004). All neutral atoms Z=1–98 are included:

H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf

## Neutron scattering lengths

Coherent neutron scattering lengths (fm) from the NIST Center for Neutron Research are provided for most stable elements up to uranium. Notable elements with large absorption cross-sections (e.g., Cd, Sm, Gd) are included but should be used with care.

## Molar masses

IUPAC 2021 standard atomic weights for Z=1–98. Used for mass density computation and box rescaling.

## Adding ionic form factors

The current tables use neutral-atom form factors. If ionic form factors are needed (e.g., Ca²⁺, O¹⁻), additional entries can be added to the `cromer_mann_params` match table in `src/xray.rs`. The International Tables tabulate coefficients for common valence states.
