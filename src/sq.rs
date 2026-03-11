use std::collections::HashMap;
use std::f64::consts::PI;

/// Compute partial S_ab(Q) from g_ab(r) via direct sine transform.
///
/// S(Q) = 1 + 4πρ₀ dr Σᵢ rᵢ [g(rᵢ)-1] W(rᵢ) sin(Q rᵢ) / Q
///
/// where W(r) is the optional Lorch window function.
pub fn compute_partial_sq(
    r: &[f64],
    gr: &[f64],
    rho0: f64,
    q: &[f64],
    lorch: bool,
) -> Vec<f64> {
    let nr = r.len();
    let dr = if nr > 1 { r[1] - r[0] } else { 1.0 };
    let r_max = r[nr - 1];
    let prefactor = 4.0 * PI * rho0 * dr;

    // Precompute f(r) = r * [g(r) - 1] * W(r)
    let f: Vec<f64> = (0..nr)
        .map(|i| {
            let mut h = gr[i] - 1.0;
            if lorch && r[i] > 0.0 {
                let arg = PI * r[i] / r_max;
                h *= arg.sin() / arg;
            }
            r[i] * h
        })
        .collect();

    q.iter()
        .map(|&qi| {
            if qi < 0.05 {
                return 1.0;
            }
            let sum: f64 = f.iter().zip(r.iter()).map(|(&fi, &ri)| fi * (qi * ri).sin()).sum();
            1.0 + prefactor * sum / qi
        })
        .collect()
}

/// Compute all partial S_ab(Q) from partial g_ab(r).
pub fn compute_all_partial_sq(
    r: &[f64],
    partials_gr: &HashMap<usize, Vec<f64>>,
    rho0: f64,
    q: &[f64],
    lorch: bool,
) -> HashMap<usize, Vec<f64>> {
    partials_gr
        .iter()
        .map(|(&pair_idx, gr)| (pair_idx, compute_partial_sq(r, gr, rho0, q, lorch)))
        .collect()
}
