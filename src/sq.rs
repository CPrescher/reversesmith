use std::collections::HashMap;
use std::f64::consts::PI;

use rustfft::num_complex::Complex;
use rustfft::FftPlanner;

/// Next power of 2 >= n.
fn next_pow2(n: usize) -> usize {
    let mut p = 1;
    while p < n {
        p <<= 1;
    }
    p
}

/// Compute partial S_ab(Q) from g_ab(r) via DST-II (FFT-based) + interpolation.
///
/// The r-grid has bin centres at r_i = (i+0.5)*dr.
/// We compute: I(Q) = dr * Σ_i  r_i * h(r_i) * sin(Q * r_i)
/// as a DST-II via FFT, yielding S(Q) on a fine uniform Q-grid,
/// then interpolate to the requested Q points.
///
/// DST-II identity:
///   Σ_{n=0}^{N-1} x_n sin(π(k+1)(2n+1)/(2N)) = -Im(exp(-iθ_k) * Y_{k+1})
/// where Y = FFT(x zero-padded to 2N), θ_k = π(k+1)/(2N).
/// Output Q-grid: Q_k = π(k+1)/(N_fft * dr).
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

    // Precompute f(r) = r * [g(r) - 1] * W(r)
    let mut f = vec![0.0; nr];
    for i in 0..nr {
        let mut h = gr[i] - 1.0;
        if lorch && r[i] > 0.0 {
            let arg = PI * r[i] / r_max;
            h *= arg.sin() / arg;
        }
        f[i] = r[i] * h;
    }

    // Zero-pad to N_fft (at least 8x nr for good interpolation resolution)
    let n_fft = next_pow2(nr * 8).max(4096);
    let fft_len = 2 * n_fft;

    // Prepare FFT input: f zero-padded to length fft_len
    let mut data: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); fft_len];
    for i in 0..nr {
        data[i] = Complex::new(f[i], 0.0);
    }

    // FFT (forward: Y_m = Σ y_n exp(-i 2π m n / M))
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(fft_len);
    fft.process(&mut data);

    // Extract DST-II coefficients and build S(Q) on the FFT Q-grid
    // Q_k = π(k+1) / (n_fft * dr), k = 0..n_fft-1
    let dq_fft = PI / (n_fft as f64 * dr);
    let prefactor = 4.0 * PI * rho0 * dr;

    // Only compute up to Q slightly beyond what we need
    let q_max_needed = q.iter().cloned().fold(0.0f64, f64::max);
    let k_max = ((q_max_needed / dq_fft).ceil() as usize + 2).min(n_fft);

    let mut q_fft = Vec::with_capacity(k_max);
    let mut sq_fft = Vec::with_capacity(k_max);

    for k in 0..k_max {
        let qi = (k + 1) as f64 * dq_fft;
        q_fft.push(qi);

        // DST-II[k] = -Im(exp(-iθ) * Y[k+1])
        let theta = PI * (k + 1) as f64 / (2.0 * n_fft as f64);
        let twiddle = Complex::new(theta.cos(), -theta.sin());
        let y = data[k + 1];
        let dst_val = -(twiddle * y).im;

        let integral = dr * dst_val;
        if qi > 0.05 {
            sq_fft.push(1.0 + prefactor * integral / (qi * dr));
        } else {
            sq_fft.push(1.0);
        }
    }

    // Interpolate onto requested Q-grid (parallel)
    q.iter()
        .map(|&qi| {
            if qi < 0.05 {
                return 1.0;
            }
            interpolate_linear(&q_fft, &sq_fft, qi)
        })
        .collect()
}

/// Linear interpolation with binary search.
fn interpolate_linear(x: &[f64], y: &[f64], xi: f64) -> f64 {
    let n = x.len();
    if n == 0 {
        return 0.0;
    }
    if xi <= x[0] {
        return y[0];
    }
    if xi >= x[n - 1] {
        return y[n - 1];
    }

    let mut lo = 0;
    let mut hi = n - 1;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if x[mid] <= xi {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    let t = (xi - x[lo]) / (x[hi] - x[lo]);
    y[lo] + t * (y[hi] - y[lo])
}

/// Compute all partial S_ab(Q) from partial g_ab(r) — pairs computed in parallel.
pub fn compute_all_partial_sq(
    r: &[f64],
    partials_gr: &HashMap<usize, Vec<f64>>,
    rho0: f64,
    q: &[f64],
    lorch: bool,
) -> HashMap<usize, Vec<f64>> {
    let pairs: Vec<(usize, &Vec<f64>)> = partials_gr.iter().map(|(&k, v)| (k, v)).collect();

    let results: Vec<(usize, Vec<f64>)> = pairs
        .iter()
        .map(|&(pair_idx, gr)| (pair_idx, compute_partial_sq(r, gr, rho0, q, lorch)))
        .collect();

    results.into_iter().collect()
}
