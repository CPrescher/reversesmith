use rustfft::num_complex::Complex64;
use rustfft::{Fft, FftPlanner};
use std::sync::Arc;

/// Chirp-Z Transform for computing the discrete sine transform:
///   Y[k] = Σ_{i=0}^{N-1} f[i] * sin(Q_k * r_i)
/// where Q_k = qmin + k * dQ and r_i are uniformly spaced,
/// using Bluestein's algorithm (O(L log L) instead of O(N * M)).
///
/// The key idea: the sum is Im(Σ f[i] * exp(j * Q_k * r_i)), which is
/// a non-uniform DFT that can be expressed as a convolution via
/// Bluestein's identity: k*i = (k² + i² - (k-i)²)/2.
pub struct CztSineTransform {
    nbins: usize,
    nq: usize,
    fft_len: usize,
    /// Input chirp: exp(j * (qmin * r_i + α * i² / 2))
    pre_chirp: Vec<Complex64>,
    /// Output phase: (cos θ_k, sin θ_k) where θ_k = α * k * (k+1) / 2
    post_chirp: Vec<(f64, f64)>,
    /// FFT of the convolution kernel
    kernel_fft: Vec<Complex64>,
    fft: Arc<dyn Fft<f64>>,
    ifft: Arc<dyn Fft<f64>>,
    inv_fft_len: f64,
    buf: Vec<Complex64>,
    scratch: Vec<Complex64>,
}

impl CztSineTransform {
    /// Create a new CZT for the given Q-grid and r-grid.
    ///
    /// Both grids must be uniformly spaced.
    pub fn new(q_grid: &[f64], r_grid: &[f64]) -> Self {
        let nbins = r_grid.len();
        let nq = q_grid.len();
        let dq = if nq > 1 { q_grid[1] - q_grid[0] } else { 1.0 };
        let dr = if nbins > 1 {
            r_grid[1] - r_grid[0]
        } else {
            1.0
        };
        let qmin = q_grid[0];
        let alpha = dq * dr;

        let fft_len = (nbins + nq - 1).next_power_of_two();

        // Pre-chirp: exp(j * (qmin * r_i + α * i² / 2))
        let pre_chirp: Vec<Complex64> = (0..nbins)
            .map(|i| {
                let phase = qmin * r_grid[i] + alpha * (i * i) as f64 / 2.0;
                Complex64::new(phase.cos(), phase.sin())
            })
            .collect();

        // Post-chirp phase: θ_k = α * k * (k + 1) / 2
        // Y[k] = Im(exp(j*θ) * c[k]) = cos(θ)*c_im + sin(θ)*c_re
        let post_chirp: Vec<(f64, f64)> = (0..nq)
            .map(|k| {
                let phase = alpha * k as f64 * (k as f64 + 1.0) / 2.0;
                (phase.cos(), phase.sin())
            })
            .collect();

        // Convolution kernel h[n] = exp(-j * α * n² / 2), stored circularly
        let mut h = vec![Complex64::new(0.0, 0.0); fft_len];
        // Positive indices: n = 0..min(nq, fft_len)
        for n in 0..nq.min(fft_len) {
            let phase = -alpha * (n * n) as f64 / 2.0;
            h[n] = Complex64::new(phase.cos(), phase.sin());
        }
        // Negative indices wrapped: h[L-m] = h(-m) = h(m) since h depends on n²
        for m in 1..nbins.min(fft_len) {
            let phase = -alpha * (m * m) as f64 / 2.0;
            h[fft_len - m] = Complex64::new(phase.cos(), phase.sin());
        }

        let mut planner = FftPlanner::<f64>::new();
        let fft = planner.plan_fft_forward(fft_len);
        let ifft = planner.plan_fft_inverse(fft_len);

        let scratch_len = fft
            .get_inplace_scratch_len()
            .max(ifft.get_inplace_scratch_len());
        let mut scratch = vec![Complex64::new(0.0, 0.0); scratch_len];

        let mut kernel_fft = h;
        fft.process_with_scratch(&mut kernel_fft, &mut scratch);

        CztSineTransform {
            nbins,
            nq,
            fft_len,
            pre_chirp,
            post_chirp,
            kernel_fft,
            fft,
            ifft,
            inv_fft_len: 1.0 / fft_len as f64,
            buf: vec![Complex64::new(0.0, 0.0); fft_len],
            scratch,
        }
    }

    /// Compute Y[k] = Σ_{i=0}^{N-1} input[i] * sin(Q_k * r_i) for k = 0..nq-1.
    pub fn transform(&mut self, input: &[f64], output: &mut [f64]) {
        debug_assert_eq!(input.len(), self.nbins);
        debug_assert_eq!(output.len(), self.nq);

        // Pre-chirp: a[i] = input[i] * pre_chirp[i]
        for i in 0..self.nbins {
            let f = input[i];
            let c = self.pre_chirp[i];
            self.buf[i] = Complex64::new(f * c.re, f * c.im);
        }
        for i in self.nbins..self.fft_len {
            self.buf[i] = Complex64::new(0.0, 0.0);
        }

        // Forward FFT
        self.fft
            .process_with_scratch(&mut self.buf, &mut self.scratch);

        // Pointwise multiply with kernel
        for i in 0..self.fft_len {
            let a = self.buf[i];
            let b = self.kernel_fft[i];
            self.buf[i] = Complex64::new(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
        }

        // Inverse FFT
        self.ifft
            .process_with_scratch(&mut self.buf, &mut self.scratch);

        // Post-chirp, normalize, extract Im
        let s = self.inv_fft_len;
        for k in 0..self.nq {
            let c_re = self.buf[k].re * s;
            let c_im = self.buf[k].im * s;
            let (p_cos, p_sin) = self.post_chirp[k];
            output[k] = p_cos * c_im + p_sin * c_re;
        }
    }
}
