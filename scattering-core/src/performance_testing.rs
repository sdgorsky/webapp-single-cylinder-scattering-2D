use crate::bessel::hankel1;
use num_complex::Complex64;
use std::f64::consts::PI;

const GRID_SIZE: usize = 512;
const HALF_GRID_SIZE: usize = GRID_SIZE / 2;
const VIEW_SIZE: f64 = 5.0;
const TOTAL_POINTS: usize = GRID_SIZE * GRID_SIZE;
const RADIUS: f64 = 0.5;

const DS: f64 = VIEW_SIZE / (GRID_SIZE as f64);
const K0: f64 = 2.0 * PI; // Wavelength = 1.0

// Test scattering coefficients (21 orders: -10 to +10)
const MAX_ORDER: i32 = 10;
const NUM_ORDERS: usize = (2 * MAX_ORDER + 1) as usize;

pub fn create_coord_vec() -> (Vec<f64>, Vec<f64>) {
    let mut x_pos = vec![0.0; TOTAL_POINTS];
    let mut y_pos = vec![0.0; TOTAL_POINTS];
    for iy in 0..GRID_SIZE / 2 {
        for ix in 0..GRID_SIZE {
            let idx = iy * GRID_SIZE + ix;
            x_pos[idx] = (ix as f64 + 0.5) * DS;
            y_pos[idx] = (iy as f64 + 0.5) * DS;
        }
    }
    (x_pos, y_pos)
}

pub fn create_coord_array_boxed() -> (Box<[f64; TOTAL_POINTS]>, Box<[f64; TOTAL_POINTS]>) {
    let mut x_pos = Box::new([0.0; TOTAL_POINTS]);
    let mut y_pos = Box::new([0.0; TOTAL_POINTS]);
    for iy in 0..GRID_SIZE / 2 {
        for ix in 0..GRID_SIZE {
            let idx = iy * GRID_SIZE + ix;
            x_pos[idx] = (ix as f64 + 0.5) * DS;
            y_pos[idx] = (iy as f64 + 0.5) * DS;
        }
    }
    (x_pos, y_pos)
}

pub fn create_coord_array() -> ([f64; TOTAL_POINTS], [f64; TOTAL_POINTS]) {
    let mut x_pos: [f64; TOTAL_POINTS] = [0.0; TOTAL_POINTS];
    let mut y_pos: [f64; TOTAL_POINTS] = [0.0; TOTAL_POINTS];
    for iy in 0..GRID_SIZE / 2 {
        for ix in 0..GRID_SIZE {
            let idx = iy * GRID_SIZE + ix;
            x_pos[idx] = (ix as f64 + 0.5) * DS;
            y_pos[idx] = (iy as f64 + 0.5) * DS;
        }
    }
    (x_pos, y_pos)
}

pub fn pw_vec_1() -> (Vec<f64>, Vec<f64>) {
    let mut field_real = vec![0.0; TOTAL_POINTS];
    let mut field_imag = vec![0.0; TOTAL_POINTS];

    for iy in 0..GRID_SIZE {
        for ix in 0..GRID_SIZE {
            let idx: usize = iy * GRID_SIZE + ix;
            let x = (ix as f64 + 0.5) * DS;

            field_real[idx] = (K0 * x).cos();
            field_imag[idx] = (K0 * x).sin();
        }
    }
    (field_real, field_imag)
}

pub fn pw_vec_2() -> (Vec<f64>, Vec<f64>) {
    let mut field_real = vec![0.0; TOTAL_POINTS];
    let mut field_imag = vec![0.0; TOTAL_POINTS];

    // Grid setup
    let half_size = VIEW_SIZE / 2.0;
    let x_min = -half_size;

    for ix in 0..GRID_SIZE {
        let x = x_min + (ix as f64 + 0.5) * DS;
        let f_real = (K0 * x).cos();
        let f_imag = (K0 * x).sin();

        for iy in 0..GRID_SIZE {
            let idx = iy * GRID_SIZE + ix;
            field_real[idx] = f_real;
            field_imag[idx] = f_imag;
        }
    }
    (field_real, field_imag)
}

pub fn pw_vec_3() -> (Vec<f64>, Vec<f64>) {
    let mut field_real = vec![0.0; TOTAL_POINTS];
    let mut field_imag = vec![0.0; TOTAL_POINTS];

    // Grid setup
    let half_size = VIEW_SIZE / 2.0;
    let x_min = -half_size;

    // Precompute trig values for all columns (using symmetry: only compute half)
    let mut col_real = vec![0.0; GRID_SIZE];
    let mut col_imag = vec![0.0; GRID_SIZE];

    for ix in 0..GRID_SIZE / 2 {
        let kox = K0 * (x_min + (ix as f64 + 0.5) * DS);
        let (s, c) = kox.sin_cos();
        col_real[ix] = c;
        col_imag[ix] = s;
        // Symmetric column
        let ix_sym = GRID_SIZE - 1 - ix;
        col_real[ix_sym] = c;
        col_imag[ix_sym] = -s;
    }

    // Single sequential pass through memory
    for iy in 0..GRID_SIZE {
        let row_offset = iy * GRID_SIZE;
        for ix in 0..GRID_SIZE {
            let idx = row_offset + ix;
            field_real[idx] = col_real[ix];
            field_imag[idx] = col_imag[ix];
        }
    }
    (field_real, field_imag)
}

/// Generate test scattering coefficients (simple exponential decay)
fn get_test_coefficients() -> Vec<Complex64> {
    let mut coeffs = Vec::with_capacity(NUM_ORDERS);
    for n in -MAX_ORDER..=MAX_ORDER {
        // Simple test coefficients that decay with order
        let mag = 0.5 * (-0.1 * (n.abs() as f64)).exp();
        coeffs.push(Complex64::new(mag, mag * 0.1));
    }
    coeffs
}

/// Original (naive) exterior field computation - for comparison
pub fn exterior_field_naive() -> (Vec<f64>, Vec<f64>) {
    let mut field_real = vec![0.0; TOTAL_POINTS];
    let mut field_imag = vec![0.0; TOTAL_POINTS];

    let half_size = VIEW_SIZE / 2.0;
    let b_n = get_test_coefficients();

    for iy in 0..GRID_SIZE {
        let y = half_size - (iy as f64 + 0.5) * DS;

        for ix in 0..GRID_SIZE {
            let x = -half_size + (ix as f64 + 0.5) * DS;
            let r = (x * x + y * y).sqrt();
            let theta = y.atan2(x);

            let idx = iy * GRID_SIZE + ix;

            if r < RADIUS {
                // Skip interior points for this benchmark
                field_real[idx] = 0.0;
                field_imag[idx] = 0.0;
                continue;
            }

            let k0_r = Complex64::new(K0 * r, 0.0);
            let k0x = K0 * r * theta.cos();
            let mut field = Complex64::new(0.0, 0.0);

            for (i, n) in (-MAX_ORDER..=MAX_ORDER).enumerate() {
                let hn = hankel1(n, k0_r);
                let n_theta = n as f64 * theta;
                let exp_intheta = Complex64::new(n_theta.cos(), n_theta.sin());
                let phiz_incident = Complex64::new(k0x.cos(), k0x.sin());
                field += phiz_incident + b_n[i] * hn * exp_intheta;
            }

            field_real[idx] = field.re;
            field_imag[idx] = field.im;
        }
    }

    (field_real, field_imag)
}

/// Optimized exterior field: precompute Hankel values, then sequential fill
pub fn exterior_field_optimized() -> (Vec<f64>, Vec<f64>) {
    let mut field_real = vec![0.0; TOTAL_POINTS];
    let mut field_imag = vec![0.0; TOTAL_POINTS];

    let half_size = VIEW_SIZE / 2.0;
    let b_n = get_test_coefficients();

    // ===== PHASE 1: Precompute Hankel values for first quadrant =====
    // Using 4-fold symmetry: (x,y), (-x,y), (x,-y), (-x,-y) have same radius
    // Store Hankel values indexed by (qx, qy) where qx, qy are quadrant indices

    let quad_size = HALF_GRID_SIZE;
    let quad_points = quad_size * quad_size;

    // For each quadrant point, store Hankel values for all orders
    // hankel_cache[qy * quad_size + qx][order_idx] = H_n(k0*r)
    let mut hankel_cache: Vec<Vec<Complex64>> = Vec::with_capacity(quad_points);
    let mut radius_cache: Vec<f64> = Vec::with_capacity(quad_points);
    let mut is_exterior: Vec<bool> = Vec::with_capacity(quad_points);

    for qy in 0..quad_size {
        // y coordinate for this quadrant row (positive y)
        let y = (qy as f64 + 0.5) * DS;

        for qx in 0..quad_size {
            // x coordinate for this quadrant column (positive x)
            let x = (qx as f64 + 0.5) * DS;
            let r = (x * x + y * y).sqrt();

            radius_cache.push(r);

            if r < RADIUS {
                is_exterior.push(false);
                hankel_cache.push(Vec::new()); // Empty for interior points
            } else {
                is_exterior.push(true);
                let k0_r = Complex64::new(K0 * r, 0.0);
                let mut hn_values = Vec::with_capacity(NUM_ORDERS);
                for n in -MAX_ORDER..=MAX_ORDER {
                    hn_values.push(hankel1(n, k0_r));
                }
                hankel_cache.push(hn_values);
            }
        }
    }

    // ===== PHASE 2: Sequential fill using precomputed values =====

    for iy in 0..GRID_SIZE {
        let y = half_size - (iy as f64 + 0.5) * DS;
        // Map to quadrant y index
        let qy = if y >= 0.0 {
            ((y / DS - 0.5).round() as usize).min(quad_size - 1)
        } else {
            ((-y / DS - 0.5).round() as usize).min(quad_size - 1)
        };

        let row_offset = iy * GRID_SIZE;

        for ix in 0..GRID_SIZE {
            let x = -half_size + (ix as f64 + 0.5) * DS;
            // Map to quadrant x index
            let qx = if x >= 0.0 {
                ((x / DS - 0.5).round() as usize).min(quad_size - 1)
            } else {
                ((-x / DS - 0.5).round() as usize).min(quad_size - 1)
            };

            let q_idx = qy * quad_size + qx;
            let idx = row_offset + ix;

            if !is_exterior[q_idx] {
                field_real[idx] = 0.0;
                field_imag[idx] = 0.0;
                continue;
            }

            let r = radius_cache[q_idx];
            let theta = y.atan2(x);
            let k0x = K0 * r * theta.cos();

            let hn_values = &hankel_cache[q_idx];
            let mut field = Complex64::new(0.0, 0.0);

            for (i, n) in (-MAX_ORDER..=MAX_ORDER).enumerate() {
                let n_theta = n as f64 * theta;
                let exp_intheta = Complex64::new(n_theta.cos(), n_theta.sin());
                let phiz_incident = Complex64::new(k0x.cos(), k0x.sin());
                field += phiz_incident + b_n[i] * hn_values[i] * exp_intheta;
            }

            field_real[idx] = field.re;
            field_imag[idx] = field.im;
        }
    }

    (field_real, field_imag)
}

// ===== Spline-based exterior field computation =====

/// Number of points in the radial spline grid
const SPLINE_GRID_SIZE: usize = 256;

/// Maximum radius on the grid (corner of the view)
const R_MAX: f64 = 3.536; // ~sqrt(2) * VIEW_SIZE / 2 = sqrt(2) * 2.5

/// Cubic spline interpolation coefficients for a single Hankel order
/// Stores coefficients for both real and imaginary parts
struct HankelSpline {
    r_min: f64,
    dr: f64,
    // Spline coefficients: for each interval [i, i+1], store a, b, c, d
    // where f(x) = a + b*(x-x_i) + c*(x-x_i)^2 + d*(x-x_i)^3
    coeffs_re: Vec<[f64; 4]>, // Real part coefficients
    coeffs_im: Vec<[f64; 4]>, // Imaginary part coefficients
}

impl HankelSpline {
    /// Create a cubic spline for H_n(k0*r) over the range [r_min, r_max]
    fn new(order: i32, k0: f64, r_min: f64, r_max: f64, n_points: usize) -> Self {
        let dr = (r_max - r_min) / (n_points - 1) as f64;

        // Evaluate Hankel function at grid points
        let mut values_re = Vec::with_capacity(n_points);
        let mut values_im = Vec::with_capacity(n_points);

        for i in 0..n_points {
            let r = r_min + i as f64 * dr;
            let k0_r = Complex64::new(k0 * r, 0.0);
            let hn = hankel1(order, k0_r);
            values_re.push(hn.re);
            values_im.push(hn.im);
        }

        // Compute cubic spline coefficients using natural spline boundary conditions
        let coeffs_re = Self::compute_spline_coefficients(&values_re, dr);
        let coeffs_im = Self::compute_spline_coefficients(&values_im, dr);

        HankelSpline {
            r_min,
            dr,
            coeffs_re,
            coeffs_im,
        }
    }

    /// Compute natural cubic spline coefficients for uniform grid
    fn compute_spline_coefficients(y: &[f64], h: f64) -> Vec<[f64; 4]> {
        let n = y.len();
        if n < 2 {
            return vec![];
        }

        // Solve tridiagonal system for second derivatives (natural spline: M_0 = M_{n-1} = 0)
        let mut m = vec![0.0; n]; // Second derivatives at each point

        if n > 2 {
            // Set up tridiagonal system: [1, 4, 1] * M = 6/h^2 * (y_{i+1} - 2*y_i + y_{i-1})
            let mut rhs = vec![0.0; n - 2];
            for i in 1..n - 1 {
                rhs[i - 1] = 6.0 / (h * h) * (y[i + 1] - 2.0 * y[i] + y[i - 1]);
            }

            // Thomas algorithm for tridiagonal [1, 4, 1] system
            let mut c_prime = vec![0.0; n - 2];
            let mut d_prime = vec![0.0; n - 2];

            // Forward sweep
            c_prime[0] = 1.0 / 4.0;
            d_prime[0] = rhs[0] / 4.0;

            for i in 1..n - 2 {
                let denom = 4.0 - c_prime[i - 1];
                c_prime[i] = 1.0 / denom;
                d_prime[i] = (rhs[i] - d_prime[i - 1]) / denom;
            }

            // Back substitution
            m[n - 2] = d_prime[n - 3];
            for i in (1..n - 2).rev() {
                m[i] = d_prime[i - 1] - c_prime[i - 1] * m[i + 1];
            }
        }

        // Compute coefficients for each interval
        let mut coeffs = Vec::with_capacity(n - 1);
        for i in 0..n - 1 {
            let a = y[i];
            let b = (y[i + 1] - y[i]) / h - h * (2.0 * m[i] + m[i + 1]) / 6.0;
            let c = m[i] / 2.0;
            let d = (m[i + 1] - m[i]) / (6.0 * h);
            coeffs.push([a, b, c, d]);
        }

        coeffs
    }

    /// Evaluate the spline at radius r, returning (real, imag) parts
    #[inline]
    fn eval(&self, r: f64) -> (f64, f64) {
        // Find interval index
        let t = (r - self.r_min) / self.dr;
        let i = (t as usize).min(self.coeffs_re.len() - 1);
        let x = r - (self.r_min + i as f64 * self.dr);

        // Evaluate cubic polynomial: a + b*x + c*x^2 + d*x^3
        let [a_re, b_re, c_re, d_re] = self.coeffs_re[i];
        let [a_im, b_im, c_im, d_im] = self.coeffs_im[i];

        let re = a_re + x * (b_re + x * (c_re + x * d_re));
        let im = a_im + x * (b_im + x * (c_im + x * d_im));

        (re, im)
    }
}

/// Optimized exterior field using spline interpolation for Hankel functions
pub fn exterior_field_spline() -> (Vec<f64>, Vec<f64>) {
    let mut field_real = vec![0.0; TOTAL_POINTS];
    let mut field_imag = vec![0.0; TOTAL_POINTS];

    let half_size = VIEW_SIZE / 2.0;
    let b_n = get_test_coefficients();

    // ===== PHASE 1: Build splines for each Hankel order =====
    let mut hankel_splines: Vec<HankelSpline> = Vec::with_capacity(NUM_ORDERS);
    for n in -MAX_ORDER..=MAX_ORDER {
        hankel_splines.push(HankelSpline::new(n, K0, RADIUS, R_MAX, SPLINE_GRID_SIZE));
    }

    // ===== PHASE 2: Compute field using spline interpolation =====
    for iy in 0..GRID_SIZE {
        let y = half_size - (iy as f64 + 0.5) * DS;
        let row_offset = iy * GRID_SIZE;

        for ix in 0..GRID_SIZE {
            let x = -half_size + (ix as f64 + 0.5) * DS;
            let r = (x * x + y * y).sqrt();
            let idx = row_offset + ix;

            if r < RADIUS {
                field_real[idx] = 0.0;
                field_imag[idx] = 0.0;
                continue;
            }

            let theta = y.atan2(x);
            let k0x = K0 * r * theta.cos();
            let mut field = Complex64::new(0.0, 0.0);

            for (i, n) in (-MAX_ORDER..=MAX_ORDER).enumerate() {
                // Interpolate Hankel function
                let (hn_re, hn_im) = hankel_splines[i].eval(r);
                let hn = Complex64::new(hn_re, hn_im);

                let n_theta = n as f64 * theta;
                let exp_intheta = Complex64::new(n_theta.cos(), n_theta.sin());
                let phiz_incident = Complex64::new(k0x.cos(), k0x.sin());
                field += phiz_incident + b_n[i] * hn * exp_intheta;
            }

            field_real[idx] = field.re;
            field_imag[idx] = field.im;
        }
    }

    (field_real, field_imag)
}

/// Optimized exterior field using spline + y-axis symmetry
/// Only computes top half (y >= 0) and mirrors to bottom half
pub fn exterior_field_spline_symmetric() -> (Vec<f64>, Vec<f64>) {
    let mut field_real = vec![0.0; TOTAL_POINTS];
    let mut field_imag = vec![0.0; TOTAL_POINTS];

    let half_size = VIEW_SIZE / 2.0;
    let b_n = get_test_coefficients();

    // ===== PHASE 1: Build splines for each Hankel order =====
    let mut hankel_splines: Vec<HankelSpline> = Vec::with_capacity(NUM_ORDERS);
    for n in -MAX_ORDER..=MAX_ORDER {
        hankel_splines.push(HankelSpline::new(n, K0, RADIUS, R_MAX, SPLINE_GRID_SIZE));
    }

    // ===== PHASE 2: Compute field for top half only (y >= 0) =====
    // iy=0 is top row (y = half_size - 0.5*DS), iy=GRID_SIZE/2-1 is near y=0
    for iy in 0..GRID_SIZE / 2 {
        let y = half_size - (iy as f64 + 0.5) * DS;
        let row_offset = iy * GRID_SIZE;

        // Symmetric row: iy_sym = GRID_SIZE - 1 - iy (mirrors about center)
        let iy_sym = GRID_SIZE - 1 - iy;
        let row_offset_sym = iy_sym * GRID_SIZE;

        for ix in 0..GRID_SIZE {
            let x = -half_size + (ix as f64 + 0.5) * DS;
            let r = (x * x + y * y).sqrt();
            let idx = row_offset + ix;
            let idx_sym = row_offset_sym + ix;

            if r < RADIUS {
                field_real[idx] = 0.0;
                field_imag[idx] = 0.0;
                field_real[idx_sym] = 0.0;
                field_imag[idx_sym] = 0.0;
                continue;
            }

            let theta = y.atan2(x);
            let k0x = K0 * r * theta.cos();
            let mut field = Complex64::new(0.0, 0.0);

            for (i, n) in (-MAX_ORDER..=MAX_ORDER).enumerate() {
                // Interpolate Hankel function
                let (hn_re, hn_im) = hankel_splines[i].eval(r);
                let hn = Complex64::new(hn_re, hn_im);

                let n_theta = n as f64 * theta;
                let exp_intheta = Complex64::new(n_theta.cos(), n_theta.sin());
                let phiz_incident = Complex64::new(k0x.cos(), k0x.sin());
                field += phiz_incident + b_n[i] * hn * exp_intheta;
            }

            // Fill both the computed point and its symmetric counterpart
            field_real[idx] = field.re;
            field_imag[idx] = field.im;
            field_real[idx_sym] = field.re; // E_z(x, -y) = E_z(x, y)
            field_imag[idx_sym] = field.im;
        }
    }

    (field_real, field_imag)
}
