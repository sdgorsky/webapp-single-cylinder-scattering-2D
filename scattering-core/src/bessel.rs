//! Complex Bessel functions for electromagnetic scattering calculations.
//!
//! Pure Rust implementation using series expansions and asymptotic formulas.

use num_complex::Complex64;
use std::f64::consts::PI;

/// Euler-Mascheroni constant
const EULER_GAMMA: f64 = 0.5772156649015329;

/// Factorial function for small integers
fn factorial(n: u32) -> u64 {
    match n {
        0 | 1 => 1,
        _ => (2..=n as u64).product(),
    }
}

/// Digamma function ψ(n) for positive integers
/// ψ(n) = -γ + Σ_{k=1}^{n-1} 1/k
fn digamma(n: u32) -> f64 {
    if n == 0 {
        return f64::NEG_INFINITY;
    }
    let mut result = -EULER_GAMMA;
    for k in 1..n {
        result += 1.0 / k as f64;
    }
    result
}

/// Bessel function of the first kind J_n(z) for complex argument.
/// Uses series expansion.
pub fn bessel_j(n: i32, z: Complex64) -> Complex64 {
    // Handle negative orders: J_{-n}(z) = (-1)^n * J_n(z)
    if n < 0 {
        let sign = if (-n) % 2 == 0 { 1.0 } else { -1.0 };
        return sign * bessel_j(-n, z);
    }

    // For very small z, use leading term
    if z.norm() < 1e-15 {
        return if n == 0 {
            Complex64::new(1.0, 0.0)
        } else {
            Complex64::new(0.0, 0.0)
        };
    }

    // Use series expansion: J_n(z) = (z/2)^n * sum_{k=0}^inf (-z^2/4)^k / (k! * (n+k)!)
    let z_half = z / 2.0;
    let z_half_sq = z_half * z_half;

    // (z/2)^n
    let mut prefix = Complex64::new(1.0, 0.0);
    for _ in 0..n {
        prefix *= z_half;
    }

    let mut sum = Complex64::new(0.0, 0.0);
    let mut factorial_k: f64 = 1.0;
    let mut factorial_n_plus_k: f64 = factorial(n as u32) as f64;

    // k=0 term
    sum += Complex64::new(1.0 / factorial_n_plus_k, 0.0);

    // Sum terms until convergence
    for k in 1..150 {
        factorial_k *= k as f64;
        factorial_n_plus_k *= n as f64 + k as f64;
        let term = (-z_half_sq).powi(k) / (factorial_k * factorial_n_plus_k);
        sum += term;

        if term.norm() < 1e-15 * sum.norm() {
            break;
        }
    }

    prefix * sum
}

/// Bessel function of the second kind Y_n(z) for complex argument.
pub fn bessel_y(n: i32, z: Complex64) -> Complex64 {
    // Handle negative orders: Y_{-n}(z) = (-1)^n * Y_n(z)
    if n < 0 {
        let sign = if (-n) % 2 == 0 { 1.0 } else { -1.0 };
        return sign * bessel_y(-n, z);
    }

    // Y_n(z) = (J_n(z) * cos(n*pi) - J_{-n}(z)) / sin(n*pi)
    // For integer n, use limit form:
    // Y_n(z) = (2/pi) * J_n(z) * ln(z/2) - (1/pi) * sum of correction terms

    let j_n = bessel_j(n, z);
    let ln_z_half = (z / 2.0).ln();

    // First part: (2/π) * J_n(z) * ln(z/2)
    let part1 = (2.0 / PI) * j_n * ln_z_half;

    // Second part: negative power series
    let z_half = z / 2.0;
    let z_half_sq = z_half * z_half;

    // (z/2)^n term
    let mut prefix = Complex64::new(1.0, 0.0);
    for _ in 0..n {
        prefix *= z_half;
    }

    // Sum for ψ(k+1) + ψ(n+k+1) terms
    let mut sum1 = Complex64::new(0.0, 0.0);
    let mut factorial_k: f64 = 1.0;
    let mut factorial_n_plus_k: f64 = factorial(n as u32) as f64;

    for k in 0..100 {
        if k > 0 {
            factorial_k *= k as f64;
            factorial_n_plus_k *= (n + k) as f64;
        }

        let psi_sum = digamma(k as u32 + 1) + digamma((n + k) as u32 + 1);
        let term = (-z_half_sq).powi(k) * psi_sum / (factorial_k * factorial_n_plus_k);
        sum1 += term;

        if k > 0 && term.norm() < 1e-15 * sum1.norm() {
            break;
        }
    }

    let part2 = -(1.0 / PI) * prefix * sum1;

    // Third part: (z/2)^{-n} * sum for k < n
    let mut part3 = Complex64::new(0.0, 0.0);
    if n > 0 {
        let mut z_half_neg_n = Complex64::new(1.0, 0.0);
        for _ in 0..n {
            z_half_neg_n /= z_half;
        }

        let mut sum2 = Complex64::new(0.0, 0.0);
        for k in 0..n {
            let factorial_n_minus_k_minus_1 = factorial((n - k - 1) as u32) as f64;
            let factorial_k_val = factorial(k as u32) as f64;
            let term = factorial_n_minus_k_minus_1 * z_half_sq.powi(k) / factorial_k_val;
            sum2 += term;
        }
        part3 = (1.0 / PI) * z_half_neg_n * sum2;
    }

    part1 + part2 - part3
}

/// Hankel function of the first kind H^(1)_n(z) = J_n(z) + i*Y_n(z)
pub fn hankel1(n: i32, z: Complex64) -> Complex64 {
    bessel_j(n, z) + Complex64::i() * bessel_y(n, z)
}

/// Derivative of Bessel function J'_n(z) using recurrence relation:
/// J'_n(z) = (J_{n-1}(z) - J_{n+1}(z)) / 2
pub fn bessel_j_derivative(n: i32, z: Complex64) -> Complex64 {
    (bessel_j(n - 1, z) - bessel_j(n + 1, z)) / 2.0
}

/// Derivative of Hankel function H^(1)'_n(z) using recurrence relation:
/// H'_n(z) = (H_{n-1}(z) - H_{n+1}(z)) / 2
pub fn hankel1_derivative(n: i32, z: Complex64) -> Complex64 {
    (hankel1(n - 1, z) - hankel1(n + 1, z)) / 2.0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: Complex64, b: Complex64, tol: f64) -> bool {
        (a - b).norm() < tol
    }

    #[test]
    fn test_j0_real() {
        // J_0(0) ≈ 1
        let j0_0 = bessel_j(0, Complex64::new(1e-15, 0.0));
        assert!(approx_eq(j0_0, Complex64::new(1.0, 0.0), 1e-6));

        // J_0(1) ≈ 0.7651976866
        let j0_1 = bessel_j(0, Complex64::new(1.0, 0.0));
        assert!(approx_eq(j0_1, Complex64::new(0.7651976866, 0.0), 1e-6));

        // J_0(5) ≈ -0.1775967713
        let j0_5 = bessel_j(0, Complex64::new(5.0, 0.0));
        assert!(approx_eq(j0_5, Complex64::new(-0.1775967713, 0.0), 1e-6));
    }

    #[test]
    fn test_j1_real() {
        // J_1(1) ≈ 0.4400505857
        let j1_1 = bessel_j(1, Complex64::new(1.0, 0.0));
        assert!(approx_eq(j1_1, Complex64::new(0.4400505857, 0.0), 1e-6));
    }

    #[test]
    fn test_negative_order() {
        // J_{-n}(z) = (-1)^n J_n(z) for integer n
        let z = Complex64::new(2.5, 0.0);

        // Even order: J_{-2} = J_2
        let j2 = bessel_j(2, z);
        let j_neg2 = bessel_j(-2, z);
        assert!(
            approx_eq(j2, j_neg2, 1e-10),
            "J_2 != J_-2: {} vs {}",
            j2,
            j_neg2
        );

        // Odd order: J_{-3} = -J_3
        let j3 = bessel_j(3, z);
        let j_neg3 = bessel_j(-3, z);
        assert!(
            approx_eq(j3, -j_neg3, 1e-10),
            "J_3 != -J_-3: {} vs {}",
            j3,
            j_neg3
        );
    }
}
