//! Central Composite Design (CCD).
//!
//! Combines a 2^k factorial core with 2k axial (star) points and center
//! runs to fit a second-order response surface model.
//!
//! # Structure
//!
//! 1. **Factorial portion** — 2^k full factorial at coded ±1 (or ±1/α for Inscribed)
//! 2. **Axial points** — ±α on each axis (2k points)
//! 3. **Center points** — `n_center` replicates at the origin
//!
//! # Axial Distance α
//!
//! | Type | α | Property |
//! |------|---|----------|
//! | Face-centered | 1.0 | Axial points on cube faces; 3-level design |
//! | Rotatable | (2^k)^(1/4) | Prediction variance uniform on spheres |
//! | Inscribed | (2^k)^(1/4) | Axial=±1, factorial scaled to ±(1/α) |
//!
//! Reference: Box, G.E.P. & Wilson, K.B. (1951). "On the Experimental
//! Attainment of Optimum Conditions". *Journal of the Royal Statistical
//! Society B*, 13(1), pp. 1–45.

use super::{factorial::full_factorial, DesignMatrix};
use crate::error::DoeError;

/// Method for choosing the axial distance α.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AlphaType {
    /// α = 1.0: axial points lie on the faces of the hypercube.
    FaceCentered,
    /// α = (2^k)^(1/4): prediction variance is constant on spheres.
    Rotatable,
    /// Axial points fixed at ±1; factorial points scaled to ±(1/α).
    Inscribed,
}

/// Compute the rotatable axial distance: (2^k)^(1/4).
///
/// # Examples
///
/// ```
/// use u_doe::design::ccd::rotatable_alpha;
/// let alpha = rotatable_alpha(2);
/// assert!((alpha - 2.0_f64.sqrt()).abs() < 1e-6);
/// ```
pub fn rotatable_alpha(k: usize) -> f64 {
    let factorial_runs = (1_usize << k) as f64;
    factorial_runs.powf(0.25)
}

/// Generate a Central Composite Design.
///
/// # Arguments
/// * `k`        — Number of factors (2..=6)
/// * `alpha`    — Method for choosing the axial distance
/// * `n_center` — Number of center point replicates (≥ 1)
///
/// # Errors
///
/// Returns `Err` if `k < 2`, `k > 6`, or `n_center == 0`.
///
/// # Examples
///
/// ```
/// use u_doe::design::ccd::{ccd, AlphaType};
/// let d = ccd(2, AlphaType::FaceCentered, 3).unwrap();
/// assert_eq!(d.run_count(), 11); // 4 + 4 + 3
/// assert_eq!(d.factor_count(), 2);
/// ```
pub fn ccd(k: usize, alpha_type: AlphaType, n_center: usize) -> Result<DesignMatrix, DoeError> {
    if !(2..=6).contains(&k) {
        return Err(DoeError::InvalidFactorCount { min: 2, max: 6, got: k });
    }
    if n_center == 0 {
        return Err(DoeError::InvalidSpecification(
            "n_center must be at least 1".to_string(),
        ));
    }

    let alpha_val = match alpha_type {
        AlphaType::FaceCentered => 1.0,
        AlphaType::Rotatable | AlphaType::Inscribed => rotatable_alpha(k),
    };

    let factorial_scale = if alpha_type == AlphaType::Inscribed {
        1.0 / alpha_val
    } else {
        1.0
    };

    let axial_val = if alpha_type == AlphaType::Inscribed {
        1.0
    } else {
        alpha_val
    };

    // 1. Factorial portion (scaled)
    let base = full_factorial(k)?;
    let mut data: Vec<Vec<f64>> = base
        .data
        .iter()
        .map(|row| row.iter().map(|&x| x * factorial_scale).collect())
        .collect();

    // 2. Axial points: ±axial_val on each axis, all others 0
    for j in 0..k {
        let mut lo = vec![0.0_f64; k];
        let mut hi = vec![0.0_f64; k];
        lo[j] = -axial_val;
        hi[j] = axial_val;
        data.push(lo);
        data.push(hi);
    }

    // 3. Center points
    for _ in 0..n_center {
        data.push(vec![0.0_f64; k]);
    }

    Ok(DesignMatrix {
        data,
        factor_names: DesignMatrix::default_names(k),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ccd_face_centered_k2_run_count() {
        // 2^2 + 2*2 + 3 = 4 + 4 + 3 = 11
        let d = ccd(2, AlphaType::FaceCentered, 3).unwrap();
        assert_eq!(d.run_count(), 11);
        assert_eq!(d.factor_count(), 2);
    }

    #[test]
    fn ccd_face_centered_k3_run_count() {
        // 2^3 + 2*3 + 1 = 8 + 6 + 1 = 15
        let d = ccd(3, AlphaType::FaceCentered, 1).unwrap();
        assert_eq!(d.run_count(), 15);
    }

    #[test]
    fn ccd_rotatable_alpha_k2() {
        // α = (2^2)^(1/4) = 4^0.25 = sqrt(2) ≈ 1.4142
        let alpha = rotatable_alpha(2);
        assert!((alpha - 2.0_f64.sqrt()).abs() < 1e-6);
    }

    #[test]
    fn ccd_rotatable_alpha_k3() {
        // α = (2^3)^(1/4) = 8^0.25 ≈ 1.6818
        let alpha = rotatable_alpha(3);
        assert!((alpha - 1.6818).abs() < 1e-3);
    }

    #[test]
    fn ccd_factorial_points_at_plus_minus_one() {
        let d = ccd(2, AlphaType::FaceCentered, 1).unwrap();
        // First 4 rows = factorial: all values ±1
        for run in 0..4 {
            for &v in &d.data[run] {
                assert!((v.abs() - 1.0).abs() < 1e-10,
                    "factorial point run {run} has value {v}");
            }
        }
    }

    #[test]
    fn ccd_axial_points_one_nonzero_each() {
        let d = ccd(3, AlphaType::FaceCentered, 1).unwrap();
        let factorial_runs = 1 << 3; // 8
        for i in 0..(2 * 3) {
            let run = factorial_runs + i;
            let nonzero: Vec<f64> = d.data[run].iter().copied()
                .filter(|&v| v.abs() > 1e-9)
                .collect();
            assert_eq!(nonzero.len(), 1,
                "axial run {i} should have exactly 1 nonzero factor");
        }
    }

    #[test]
    fn ccd_center_points_all_zero() {
        let d = ccd(2, AlphaType::FaceCentered, 3).unwrap();
        let center_start = 4 + 4; // factorial + axial
        for i in 0..3 {
            for &v in &d.data[center_start + i] {
                assert!(v.abs() < 1e-10, "center point not zero: {v}");
            }
        }
    }

    #[test]
    fn ccd_k1_error() {
        assert!(ccd(1, AlphaType::FaceCentered, 1).is_err());
    }

    #[test]
    fn ccd_k7_error() {
        assert!(ccd(7, AlphaType::FaceCentered, 1).is_err());
    }

    #[test]
    fn ccd_n_center_zero_error() {
        assert!(ccd(2, AlphaType::FaceCentered, 0).is_err());
    }

    #[test]
    fn ccd_inscribed_axial_at_one() {
        let d = ccd(2, AlphaType::Inscribed, 1).unwrap();
        // Axial points for Inscribed are ±1
        let factorial_runs = 4;
        for i in 0..(2 * 2) {
            let run = factorial_runs + i;
            let nonzero: Vec<f64> = d.data[run].iter().copied()
                .filter(|&v| v.abs() > 1e-9)
                .collect();
            assert_eq!(nonzero.len(), 1);
            assert!((nonzero[0].abs() - 1.0).abs() < 1e-10,
                "inscribed axial should be ±1");
        }
    }
}
