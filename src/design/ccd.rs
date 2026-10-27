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

use super::{factorial::full_factorial, ClassifiedDesign, DesignMatrix, PointType};
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
/// * `n_center` — Number of center point replicates (1 to [`MAX_CENTER_POINTS`](crate::design::MAX_CENTER_POINTS))
///
/// # Errors
///
/// Returns `Err` if `k < 2`, `k > 6`, or `n_center` is not in `1..=MAX_CENTER_POINTS`.
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
    ccd_design(k, alpha_type, n_center).map(|d| d.matrix)
}

/// Generate a Central Composite Design with the block of each run: the
/// factorial rows first, then the axial pair of each factor in turn, then the
/// center rows -- the order of the matrix, which is how a run sheet or a plot
/// tells the three apart without guessing from the coordinates (an inscribed
/// design's factorial level is not an integer).
///
/// Arguments and errors as for [`ccd`].
///
/// # Examples
///
/// ```
/// use u_doe::design::ccd::{ccd_design, AlphaType};
/// use u_doe::design::PointType;
/// let d = ccd_design(2, AlphaType::FaceCentered, 3).unwrap();
/// assert_eq!(d.point_types.iter().filter(|&&t| t == PointType::Axial).count(), 4);
/// assert_eq!(d.point_types.len(), d.matrix.run_count());
/// ```
pub fn ccd_design(
    k: usize,
    alpha_type: AlphaType,
    n_center: usize,
) -> Result<ClassifiedDesign, DoeError> {
    if !(2..=6).contains(&k) {
        return Err(DoeError::count_out_of_range("k", 2..=6, k));
    }
    if !(1..=super::MAX_CENTER_POINTS).contains(&n_center) {
        return Err(DoeError::count_out_of_range(
            "n_center",
            1..=super::MAX_CENTER_POINTS,
            n_center,
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

    let mut point_types = vec![PointType::Factorial; data.len()];

    // 2. Axial points: ±axial_val on each axis, all others 0
    for j in 0..k {
        let mut lo = vec![0.0_f64; k];
        let mut hi = vec![0.0_f64; k];
        lo[j] = -axial_val;
        hi[j] = axial_val;
        data.push(lo);
        data.push(hi);
        point_types.extend([PointType::Axial; 2]);
    }

    // 3. Center points
    for _ in 0..n_center {
        data.push(vec![0.0_f64; k]);
        point_types.push(PointType::Center);
    }

    Ok(ClassifiedDesign {
        matrix: DesignMatrix {
            data,
            factor_names: DesignMatrix::default_names(k),
        },
        point_types,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `n_center` is a number of rows to build: unbounded, it exhausted memory
    /// before anything refused it.
    #[test]
    fn n_center_is_bounded() {
        use crate::design::MAX_CENTER_POINTS;
        let too_many = MAX_CENTER_POINTS + 1;
        assert_eq!(
            ccd(2, AlphaType::FaceCentered, too_many).unwrap_err(),
            DoeError::count_out_of_range("n_center", 1..=MAX_CENTER_POINTS, too_many)
        );
        assert_eq!(
            crate::design::box_behnken::box_behnken_design(3, usize::MAX).unwrap_err(),
            DoeError::count_out_of_range("n_center", 1..=MAX_CENTER_POINTS, usize::MAX)
        );
        let d = ccd(2, AlphaType::FaceCentered, MAX_CENTER_POINTS).unwrap();
        assert_eq!(d.data.len(), 4 + 4 + MAX_CENTER_POINTS);
    }

    /// The type of every run agrees with its coordinates: a center run is all
    /// zeros, an axial run has exactly one non-zero coordinate, a factorial
    /// run has none at zero -- for every alpha type and size.
    #[test]
    fn point_types_agree_with_the_coordinates() {
        for k in 2..=6 {
            for alpha in [
                AlphaType::FaceCentered,
                AlphaType::Rotatable,
                AlphaType::Inscribed,
            ] {
                let d = ccd_design(k, alpha, 3).unwrap();
                assert_eq!(d.point_types.len(), d.matrix.run_count());
                let counts = |t| d.point_types.iter().filter(|&&p| p == t).count();
                assert_eq!(counts(PointType::Factorial), 1 << k, "k={k}");
                assert_eq!(counts(PointType::Axial), 2 * k, "k={k}");
                assert_eq!(counts(PointType::Center), 3, "k={k}");
                for (row, t) in d.matrix.data.iter().zip(&d.point_types) {
                    let zeros = row.iter().filter(|&&x| x == 0.0).count();
                    match t {
                        PointType::Center => assert_eq!(zeros, k),
                        PointType::Axial => assert_eq!(zeros, k - 1),
                        PointType::Factorial => assert_eq!(zeros, 0),
                        PointType::Edge => unreachable!("a CCD has no edge points"),
                    }
                }
            }
        }
    }

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
                assert!(
                    (v.abs() - 1.0).abs() < 1e-10,
                    "factorial point run {run} has value {v}"
                );
            }
        }
    }

    #[test]
    fn ccd_axial_points_one_nonzero_each() {
        let d = ccd(3, AlphaType::FaceCentered, 1).unwrap();
        let factorial_runs = 1 << 3; // 8
        for i in 0..(2 * 3) {
            let run = factorial_runs + i;
            let nonzero: Vec<f64> = d.data[run]
                .iter()
                .copied()
                .filter(|&v| v.abs() > 1e-9)
                .collect();
            assert_eq!(
                nonzero.len(),
                1,
                "axial run {i} should have exactly 1 nonzero factor"
            );
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
            let nonzero: Vec<f64> = d.data[run]
                .iter()
                .copied()
                .filter(|&v| v.abs() > 1e-9)
                .collect();
            assert_eq!(nonzero.len(), 1);
            assert!(
                (nonzero[0].abs() - 1.0).abs() < 1e-10,
                "inscribed axial should be ±1"
            );
        }
    }
}
