//! Box-Behnken Design (BBD).
//!
//! Three-level design for response surface modeling that avoids vertex
//! (corner) points. Uses balanced incomplete block design (BIB) structure
//! so each factor pair appears at ±1 while remaining factors are 0.
//!
//! Supported sizes:
//! - k=3: 15 runs (12 edge midpoints + 3 center)
//! - k=4: 27 runs (24 edge midpoints + 3 center)
//! - k=5: 46 runs (40 edge midpoints + 6 center)
//!
//! Reference: Box, G.E.P. & Behnken, D.W. (1960). "Some New Three Level
//! Designs for the Study of Quantitative Variables".
//! *Technometrics* 2(4), pp. 455–475.

use super::{ClassifiedDesign, DesignMatrix, PointType};
use crate::error::DoeError;

/// Generate a Box-Behnken design for `k` factors (k = 3, 4, or 5).
///
/// # Arguments
/// * `k`        — Number of factors (3, 4, or 5)
/// * `n_center` — Number of center point replicates (1 to [`MAX_CENTER_POINTS`](crate::design::MAX_CENTER_POINTS))
///
/// # Errors
///
/// Returns `Err` if `k` is not 3, 4, or 5, or if `n_center` is not in `1..=MAX_CENTER_POINTS`.
///
/// # Examples
///
/// ```
/// use u_doe::design::box_behnken::box_behnken;
/// let d = box_behnken(3, 3).unwrap();
/// assert_eq!(d.run_count(), 15);
/// assert_eq!(d.factor_count(), 3);
/// ```
pub fn box_behnken(k: usize, n_center: usize) -> Result<DesignMatrix, DoeError> {
    box_behnken_design(k, n_center).map(|d| d.matrix)
}

/// Generate a Box-Behnken design with the block of each run: the edge
/// midpoints first, then the center rows -- the order of the matrix.
///
/// Arguments and errors as for [`box_behnken`].
///
/// # Examples
///
/// ```
/// use u_doe::design::box_behnken::box_behnken_design;
/// use u_doe::design::PointType;
/// let d = box_behnken_design(3, 3).unwrap();
/// assert_eq!(d.point_types.iter().filter(|&&t| t == PointType::Edge).count(), 12);
/// assert_eq!(d.point_types.iter().filter(|&&t| t == PointType::Center).count(), 3);
/// ```
pub fn box_behnken_design(k: usize, n_center: usize) -> Result<ClassifiedDesign, DoeError> {
    if !matches!(k, 3..=5) {
        return Err(DoeError::count_out_of_range("k", 3..=5, k));
    }
    if !(1..=super::MAX_CENTER_POINTS).contains(&n_center) {
        return Err(DoeError::count_out_of_range(
            "n_center",
            1..=super::MAX_CENTER_POINTS,
            n_center,
        ));
    }

    let pairs = bib_pairs(k);
    let mut data: Vec<Vec<f64>> = Vec::new();

    // For each pair (i, j), generate 4 runs: [±1, ±1] at positions i and j, 0 elsewhere
    for &(i, j) in &pairs {
        for &si in &[-1.0_f64, 1.0] {
            for &sj in &[-1.0_f64, 1.0] {
                let mut row = vec![0.0_f64; k];
                row[i] = si;
                row[j] = sj;
                data.push(row);
            }
        }
    }

    let mut point_types = vec![PointType::Edge; data.len()];

    // Center points
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

/// BIB (Balanced Incomplete Block) pairs for Box-Behnken designs.
///
/// Returns the factor index pairs to use for edge-midpoint runs.
/// Source: Box & Behnken (1960), Tables 1–3.
fn bib_pairs(k: usize) -> Vec<(usize, usize)> {
    match k {
        3 => vec![(0, 1), (0, 2), (1, 2)],
        4 => vec![(0, 1), (2, 3), (0, 2), (1, 3), (0, 3), (1, 2)],
        5 => {
            // All C(5,2) = 10 pairs
            let mut pairs = Vec::new();
            for i in 0..5 {
                for j in (i + 1)..5 {
                    pairs.push((i, j));
                }
            }
            pairs
        }
        _ => unreachable!("k already validated as 3, 4, or 5"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every run's type agrees with its coordinates: an edge midpoint has
    /// exactly two non-zero coordinates, a center run none.
    #[test]
    fn point_types_agree_with_the_coordinates() {
        for k in 3..=5 {
            let d = box_behnken_design(k, 3).unwrap();
            assert_eq!(d.point_types.len(), d.matrix.run_count());
            let centers = d
                .point_types
                .iter()
                .filter(|&&t| t == PointType::Center)
                .count();
            assert_eq!(centers, 3, "k={k}");
            for (row, t) in d.matrix.data.iter().zip(&d.point_types) {
                let nonzero = row.iter().filter(|&&x| x != 0.0).count();
                match t {
                    PointType::Center => assert_eq!(nonzero, 0),
                    PointType::Edge => assert_eq!(nonzero, 2),
                    other => panic!("a Box-Behnken design has no {other:?} points"),
                }
            }
        }
    }

    #[test]
    fn bbd_k3_runs() {
        let d = box_behnken(3, 3).unwrap();
        assert_eq!(d.run_count(), 15);
        assert_eq!(d.factor_count(), 3);
    }

    #[test]
    fn bbd_k4_runs() {
        let d = box_behnken(4, 3).unwrap();
        assert_eq!(d.run_count(), 27);
    }

    #[test]
    fn bbd_k5_runs() {
        let d = box_behnken(5, 6).unwrap();
        assert_eq!(d.run_count(), 46);
    }

    #[test]
    fn bbd_no_corner_points() {
        // No run should have ALL factors simultaneously at ±1
        let d = box_behnken(3, 1).unwrap();
        for row in &d.data {
            let at_extreme = row.iter().filter(|&&v| v.abs() > 0.9).count();
            assert!(at_extreme <= 2, "corner point detected: {row:?}");
        }
    }

    #[test]
    fn bbd_values_in_range() {
        let d = box_behnken(4, 1).unwrap();
        for row in &d.data {
            for &v in row {
                assert!(v.abs() <= 1.0 + 1e-10, "value {v} out of [-1,1]");
            }
        }
    }

    #[test]
    fn bbd_center_points_all_zero() {
        let d = box_behnken(3, 3).unwrap();
        // Last 3 rows are center points
        let n = d.run_count();
        for i in (n - 3)..n {
            for &v in &d.data[i] {
                assert!(v.abs() < 1e-10, "center point not zero: {v}");
            }
        }
    }

    #[test]
    fn bbd_k2_error() {
        assert!(box_behnken(2, 1).is_err());
    }

    #[test]
    fn bbd_k6_error() {
        assert!(box_behnken(6, 1).is_err());
    }
}
