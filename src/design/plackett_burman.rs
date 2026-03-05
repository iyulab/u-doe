//! Plackett-Burman screening designs.
//!
//! Generates orthogonal screening designs for up to N-1 factors in N runs,
//! where N is a multiple of 4. Uses cyclic construction from a base row.
//!
//! # Algorithm
//!
//! 1. Start with a base row of length N-1.
//! 2. Generate N-1 rows by successive cyclic shifts of the base row.
//! 3. Append an all-(−1) row to complete the N×(N-1) Hadamard matrix.
//! 4. Take the first k columns for k factors.
//!
//! # References
//!
//! Plackett, R.L. & Burman, J.P. (1946). "The Design of Optimum
//! Multifactorial Experiments". *Biometrika* 33(4), pp. 305–325.

use super::DesignMatrix;
use crate::error::DoeError;

/// (N, base_row) pairs. Up to N-1 factors can be screened in N runs.
/// Source: Plackett & Burman (1946), Table 1.
const PB_BASES: &[(usize, &[i8])] = &[
    (8,  &[ 1,  1,  1, -1,  1, -1, -1]),
    (12, &[ 1,  1, -1,  1,  1,  1, -1, -1, -1,  1, -1]),
    (16, &[ 1,  1,  1,  1, -1,  1, -1,  1,  1, -1, -1,  1, -1, -1, -1]),
    (20, &[ 1,  1, -1, -1,  1,  1,  1,  1, -1,  1, -1,  1, -1, -1, -1, -1,  1,  1, -1]),
];

/// Generate a Plackett-Burman design for `k` factors (1 ≤ k ≤ 19).
///
/// Automatically selects the smallest N (multiple of 4) such that N − 1 ≥ k.
///
/// # Errors
///
/// Returns `Err` if `k == 0` or `k > 19`.
///
/// # Examples
///
/// ```
/// use u_doe::design::plackett_burman::plackett_burman;
/// let d = plackett_burman(7).unwrap();
/// assert_eq!(d.run_count(), 8);
/// assert_eq!(d.factor_count(), 7);
/// ```
pub fn plackett_burman(k: usize) -> Result<DesignMatrix, DoeError> {
    if k == 0 || k > 19 {
        return Err(DoeError::InvalidFactorCount { min: 1, max: 19, got: k });
    }

    let (n, base) = PB_BASES
        .iter()
        .find(|(n, _)| *n > k)
        .map(|(n, b)| (*n, *b))
        .ok_or(DoeError::InvalidFactorCount { min: 1, max: 19, got: k })?;

    let m = n - 1; // number of base columns

    // Build N rows: N-1 cyclic shifts + all-(-1) row
    let mut rows: Vec<Vec<f64>> = (0..m)
        .map(|shift| (0..m).map(|j| base[(j + shift) % m] as f64).collect())
        .collect();
    rows.push(vec![-1.0; m]);

    // Keep only first k columns
    let data: Vec<Vec<f64>> = rows
        .into_iter()
        .map(|mut row| {
            row.truncate(k);
            row
        })
        .collect();

    Ok(DesignMatrix {
        data,
        factor_names: DesignMatrix::default_names(k),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pb_8_run_7_factors() {
        let d = plackett_burman(7).unwrap();
        assert_eq!(d.run_count(), 8);
        assert_eq!(d.factor_count(), 7);
    }

    #[test]
    fn pb_12_run_11_factors() {
        let d = plackett_burman(11).unwrap();
        assert_eq!(d.run_count(), 12);
    }

    #[test]
    fn pb_16_run_15_factors() {
        let d = plackett_burman(15).unwrap();
        assert_eq!(d.run_count(), 16);
    }

    #[test]
    fn pb_20_run_19_factors() {
        let d = plackett_burman(19).unwrap();
        assert_eq!(d.run_count(), 20);
    }

    #[test]
    fn pb_auto_selects_smallest_n() {
        // 8 factors needs N=12 (since N=8 only fits 7 factors)
        let d = plackett_burman(8).unwrap();
        assert_eq!(d.run_count(), 12);
    }

    #[test]
    fn pb_orthogonality_8run() {
        let d = plackett_burman(7).unwrap();
        for i in 0..7 {
            for j in (i + 1)..7 {
                let dot: f64 = (0..8).map(|r| d.data[r][i] * d.data[r][j]).sum();
                assert!(dot.abs() < 1e-10, "cols {i},{j} dot={dot}");
            }
        }
    }

    #[test]
    fn pb_orthogonality_12run() {
        let d = plackett_burman(11).unwrap();
        for i in 0..11 {
            for j in (i + 1)..11 {
                let dot: f64 = (0..12).map(|r| d.data[r][i] * d.data[r][j]).sum();
                assert!(dot.abs() < 1e-10, "cols {i},{j} dot={dot}");
            }
        }
    }

    #[test]
    fn pb_values_are_plus_minus_one() {
        let d = plackett_burman(11).unwrap();
        for row in &d.data {
            for &v in row {
                assert!(
                    (v - 1.0).abs() < 1e-10 || (v + 1.0).abs() < 1e-10,
                    "unexpected value {v}"
                );
            }
        }
    }

    #[test]
    fn pb_k0_error() {
        assert!(plackett_burman(0).is_err());
    }

    #[test]
    fn pb_k20_error() {
        assert!(plackett_burman(20).is_err());
    }
}
