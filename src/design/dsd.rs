//! Definitive Screening Design (DSD).
//!
//! Generates a DSD for k factors at three levels {-1, 0, +1}: 2k+1 runs for
//! even k, 2k+3 for odd k. Main effects are pairwise orthogonal (column dot
//! products = 0), achieved via Paley conference matrices. For odd k the
//! conference matrix of order k+1 is used and its last column dropped -- the
//! "fake factor" construction of Jones & Nachtsheim (2011), which keeps the
//! columns orthogonal at the price of two more runs.
//!
//! # Reference
//! Jones, B. & Nachtsheim, C.J. (2011). "A Class of Three-Level Designs for
//! Definitive Screening in the Presence of Second-Order Effects."
//! *Journal of Quality Technology*, 43(1), 1–15.

use super::DesignMatrix;
use crate::error::DoeError;

/// Build a Definitive Screening Design for `k` continuous factors.
///
/// Supports `k = 2..=12`. Returns a [`DesignMatrix`] with `k` columns and
/// `2k+1` rows for even `k`, `2k+3` rows for odd `k`. All values are in
/// `{-1, 0, +1}`. The last row is the centre point (all zeros). The first
/// `m` rows form a conference matrix `C` (`m = k` for even `k`, `k+1` for odd
/// `k`, whose last column is then dropped); the next `m` rows form its
/// fold-over `-C`.
///
/// The columns are pairwise orthogonal over the whole design, and every main
/// effect is orthogonal to every quadratic and two-factor-interaction column.
///
/// # Errors
/// Returns [`DoeError::ParameterOutOfRange`] when `k < 2` or `k > 12`.
pub fn definitive_screening(k: usize) -> Result<DesignMatrix, DoeError> {
    if !(2..=12).contains(&k) {
        return Err(DoeError::count_out_of_range("k", 2..=12, k));
    }

    // An odd k borrows the next even order and drops the surplus column; the
    // rows that had a 0 there keep their other entries at +-1.
    let m = k + k % 2;
    let conf = conference_matrix(m);
    let mut data: Vec<Vec<f64>> = Vec::with_capacity(2 * m + 1);

    for row in &conf {
        data.push(row[..k].iter().map(|&v| v as f64).collect());
    }
    for row in &conf {
        data.push(row[..k].iter().map(|&v| -(v as f64)).collect());
    }
    data.push(vec![0.0; k]);

    Ok(DesignMatrix {
        data,
        factor_names: DesignMatrix::default_names(k),
    })
}

// ---------------------------------------------------------------------------
// Conference matrix construction
// ---------------------------------------------------------------------------

/// Returns an even-order k×k conference matrix (C'C = (k-1)I, columns
/// pairwise orthogonal) with entries in {-1, 0, +1}.
fn conference_matrix(k: usize) -> Vec<Vec<i8>> {
    match k {
        // --- k=2: trivial conference matrix C(2) ---
        2 => vec![vec![0, 1], vec![1, 0]],
        // --- k=4: C(4) via Paley, column-orthogonal ---
        4 => vec![
            vec![0, -1, -1, -1],
            vec![-1, 0, -1, 1],
            vec![-1, 1, 0, -1],
            vec![-1, -1, 1, 0],
        ],
        // --- k=6: C(6) via Paley type I (q=5), column-orthogonal ---
        6 => vec![
            vec![0, 1, 1, 1, 1, 1],
            vec![1, 0, 1, -1, -1, 1],
            vec![1, 1, 0, 1, -1, -1],
            vec![1, -1, 1, 0, 1, -1],
            vec![1, -1, -1, 1, 0, 1],
            vec![1, 1, -1, -1, 1, 0],
        ],
        // --- k=8: C(8) via Paley (q=7, QR={1,2,4}), column-orthogonal ---
        8 => vec![
            vec![0, 1, 1, 1, 1, 1, 1, 1],
            vec![1, 0, 1, 1, -1, 1, -1, -1],
            vec![1, -1, 0, 1, 1, -1, 1, -1],
            vec![1, -1, -1, 0, 1, 1, -1, 1],
            vec![1, 1, -1, -1, 0, 1, 1, -1],
            vec![1, -1, 1, -1, -1, 0, 1, 1],
            vec![1, 1, -1, 1, -1, -1, 0, 1],
            vec![1, 1, 1, -1, 1, -1, -1, 0],
        ],
        // --- k=10: C(10) via GF(9) Paley, column-orthogonal ---
        10 => vec![
            vec![0, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            vec![1, 0, 1, 1, 1, -1, -1, 1, -1, -1],
            vec![1, 1, 0, 1, -1, 1, -1, -1, 1, -1],
            vec![1, 1, 1, 0, -1, -1, 1, -1, -1, 1],
            vec![1, 1, -1, -1, 0, 1, 1, 1, -1, -1],
            vec![1, -1, 1, -1, 1, 0, 1, -1, 1, -1],
            vec![1, -1, -1, 1, 1, 1, 0, -1, -1, 1],
            vec![1, 1, -1, -1, 1, -1, -1, 0, 1, 1],
            vec![1, -1, 1, -1, -1, 1, -1, 1, 0, 1],
            vec![1, -1, -1, 1, -1, -1, 1, 1, 1, 0],
        ],
        // --- k=12: C(12) via Paley-skew (q=11), column-orthogonal ---
        12 => vec![
            vec![0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            vec![-1, 0, 1, -1, 1, 1, 1, -1, -1, -1, 1, -1],
            vec![-1, -1, 0, 1, -1, 1, 1, 1, -1, -1, -1, 1],
            vec![-1, 1, -1, 0, 1, -1, 1, 1, 1, -1, -1, -1],
            vec![-1, -1, 1, -1, 0, 1, -1, 1, 1, 1, -1, -1],
            vec![-1, -1, -1, 1, -1, 0, 1, -1, 1, 1, 1, -1],
            vec![-1, -1, -1, -1, 1, -1, 0, 1, -1, 1, 1, 1],
            vec![-1, 1, -1, -1, -1, 1, -1, 0, 1, -1, 1, 1],
            vec![-1, 1, 1, -1, -1, -1, 1, -1, 0, 1, -1, 1],
            vec![-1, 1, 1, 1, -1, -1, -1, 1, -1, 0, 1, -1],
            vec![-1, -1, 1, 1, 1, -1, -1, -1, 1, -1, 0, 1],
            vec![-1, 1, -1, 1, 1, 1, -1, -1, -1, 1, -1, 0],
        ],
        _ => unreachable!("k={k} out of range 2..=12"),
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_counts_2k_plus_1_or_3() {
        for k in 2..=8 {
            let d = definitive_screening(k).unwrap();
            assert_eq!(d.run_count(), 2 * k + 1 + 2 * (k % 2), "k={k}");
            assert_eq!(d.factor_count(), k);
        }
    }

    #[test]
    fn three_levels_only() {
        for k in 2..=6 {
            let d = definitive_screening(k).unwrap();
            for r in 0..d.run_count() {
                for c in 0..k {
                    let v = d.get(r, c);
                    assert!(
                        (v - 0.0).abs() < 1e-10
                            || (v - 1.0).abs() < 1e-10
                            || (v + 1.0).abs() < 1e-10,
                        "k={k} r={r} c={c} v={v}"
                    );
                }
            }
        }
    }

    #[test]
    fn center_row_all_zeros() {
        for k in 2..=6 {
            let d = definitive_screening(k).unwrap();
            let last = d.run_count() - 1;
            for c in 0..k {
                assert!((d.get(last, c)).abs() < 1e-10, "k={k} c={c}");
            }
        }
    }

    #[test]
    fn foldover_symmetry() {
        for k in 2..=7 {
            let d = definitive_screening(k).unwrap();
            let m = k + k % 2;
            for i in 0..m {
                for c in 0..k {
                    assert!(
                        (d.get(i, c) + d.get(i + m, c)).abs() < 1e-10,
                        "k={k} i={i} c={c}: {} + {} != 0",
                        d.get(i, c),
                        d.get(i + m, c)
                    );
                }
            }
        }
    }

    #[test]
    fn main_effects_orthogonal_k4() {
        let d = definitive_screening(4).unwrap();
        let n = d.run_count();
        for j1 in 0..4 {
            for j2 in (j1 + 1)..4 {
                let dot: f64 = (0..n).map(|r| d.get(r, j1) * d.get(r, j2)).sum();
                assert!(
                    dot.abs() < 1e-10,
                    "cols {j1},{j2} not orthogonal: dot={dot}"
                );
            }
        }
    }

    #[test]
    fn invalid_k_error() {
        assert!(definitive_screening(0).is_err());
        assert!(definitive_screening(1).is_err());
        assert!(definitive_screening(13).is_err());
    }

    // -----------------------------------------------------------------------
    // Reference-value tests (Task 11)
    // Jones & Nachtsheim (2011), §2-3
    // -----------------------------------------------------------------------

    /// Run count 2k+1 for even k, 2k+3 for odd k (Jones & Nachtsheim 2011:
    /// the conference matrix of the next even order, one column dropped).
    #[test]
    fn run_counts_all_k() {
        for k in 2..=12 {
            let d = definitive_screening(k).unwrap();
            let expected = 2 * (k + k % 2) + 1;
            assert_eq!(d.run_count(), expected, "k={k}");
        }
        // Spot-check the task-specified values
        assert_eq!(definitive_screening(3).unwrap().run_count(), 9);
        assert_eq!(definitive_screening(4).unwrap().run_count(), 9);
        assert_eq!(definitive_screening(6).unwrap().run_count(), 13);
        assert_eq!(definitive_screening(12).unwrap().run_count(), 25);
    }

    /// Center point (last row) must be all zeros for every supported k.
    #[test]
    fn center_row_all_zeros_full_range() {
        for k in 2..=12 {
            let d = definitive_screening(k).unwrap();
            let last = d.run_count() - 1;
            for c in 0..k {
                assert!(
                    d.get(last, c).abs() < 1e-10,
                    "k={k} center row col={c} is not zero: {}",
                    d.get(last, c)
                );
            }
        }
    }

    /// Proper Paley conference matrices (an odd k borrows the next even order
    /// and drops a column) guarantee pairwise orthogonal columns.
    #[test]
    fn main_effects_orthogonal_every_k() {
        for k in 2..=12usize {
            let d = definitive_screening(k).unwrap();
            let n = d.run_count();
            for j1 in 0..k {
                for j2 in (j1 + 1)..k {
                    let dot: f64 = (0..n).map(|r| d.get(r, j1) * d.get(r, j2)).sum();
                    assert!(
                        dot.abs() < 1e-10,
                        "k={k} cols {j1},{j2} not orthogonal: dot={dot}"
                    );
                }
            }
        }
    }
}
