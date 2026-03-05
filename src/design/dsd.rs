//! Definitive Screening Design (DSD).
//!
//! Generates a DSD with 2k+1 runs for k factors at three levels {-1, 0, +1}.
//! For even k, main effects are pairwise orthogonal (column dot products = 0),
//! achieved via Paley conference matrices. For odd k, the design preserves the
//! fold-over structure and 3-level property but strict column orthogonality
//! requires a larger conference matrix (not implemented here).
//!
//! # Reference
//! Jones, B. & Nachtsheim, C.J. (2011). "A Class of Three-Level Designs for
//! Definitive Screening in the Presence of Second-Order Effects."
//! *Journal of Quality Technology*, 43(1), 1–15.

use super::DesignMatrix;
use crate::error::DoeError;

/// Build a Definitive Screening Design for `k` continuous factors.
///
/// Supports `k = 2..=12`. Returns a [`DesignMatrix`] with `2k+1` rows and `k`
/// columns. All values are in `{-1, 0, +1}`. The last row is the centre point
/// (all zeros). The first `k` rows form a conference sub-matrix `C`; the next
/// `k` rows form its fold-over `-C`.
///
/// For even `k`, columns are pairwise orthogonal over the full `2k+1` run
/// design. For odd `k`, the fold-over balance and 3-level structure are
/// guaranteed but full column orthogonality requires a larger conference order.
///
/// # Errors
/// Returns [`DoeError::InvalidFactorCount`] when `k < 2` or `k > 12`.
pub fn definitive_screening(k: usize) -> Result<DesignMatrix, DoeError> {
    if !(2..=12).contains(&k) {
        return Err(DoeError::InvalidFactorCount {
            min: 2,
            max: 12,
            got: k,
        });
    }

    let conf = conference_matrix(k);
    let mut data: Vec<Vec<f64>> = Vec::with_capacity(2 * k + 1);

    for row in &conf {
        data.push(row.iter().map(|&v| v as f64).collect());
    }
    for row in &conf {
        data.push(row.iter().map(|&v| -(v as f64)).collect());
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

/// Returns a k×k matrix with entries in {-1, 0, +1} used as the DSD block.
/// For even k, this is a proper conference matrix (C'C = (k-1)I, columns
/// pairwise orthogonal). For odd k, a balanced approximation is used.
fn conference_matrix(k: usize) -> Vec<Vec<i8>> {
    match k {
        // --- k=2: trivial conference matrix C(2) ---
        2 => vec![vec![0, 1], vec![1, 0]],
        // --- k=3: best 3×3 fold-over block (not fully column-orthogonal) ---
        3 => vec![vec![0, 1, 1], vec![1, 0, -1], vec![-1, 1, 0]],
        // --- k=4: C(4) via Paley, column-orthogonal ---
        4 => vec![
            vec![0, -1, -1, -1],
            vec![-1, 0, -1, 1],
            vec![-1, 1, 0, -1],
            vec![-1, -1, 1, 0],
        ],
        // --- k=5: best 5×5 fold-over block (not fully column-orthogonal) ---
        5 => vec![
            vec![0, 1, 1, 1, 1],
            vec![1, 0, 1, -1, -1],
            vec![-1, 1, 0, -1, 1],
            vec![-1, -1, 1, 0, 1],
            vec![1, -1, -1, 1, 0],
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
        // --- k=7: best 7×7 fold-over block (not fully column-orthogonal) ---
        7 => cyclic_conference_matrix(7, &[0, 1, 1, 1, -1, -1, 1]),
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
        // --- k=9: best 9×9 fold-over block (not fully column-orthogonal) ---
        9 => cyclic_conference_matrix(9, &[0, 1, 1, -1, 1, -1, 1, 1, -1]),
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
        // --- k=11: best 11×11 fold-over block (not fully column-orthogonal) ---
        11 => cyclic_conference_matrix(11, &[0, 1, 1, -1, 1, -1, -1, 1, 1, 1, -1]),
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

/// Cyclic construction for odd k fallback:
/// `mat[i][j] = base[(j + k - i) % k]`, then `mat[i][i] = 0`.
fn cyclic_conference_matrix(k: usize, base: &[i8]) -> Vec<Vec<i8>> {
    let mut mat = vec![vec![0i8; k]; k];
    for i in 0..k {
        for j in 0..k {
            mat[i][j] = base[(j + k - i) % k];
        }
        mat[i][i] = 0;
    }
    mat
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_counts_2k_plus_1() {
        for k in 2..=8 {
            let d = definitive_screening(k).unwrap();
            assert_eq!(d.run_count(), 2 * k + 1, "k={k}");
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
        for k in 2..=6 {
            let d = definitive_screening(k).unwrap();
            for i in 0..k {
                for c in 0..k {
                    assert!(
                        (d.get(i, c) + d.get(i + k, c)).abs() < 1e-10,
                        "k={k} i={i} c={c}: {} + {} != 0",
                        d.get(i, c),
                        d.get(i + k, c)
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

    /// For even k, proper Paley conference matrices guarantee pairwise orthogonal
    /// columns in the full 2k+1 run design.
    #[test]
    fn main_effects_orthogonal_even_k() {
        for k in [2usize, 4, 6, 8, 10, 12] {
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
