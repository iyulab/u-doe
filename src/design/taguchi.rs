//! Taguchi orthogonal arrays.
//!
//! Provides standard Taguchi orthogonal array designs as lookup tables.
//! Coded values: ±1 for two-level factors, −1/0/+1 for three-level factors.
//!
//! Supported arrays: L4, L8, L9, L12, L16, L18, L27.
//!
//! Reference: Taguchi, G. (1987). *System of Experimental Design*.
//! UNIPUB/Kraus International.

use super::DesignMatrix;
use crate::error::DoeError;

struct ArraySpec {
    name: &'static str,
    runs: usize,
    max_factors: usize,
    /// Stored as 1-based level indices; converted to coded values on output.
    /// Level count is inferred per-column from the maximum value in that column:
    ///   max == 2 → two-level:   1 → −1,  2 → +1
    ///   max == 3 → three-level: 1 → −1,  2 →  0,  3 → +1
    columns: &'static [&'static [u8]],
}

// -----------------------------------------------------------------------
// Array definitions (column-major: each inner slice is one column)
// Source: Taguchi (1987)
// -----------------------------------------------------------------------

// L4(2^3): 4 runs, 3 two-level columns
const L4_COLS: &[&[u8]] = &[&[1, 1, 2, 2], &[1, 2, 1, 2], &[1, 2, 2, 1]];

// L8(2^7): 8 runs, 7 two-level columns
const L8_COLS: &[&[u8]] = &[
    &[1, 1, 1, 1, 2, 2, 2, 2],
    &[1, 1, 2, 2, 1, 1, 2, 2],
    &[1, 1, 2, 2, 2, 2, 1, 1],
    &[1, 2, 1, 2, 1, 2, 1, 2],
    &[1, 2, 1, 2, 2, 1, 2, 1],
    &[1, 2, 2, 1, 1, 2, 2, 1],
    &[1, 2, 2, 1, 2, 1, 1, 2],
];

// L9(3^4): 9 runs, 4 three-level columns
const L9_COLS: &[&[u8]] = &[
    &[1, 1, 1, 2, 2, 2, 3, 3, 3],
    &[1, 2, 3, 1, 2, 3, 1, 2, 3],
    &[1, 2, 3, 2, 3, 1, 3, 1, 2],
    &[1, 2, 3, 3, 1, 2, 2, 3, 1],
];

// L12(2^11): 12 runs, 11 two-level columns (Plackett-Burman)
const L12_COLS: &[&[u8]] = &[
    &[1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2],
    &[1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2],
    &[1, 1, 2, 1, 2, 2, 1, 2, 2, 1, 1, 2],
    &[1, 1, 2, 2, 1, 2, 2, 1, 2, 1, 2, 1],
    &[1, 1, 2, 2, 2, 1, 2, 2, 1, 2, 1, 1],
    &[1, 2, 1, 1, 2, 2, 1, 2, 2, 2, 1, 1],
    &[1, 2, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1],
    &[1, 2, 1, 2, 2, 1, 2, 2, 1, 1, 1, 2],
    &[1, 2, 2, 1, 1, 2, 2, 2, 1, 1, 2, 1],
    &[1, 2, 2, 1, 2, 1, 1, 1, 2, 2, 2, 1],
    &[1, 2, 2, 2, 1, 1, 1, 2, 1, 2, 1, 2],
];

// L16(2^15): 16 runs, 15 two-level columns (full 2^4 factorial + interactions)
const L16_COLS: &[&[u8]] = &[
    &[1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2],
    &[1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2],
    &[1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1],
    &[1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2],
    &[1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1],
    &[1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1],
    &[1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2],
    &[1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2],
    &[1, 2, 1, 2, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1],
    &[1, 2, 1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1],
    &[1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2],
    &[1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1],
    &[1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2],
    &[1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 2],
    &[1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 2, 2, 1],
];

// L18(2^1 × 3^7): 18 runs, 8 columns.
// Column 0: two-level (max = 2 → coded as ±1).
// Columns 1-7: three-level (max = 3 → coded as −1/0/+1).
// Level count is automatically inferred per column via `infer_levels`.
const L18_COLS: &[&[u8]] = &[
    &[1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2], // 2-level
    &[1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3],
    &[1, 1, 2, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3],
    &[1, 2, 3, 1, 2, 3, 2, 3, 1, 3, 1, 2, 2, 3, 1, 3, 1, 2],
    &[1, 2, 3, 2, 3, 1, 3, 1, 2, 2, 3, 1, 1, 2, 3, 3, 1, 2],
    &[1, 2, 3, 3, 1, 2, 1, 2, 3, 3, 1, 2, 3, 1, 2, 2, 3, 1],
    &[1, 3, 2, 1, 3, 2, 2, 1, 3, 2, 1, 3, 3, 2, 1, 1, 3, 2],
    &[1, 3, 2, 3, 2, 1, 3, 2, 1, 1, 3, 2, 2, 1, 3, 2, 1, 3],
];

// L27(3^13): 27 runs, 13 three-level columns
const L27_COLS: &[&[u8]] = &[
    &[
        1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    ],
    &[
        1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3,
    ],
    &[
        1, 1, 1, 2, 2, 2, 3, 3, 3, 2, 2, 2, 3, 3, 3, 1, 1, 1, 3, 3, 3, 1, 1, 1, 2, 2, 2,
    ],
    &[
        1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3,
    ],
    &[
        1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3, 1, 2, 3, 1, 2, 3, 1, 3, 1, 2, 3, 1, 2, 3, 1, 2,
    ],
    &[
        1, 2, 3, 2, 3, 1, 3, 1, 2, 1, 2, 3, 2, 3, 1, 3, 1, 2, 2, 3, 1, 3, 1, 2, 1, 2, 3,
    ],
    &[
        1, 2, 3, 2, 3, 1, 3, 1, 2, 2, 3, 1, 3, 1, 2, 1, 2, 3, 3, 1, 2, 1, 2, 3, 2, 3, 1,
    ],
    &[
        1, 2, 3, 3, 1, 2, 2, 3, 1, 1, 2, 3, 3, 1, 2, 2, 3, 1, 1, 2, 3, 3, 1, 2, 2, 3, 1,
    ],
    &[
        1, 2, 3, 3, 1, 2, 2, 3, 1, 2, 3, 1, 1, 2, 3, 3, 1, 2, 2, 3, 1, 1, 2, 3, 3, 1, 2,
    ],
    &[
        1, 2, 3, 3, 1, 2, 2, 3, 1, 3, 1, 2, 2, 3, 1, 1, 2, 3, 3, 1, 2, 2, 3, 1, 1, 2, 3,
    ],
    &[
        1, 3, 2, 1, 3, 2, 1, 3, 2, 2, 1, 3, 2, 1, 3, 2, 1, 3, 3, 2, 1, 3, 2, 1, 3, 2, 1,
    ],
    &[
        1, 3, 2, 2, 1, 3, 3, 2, 1, 1, 3, 2, 2, 1, 3, 3, 2, 1, 1, 3, 2, 2, 1, 3, 3, 2, 1,
    ],
    &[
        1, 3, 2, 3, 2, 1, 2, 1, 3, 1, 3, 2, 3, 2, 1, 2, 1, 3, 2, 1, 3, 1, 3, 2, 3, 2, 1,
    ],
];

const ARRAYS: &[ArraySpec] = &[
    ArraySpec {
        name: "L4",
        runs: 4,
        max_factors: 3,
        columns: L4_COLS,
    },
    ArraySpec {
        name: "L8",
        runs: 8,
        max_factors: 7,
        columns: L8_COLS,
    },
    ArraySpec {
        name: "L9",
        runs: 9,
        max_factors: 4,
        columns: L9_COLS,
    },
    ArraySpec {
        name: "L12",
        runs: 12,
        max_factors: 11,
        columns: L12_COLS,
    },
    ArraySpec {
        name: "L16",
        runs: 16,
        max_factors: 15,
        columns: L16_COLS,
    },
    ArraySpec {
        name: "L18",
        runs: 18,
        max_factors: 8,
        columns: L18_COLS,
    },
    ArraySpec {
        name: "L27",
        runs: 27,
        max_factors: 13,
        columns: L27_COLS,
    },
];

/// Infer the number of levels for a column from its maximum stored value.
///
/// Two-level columns use values 1 and 2 (max == 2).
/// Three-level columns use values 1, 2, and 3 (max == 3).
fn infer_levels(col_data: &[u8]) -> usize {
    let max = col_data
        .iter()
        .copied()
        .max()
        .expect("column must be non-empty");
    if max <= 2 {
        2
    } else {
        3
    }
}

/// Convert a 1-based level index to a coded value.
///
/// Two-level:   1 → −1.0,  2 → +1.0
/// Three-level: 1 → −1.0,  2 →  0.0,  3 → +1.0
fn to_coded(level: u8, n_levels: usize) -> f64 {
    match n_levels {
        2 => {
            if level == 1 {
                -1.0
            } else {
                1.0
            }
        }
        3 => match level {
            1 => -1.0,
            2 => 0.0,
            _ => 1.0,
        },
        _ => unreachable!("only 2- and 3-level arrays are supported"),
    }
}

/// Generate a Taguchi orthogonal array design.
///
/// # Arguments
/// * `name` — Array name: `"L4"`, `"L8"`, `"L9"`, `"L12"`, `"L16"`, `"L18"`, `"L27"`
/// * `k`    — Number of factors to include (1 ≤ k ≤ max_factors for the array)
///
/// # Errors
///
/// Returns [`DoeError::UnsupportedDesign`] if `name` is not recognised.
/// Returns [`DoeError::InvalidFactorCount`] if `k == 0` or `k > max_factors`.
///
/// # Examples
///
/// ```
/// use u_doe::design::taguchi::taguchi_array;
/// let d = taguchi_array("L8", 7).unwrap();
/// assert_eq!(d.run_count(), 8);
/// assert_eq!(d.factor_count(), 7);
/// ```
pub fn taguchi_array(name: &str, k: usize) -> Result<DesignMatrix, DoeError> {
    let spec = ARRAYS.iter().find(|a| a.name == name).ok_or_else(|| {
        DoeError::UnsupportedDesign(format!(
            "unknown Taguchi array '{name}'; supported: L4, L8, L9, L12, L16, L18, L27"
        ))
    })?;

    if k == 0 || k > spec.max_factors {
        return Err(DoeError::InvalidFactorCount {
            min: 1,
            max: spec.max_factors,
            got: k,
        });
    }

    // Build run-major data matrix from column-major storage.
    // Level count is inferred per-column to handle mixed designs (e.g. L18).
    let data: Vec<Vec<f64>> = (0..spec.runs)
        .map(|run| {
            (0..k)
                .map(|col| {
                    let n_levels = infer_levels(spec.columns[col]);
                    to_coded(spec.columns[col][run], n_levels)
                })
                .collect()
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
    fn l4_structure() {
        let d = taguchi_array("L4", 3).unwrap();
        assert_eq!(d.run_count(), 4);
        assert_eq!(d.factor_count(), 3);
    }

    #[test]
    fn l8_structure() {
        let d = taguchi_array("L8", 7).unwrap();
        assert_eq!(d.run_count(), 8);
        assert_eq!(d.factor_count(), 7);
    }

    #[test]
    fn l9_three_level_values() {
        let d = taguchi_array("L9", 4).unwrap();
        assert_eq!(d.run_count(), 9);
        for row in &d.data {
            for &v in row {
                assert!(
                    (v + 1.0).abs() < 1e-9 || v.abs() < 1e-9 || (v - 1.0).abs() < 1e-9,
                    "unexpected value {v}"
                );
            }
        }
    }

    #[test]
    fn l8_orthogonality() {
        let d = taguchi_array("L8", 7).unwrap();
        for i in 0..7 {
            for j in (i + 1)..7 {
                let dot: f64 = (0..8).map(|r| d.data[r][i] * d.data[r][j]).sum();
                assert!(dot.abs() < 1e-9, "L8 cols {i},{j} not orthogonal");
            }
        }
    }

    #[test]
    fn l9_orthogonality() {
        let d = taguchi_array("L9", 4).unwrap();
        // For 3-level orthogonal arrays, each pair of columns has zero inner product
        for i in 0..4 {
            for j in (i + 1)..4 {
                let dot: f64 = (0..9).map(|r| d.data[r][i] * d.data[r][j]).sum();
                assert!(dot.abs() < 1e-9, "L9 cols {i},{j} not orthogonal");
            }
        }
    }

    #[test]
    fn l16_structure() {
        let d = taguchi_array("L16", 15).unwrap();
        assert_eq!(d.run_count(), 16);
        assert_eq!(d.factor_count(), 15);
    }

    #[test]
    fn l27_structure() {
        let d = taguchi_array("L27", 13).unwrap();
        assert_eq!(d.run_count(), 27);
        assert_eq!(d.factor_count(), 13);
    }

    #[test]
    fn too_many_factors_error() {
        // L4 max 3 factors
        assert!(taguchi_array("L4", 4).is_err());
    }

    #[test]
    fn zero_factors_error() {
        assert!(taguchi_array("L8", 0).is_err());
    }

    #[test]
    fn unknown_array_error() {
        assert!(taguchi_array("L6", 3).is_err());
    }

    #[test]
    fn truncate_to_k_factors() {
        // Request fewer factors than available columns
        let d = taguchi_array("L8", 3).unwrap();
        assert_eq!(d.factor_count(), 3);
        assert_eq!(d.run_count(), 8);
    }

    #[test]
    fn l18_mixed_levels() {
        let d = taguchi_array("L18", 8).unwrap();
        assert_eq!(d.run_count(), 18);
        assert_eq!(d.factor_count(), 8);
        // Column 0 must be two-level: only ±1
        for r in 0..18 {
            let v = d.data[r][0];
            assert!(
                (v + 1.0).abs() < 1e-9 || (v - 1.0).abs() < 1e-9,
                "L18 col 0 run {r}: expected ±1, got {v}"
            );
        }
        // Columns 1-7 must be three-level: only −1, 0, +1
        for col in 1..8 {
            for r in 0..18 {
                let v = d.data[r][col];
                assert!(
                    (v + 1.0).abs() < 1e-9 || v.abs() < 1e-9 || (v - 1.0).abs() < 1e-9,
                    "L18 col {col} run {r}: unexpected value {v}"
                );
            }
        }
    }

    #[test]
    fn l12_structure() {
        let d = taguchi_array("L12", 11).unwrap();
        assert_eq!(d.run_count(), 12);
        assert_eq!(d.factor_count(), 11);
    }
}
