//! Full and fractional two-level factorial designs.
//!
//! # Full Factorial (2^k)
//!
//! Generates all 2^k combinations of k factors at ±1 levels in Yates order.
//! Factor j alternates sign every `2^j` runs (0-indexed j).
//!
//! Reference: Montgomery (2019), *Introduction to Statistical Quality Control*, 8th ed., Section 6.2
//!
//! # Fractional Factorial (2^(k−p))
//!
//! Uses standard generators from Montgomery (2019) Table 8.14.
//! Supported: k=4..7, p=1..3 (standard combinations only).
//!
//! # Resolution
//!
//! - `III`: no main effect aliased with another main effect, but 2FI aliased with main effects
//! - `IV`:  main effects clean; 2FI aliased with other 2FI
//! - `V`:   main effects and all 2FI estimable without aliasing

use super::DesignMatrix;
use crate::error::DoeError;

/// Resolution of a fractional factorial design.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Resolution {
    /// Resolution III: main effects aliased with two-factor interactions.
    III,
    /// Resolution IV: two-factor interactions aliased with each other.
    IV,
    /// Resolution V: main effects and all 2FI estimable.
    V,
}

impl core::fmt::Display for Resolution {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Resolution::III => "III",
            Resolution::IV => "IV",
            Resolution::V => "V",
        })
    }
}

/// Metadata for a fractional factorial design.
pub struct FractionalInfo {
    /// Total factor count.
    pub k: usize,
    /// Number of generators.
    pub p: usize,
    /// Design resolution.
    pub resolution: Resolution,
    /// Defining relation (compact notation, e.g. `"I=ABCD"`).
    pub defining_relation: String,
    /// Generator equations in factor letters, e.g. `["E=ABC", "F=BCD", "G=ACD"]`.
    ///
    /// One entry per derived factor, in factor order. Lets consumers derive
    /// the alias structure of the fraction they actually received instead of
    /// pairing the design with an external published table.
    pub generators: Vec<String>,
}

// (k, p, base_generators[generator_index][source_col_indices], resolution, defining_relation)
// Generator column is the product of the listed base columns.
// Source: Montgomery (2019), Table 8.14
struct GeneratorEntry {
    k: usize,
    p: usize,
    generators: &'static [&'static [usize]], // each slice = product of base cols at those indices
    resolution: Resolution,
    defining_relation: &'static str,
}

const GENERATOR_TABLE: &[GeneratorEntry] = &[
    // 2^(4-1): 8 runs, D = ABC → Resolution IV
    GeneratorEntry {
        k: 4,
        p: 1,
        generators: &[&[0, 1, 2]],
        resolution: Resolution::IV,
        defining_relation: "I=ABCD",
    },
    // 2^(5-1): 16 runs, E = ABCD → Resolution V
    GeneratorEntry {
        k: 5,
        p: 1,
        generators: &[&[0, 1, 2, 3]],
        resolution: Resolution::V,
        defining_relation: "I=ABCDE",
    },
    // 2^(5-2): 8 runs, D=AB, E=AC → Resolution III
    GeneratorEntry {
        k: 5,
        p: 2,
        generators: &[&[0, 1], &[0, 2]],
        resolution: Resolution::III,
        defining_relation: "I=ABD=ACE=BCDE",
    },
    // 2^(6-1): 32 runs, F = ABCDE → Resolution VI (treated as V in enum)
    GeneratorEntry {
        k: 6,
        p: 1,
        generators: &[&[0, 1, 2, 3, 4]],
        resolution: Resolution::V,
        defining_relation: "I=ABCDEF",
    },
    // 2^(6-2): 16 runs, E=ABC, F=BCD → Resolution IV
    GeneratorEntry {
        k: 6,
        p: 2,
        generators: &[&[0, 1, 2], &[1, 2, 3]],
        resolution: Resolution::IV,
        defining_relation: "I=ABCE=BCDF=ADEF",
    },
    // 2^(6-3): 8 runs, D=AB, E=AC, F=BC → Resolution III
    // Standard generators per Montgomery (2019) Table 8.14 / NIST e-Handbook §5.3.3.4.7.
    GeneratorEntry {
        k: 6,
        p: 3,
        generators: &[&[0, 1], &[0, 2], &[1, 2]],
        resolution: Resolution::III,
        defining_relation: "I=ABD=ACE=BCF=BCDE=ACDF=ABEF=DEF",
    },
    // 2^(7-1): 64 runs, G = ABCDEF → Resolution VII (treated as V)
    GeneratorEntry {
        k: 7,
        p: 1,
        generators: &[&[0, 1, 2, 3, 4, 5]],
        resolution: Resolution::V,
        defining_relation: "I=ABCDEFG",
    },
    // 2^(7-2): 32 runs, F=ABCD, G=ABDE → Resolution IV
    GeneratorEntry {
        k: 7,
        p: 2,
        generators: &[&[0, 1, 2, 3], &[0, 1, 3, 4]],
        resolution: Resolution::IV,
        defining_relation: "I=ABCDF=ABDEG=CEFG",
    },
    // 2^(7-3): 16 runs, E=ABC, F=BCD, G=ACD → Resolution IV
    // Standard generators per Montgomery (2019) Table 8.14 / NIST e-Handbook §5.3.3.4.7.
    GeneratorEntry {
        k: 7,
        p: 3,
        generators: &[&[0, 1, 2], &[1, 2, 3], &[0, 2, 3]],
        resolution: Resolution::IV,
        defining_relation: "I=ABCE=BCDF=ACDG=ADEF=BDEG=ABFG=CEFG",
    },
];

/// Generate a 2^k full factorial design in Yates order.
///
/// # Arguments
/// * `k` — Number of factors (2..=7). Each factor has two levels: −1 and +1.
///
/// # Errors
/// Returns `Err` if `k < 2` or `k > 7`.
///
/// # Examples
///
/// ```
/// use u_doe::design::factorial::full_factorial;
/// let d = full_factorial(3).unwrap();
/// assert_eq!(d.run_count(), 8);
/// assert_eq!(d.factor_count(), 3);
/// ```
pub fn full_factorial(k: usize) -> Result<DesignMatrix, DoeError> {
    if !(2..=7).contains(&k) {
        return Err(DoeError::InvalidFactorCount {
            min: 2,
            max: 7,
            got: k,
        });
    }
    let n = 1_usize << k; // 2^k
    let mut data = vec![vec![0.0_f64; k]; n];
    for (run, row) in data.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            // Yates order: factor j alternates sign every 2^j runs
            *cell = if (run >> j) & 1 == 0 { -1.0 } else { 1.0 };
        }
    }
    Ok(DesignMatrix {
        data,
        factor_names: DesignMatrix::default_names(k),
    })
}

/// Generate a 2^(k-p) fractional factorial design.
///
/// Uses standard generators from Montgomery (2019) Table 8.14.
/// v0.1 supports: k=4..7, p=1..3 (standard combinations only).
///
/// # Arguments
/// * `k` — Total number of factors
/// * `p` — Number of generators (fraction = 2^(-p))
///
/// # Errors
/// Returns `Err` if the (k, p) combination is not in the standard table.
///
/// # Examples
///
/// ```
/// use u_doe::design::factorial::fractional_factorial;
/// let d = fractional_factorial(4, 1).unwrap();
/// assert_eq!(d.run_count(), 8);
/// assert_eq!(d.factor_count(), 4);
/// ```
pub fn fractional_factorial(k: usize, p: usize) -> Result<DesignMatrix, DoeError> {
    let entry = find_entry(k, p)?;
    let base_k = k - p;
    let base = full_factorial(base_k)?;

    // Build extended data: start with base columns, then append one column per
    // generator.  Each generated column is the element-wise product of the
    // listed base columns, computed directly on the mutable row iterator.
    let mut data: Vec<Vec<f64>> = base.data.clone();
    for gen in entry.generators {
        for row in data.iter_mut() {
            // Generated column = product of the specified base columns
            let val: f64 = gen.iter().map(|&col| row[col]).product();
            row.push(val);
        }
    }

    Ok(DesignMatrix {
        data,
        factor_names: DesignMatrix::default_names(k),
    })
}

/// Retrieve metadata (resolution, defining relation, generators) for a 2^(k-p) design.
///
/// # Errors
/// Returns `Err` if (k, p) is not in the standard table.
///
/// # Examples
///
/// ```
/// use u_doe::design::factorial::{fractional_factorial_info, Resolution};
/// let info = fractional_factorial_info(4, 1).unwrap();
/// assert_eq!(info.resolution, Resolution::IV);
/// assert_eq!(info.defining_relation, "I=ABCD");
/// assert_eq!(info.generators, vec!["D=ABC"]);
/// ```
pub fn fractional_factorial_info(k: usize, p: usize) -> Result<FractionalInfo, DoeError> {
    let entry = find_entry(k, p)?;
    // Derive the generator equations from the numeric column specs so the
    // strings can never drift from the matrix construction (single source).
    let base_k = k - p;
    let generators = entry
        .generators
        .iter()
        .enumerate()
        .map(|(i, cols)| {
            let derived = factor_letter(base_k + i);
            let word: String = cols.iter().map(|&c| factor_letter(c)).collect();
            format!("{derived}={word}")
        })
        .collect();
    Ok(FractionalInfo {
        k,
        p,
        resolution: entry.resolution,
        defining_relation: entry.defining_relation.to_string(),
        generators,
    })
}

/// Factor letter for a 0-based column index (0 → 'A'). Valid for k ≤ 7 designs.
fn factor_letter(index: usize) -> char {
    debug_assert!(index < 26, "factor index out of letter range");
    (b'A' + index as u8) as char
}

fn find_entry(k: usize, p: usize) -> Result<&'static GeneratorEntry, DoeError> {
    GENERATOR_TABLE
        .iter()
        .find(|e| e.k == k && e.p == p)
        .ok_or_else(|| {
            DoeError::UnsupportedDesign(format!(
                "2^({k}-{p}) not in standard table; supported (k,p): \
             (4,1),(5,1),(5,2),(6,1),(6,2),(6,3),(7,1),(7,2),(7,3)"
            ))
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn full_factorial_2k2_layout() {
        let d = full_factorial(2).unwrap();
        assert_eq!(d.run_count(), 4);
        assert_eq!(d.factor_count(), 2);
        // Yates order
        assert_eq!(d.data[0], vec![-1.0, -1.0]);
        assert_eq!(d.data[1], vec![1.0, -1.0]);
        assert_eq!(d.data[2], vec![-1.0, 1.0]);
        assert_eq!(d.data[3], vec![1.0, 1.0]);
    }

    #[test]
    fn full_factorial_run_counts() {
        assert_eq!(full_factorial(2).unwrap().run_count(), 4);
        assert_eq!(full_factorial(3).unwrap().run_count(), 8);
        assert_eq!(full_factorial(4).unwrap().run_count(), 16);
        assert_eq!(full_factorial(5).unwrap().run_count(), 32);
        assert_eq!(full_factorial(6).unwrap().run_count(), 64);
        assert_eq!(full_factorial(7).unwrap().run_count(), 128);
    }

    #[test]
    fn full_factorial_k1_error() {
        assert!(full_factorial(1).is_err());
    }

    #[test]
    fn full_factorial_k8_error() {
        assert!(full_factorial(8).is_err());
    }

    #[test]
    fn full_factorial_orthogonality() {
        // All column dot products must be 0 for 2^k
        let d = full_factorial(3).unwrap();
        for i in 0..3 {
            for j in (i + 1)..3 {
                let dot: f64 = (0..8).map(|r| d.data[r][i] * d.data[r][j]).sum();
                assert!(dot.abs() < 1e-10, "cols {i},{j} not orthogonal");
            }
        }
    }

    #[test]
    fn full_factorial_default_names() {
        let d = full_factorial(3).unwrap();
        assert_eq!(d.factor_names, vec!["A", "B", "C"]);
    }

    #[test]
    fn fractional_factorial_4_1_runs() {
        // 2^(4-1) = 8 runs
        let d = fractional_factorial(4, 1).unwrap();
        assert_eq!(d.run_count(), 8);
        assert_eq!(d.factor_count(), 4);
    }

    #[test]
    fn fractional_factorial_4_1_resolution_iv() {
        let info = fractional_factorial_info(4, 1).unwrap();
        assert_eq!(info.resolution, Resolution::IV);
    }

    #[test]
    fn fractional_factorial_5_2_runs() {
        // 2^(5-2) = 8 runs
        let d = fractional_factorial(5, 2).unwrap();
        assert_eq!(d.run_count(), 8);
        assert_eq!(d.factor_count(), 5);
    }

    #[test]
    fn fractional_factorial_5_2_resolution_iii() {
        let info = fractional_factorial_info(5, 2).unwrap();
        assert_eq!(info.resolution, Resolution::III);
    }

    #[test]
    fn fractional_factorial_unsupported_error() {
        assert!(fractional_factorial(3, 1).is_err()); // not in standard table
        assert!(fractional_factorial(4, 2).is_err()); // not supported
    }

    #[test]
    fn fractional_factorial_6_3_standard_generators() {
        // 2^(6-3)_III, 8 runs: D=AB, E=AC, F=BC (NIST §5.3.3.4.7, Montgomery
        // Table 8.14).
        let d = fractional_factorial(6, 3).unwrap();
        assert_eq!(d.run_count(), 8);
        assert_eq!(d.factor_count(), 6);
        for row in &d.data {
            assert_eq!(row[3], row[0] * row[1], "D != AB");
            assert_eq!(row[4], row[0] * row[2], "E != AC");
            assert_eq!(row[5], row[1] * row[2], "F != BC");
        }
        let info = fractional_factorial_info(6, 3).unwrap();
        assert_eq!(info.resolution, Resolution::III);
        assert_eq!(info.generators, vec!["D=AB", "E=AC", "F=BC"]);
    }

    #[test]
    fn fractional_factorial_orthogonality() {
        // 2^(4-1): first 3 base cols must be orthogonal
        let d = fractional_factorial(4, 1).unwrap();
        for i in 0..3 {
            for j in (i + 1)..3 {
                let dot: f64 = (0..8).map(|r| d.data[r][i] * d.data[r][j]).sum();
                assert!(dot.abs() < 1e-10, "base cols {i},{j} not orthogonal");
            }
        }
    }

    /// Words of a defining relation, e.g. "I=ABCE=BCDF" → ["ABCE", "BCDF"].
    fn relation_words(relation: &str) -> Vec<&str> {
        let mut parts = relation.split('=');
        assert_eq!(parts.next(), Some("I"), "relation must start with I");
        parts.collect()
    }

    /// Product of the word's factor columns for one run row (letters map A→col 0, …).
    fn word_sign(row: &[f64], word: &str) -> f64 {
        word.chars()
            .map(|ch| row[(ch as usize) - ('A' as usize)])
            .product()
    }

    #[test]
    fn defining_relation_words_hold_on_design_matrix() {
        // The design matrix is its own oracle: every word W of the defining
        // relation satisfies I = W, i.e. the product of its factor columns
        // is +1 on every run. This pins generators ↔ defining_relation.
        for entry in GENERATOR_TABLE {
            let d = fractional_factorial(entry.k, entry.p).unwrap();
            for word in relation_words(entry.defining_relation) {
                for (run, row) in d.data.iter().enumerate() {
                    assert_eq!(
                        word_sign(row, word),
                        1.0,
                        "2^({}-{}): word {word} violated at run {run}",
                        entry.k,
                        entry.p
                    );
                }
            }
        }
    }

    #[test]
    fn defining_relation_word_count_and_closure() {
        // A 2^(k-p) defining relation is a group of 2^p words (incl. I):
        // exactly 2^p − 1 non-identity words, closed under symmetric-
        // difference products. Closure of every word pair is implied by the
        // matrix oracle above; here we pin the count so a truncated or
        // padded relation string cannot pass.
        for entry in GENERATOR_TABLE {
            let words = relation_words(entry.defining_relation);
            assert_eq!(
                words.len(),
                (1 << entry.p) - 1,
                "2^({}-{}): expected {} words, got {}",
                entry.k,
                entry.p,
                (1 << entry.p) - 1,
                words.len()
            );
            let distinct: std::collections::HashSet<&str> = words.iter().copied().collect();
            assert_eq!(
                distinct.len(),
                words.len(),
                "2^({}-{}): duplicate word in defining relation",
                entry.k,
                entry.p
            );
        }
    }

    #[test]
    fn resolution_matches_shortest_defining_word() {
        // Resolution = length of the shortest word in the defining relation.
        // The enum caps at V, so designs of true resolution ≥ V assert ≥ 5.
        for entry in GENERATOR_TABLE {
            let min_len = relation_words(entry.defining_relation)
                .iter()
                .map(|w| w.len())
                .min()
                .expect("non-empty relation");
            match entry.resolution {
                Resolution::III => assert_eq!(min_len, 3, "2^({}-{})", entry.k, entry.p),
                Resolution::IV => assert_eq!(min_len, 4, "2^({}-{})", entry.k, entry.p),
                Resolution::V => assert!(min_len >= 5, "2^({}-{})", entry.k, entry.p),
            }
        }
    }

    #[test]
    fn fractional_factorial_7_3_standard_generators() {
        // Pin the published 2^(7-3)_IV fraction (NIST §5.3.3.4.7, Montgomery
        // Table 8.14): E=ABC, F=BCD, G=ACD. Regression guard for upstream-018.
        let d = fractional_factorial(7, 3).unwrap();
        for row in &d.data {
            assert_eq!(row[4], row[0] * row[1] * row[2], "E != ABC");
            assert_eq!(row[5], row[1] * row[2] * row[3], "F != BCD");
            assert_eq!(row[6], row[0] * row[2] * row[3], "G != ACD");
        }
        let info = fractional_factorial_info(7, 3).unwrap();
        assert_eq!(info.resolution, Resolution::IV);
        assert_eq!(
            info.defining_relation,
            "I=ABCE=BCDF=ACDG=ADEF=BDEG=ABFG=CEFG"
        );
        assert_eq!(info.generators, vec!["E=ABC", "F=BCD", "G=ACD"]);
    }

    #[test]
    fn generator_equations_appear_in_defining_relation() {
        // Every generator "X=W" contributes the defining word sort(W ∪ {X})
        // (relation words are written in alphabetical order). Ties the derived
        // generator strings to the hand-written relation constants.
        for entry in GENERATOR_TABLE {
            let info = fractional_factorial_info(entry.k, entry.p).unwrap();
            let words = relation_words(entry.defining_relation);
            for gen in &info.generators {
                let (derived, word) = gen.split_once('=').expect("generator format X=W");
                let mut letters: Vec<char> = word.chars().chain(derived.chars()).collect();
                letters.sort_unstable();
                let expected: String = letters.into_iter().collect();
                assert!(
                    words.contains(&expected.as_str()),
                    "2^({}-{}): generator {gen} → word {expected} not in relation",
                    entry.k,
                    entry.p
                );
            }
        }
    }

    #[test]
    fn generators_hold_on_design_matrix() {
        // Each generator equation X=W must hold column-wise on the emitted
        // matrix: the derived column equals the product of the word's columns.
        for entry in GENERATOR_TABLE {
            let d = fractional_factorial(entry.k, entry.p).unwrap();
            let info = fractional_factorial_info(entry.k, entry.p).unwrap();
            for gen in &info.generators {
                let (derived, word) = gen.split_once('=').expect("generator format X=W");
                let derived_col =
                    (derived.chars().next().expect("derived letter") as usize) - ('A' as usize);
                for (run, row) in d.data.iter().enumerate() {
                    assert_eq!(
                        row[derived_col],
                        word_sign(row, word),
                        "2^({}-{}): generator {gen} violated at run {run}",
                        entry.k,
                        entry.p
                    );
                }
            }
        }
    }
}
