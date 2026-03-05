//! Taguchi Signal-to-Noise (SN) ratio analysis.
//!
//! Provides the three standard Taguchi quality objectives and utilities for
//! computing per-run SN ratios and factor-level effect tables.
//!
//! ## Quality Objectives
//!
//! | Objective        | Japanese term | Formula |
//! |-----------------|---------------|---------|
//! | `LargerIsBetter` | 望大 | `−10·log₁₀(mean(1/y²))` |
//! | `SmallerIsBetter`| 望小 | `−10·log₁₀(mean(y²))` |
//! | `NominalIsBest`  | 望目 | `10·log₁₀(ȳ²/s²)` |
//!
//! Reference: Taguchi, G. (1987). *System of Experimental Design*.
//! UNIPUB/Kraus International.

use std::collections::BTreeMap;

use crate::design::DesignMatrix;
use crate::error::DoeError;

/// Taguchi quality objective — determines the SN ratio formula.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SnGoal {
    /// 望大 (larger-is-better): `SN = −10·log₁₀(mean(1/y²))`
    ///
    /// Requires all response values to be strictly positive.
    LargerIsBetter,
    /// 望小 (smaller-is-better): `SN = −10·log₁₀(mean(y²))`
    SmallerIsBetter,
    /// 望目 (nominal-is-best): `SN = 10·log₁₀(ȳ²/s²)`
    ///
    /// Requires at least 2 replicates per run to estimate variance.
    NominalIsBest,
}

/// Compute Taguchi Signal-to-Noise ratios for a set of experimental runs.
///
/// # Arguments
/// * `responses` — One `Vec<f64>` per run; each inner vec contains the replicated
///   measurements for that run.
/// * `goal` — Quality objective selecting the SN formula.
///
/// # Returns
/// A `Vec<f64>` of SN values in dB, one per run, in the same order as `responses`.
///
/// # Errors
/// * [`DoeError::InvalidSpecification`] if:
///   - `responses` is empty.
///   - Any run has zero replicates.
///   - `goal` is `LargerIsBetter` and any measurement is ≤ 0.
///   - `goal` is `NominalIsBest` and any run has fewer than 2 replicates.
///
/// # Examples
///
/// ```
/// use u_doe::analysis::taguchi_sn::{signal_to_noise, SnGoal};
///
/// // y = [10, 10] → SN = −10·log₁₀(mean(1/100)) = 20 dB
/// let sn = signal_to_noise(&[vec![10.0, 10.0]], SnGoal::LargerIsBetter).unwrap();
/// assert!((sn[0] - 20.0).abs() < 1e-9);
/// ```
pub fn signal_to_noise(responses: &[Vec<f64>], goal: SnGoal) -> Result<Vec<f64>, DoeError> {
    if responses.is_empty() {
        return Err(DoeError::InvalidSpecification(
            "responses must not be empty".into(),
        ));
    }

    let mut sn_values = Vec::with_capacity(responses.len());

    for (run_idx, run) in responses.iter().enumerate() {
        if run.is_empty() {
            return Err(DoeError::InvalidSpecification(format!(
                "run {run_idx}: no replicates provided"
            )));
        }

        let r = run.len() as f64;

        match goal {
            SnGoal::LargerIsBetter => {
                // Validate: all values must be strictly positive
                for &y in run {
                    if y <= 0.0 {
                        return Err(DoeError::InvalidSpecification(format!(
                            "run {run_idx}: LargerIsBetter requires y > 0, got {y}"
                        )));
                    }
                }
                let mean_inv_sq: f64 = run.iter().map(|&y| 1.0 / (y * y)).sum::<f64>() / r;
                sn_values.push(-10.0 * mean_inv_sq.log10());
            }

            SnGoal::SmallerIsBetter => {
                let mean_sq: f64 = run.iter().map(|&y| y * y).sum::<f64>() / r;
                sn_values.push(-10.0 * mean_sq.log10());
            }

            SnGoal::NominalIsBest => {
                if run.len() < 2 {
                    return Err(DoeError::InvalidSpecification(format!(
                        "run {run_idx}: NominalIsBest requires >= 2 replicates, got {}",
                        run.len()
                    )));
                }
                let mean = run.iter().sum::<f64>() / r;
                let var = run.iter().map(|&y| (y - mean).powi(2)).sum::<f64>() / (r - 1.0);

                let sn = if var < 1e-15 {
                    // Degenerate: zero variance — return 0.0 dB
                    0.0
                } else {
                    10.0 * (mean * mean / var).log10()
                };
                sn_values.push(sn);
            }
        }
    }

    Ok(sn_values)
}

/// Mean SN per level for a single factor, used to identify the optimal setting.
#[derive(Debug, Clone)]
pub struct SnFactorEffect {
    /// Zero-based index of this factor in the design matrix.
    pub factor_index: usize,
    /// Mean SN value for each level (ordered by level key, ascending).
    pub level_sn_means: Vec<f64>,
    /// Index into `level_sn_means` for the level with the **highest** mean SN.
    pub optimal_level: usize,
    /// `max(level_sn_means) − min(level_sn_means)` — sensitivity measure.
    pub delta: f64,
}

/// Compute factor-level SN effects from a Taguchi design and its SN values.
///
/// For each factor, runs are grouped by their coded level value and the mean SN
/// within each group is computed.  The result is used to select the optimal level
/// and to rank factors by their `delta` (range of means).
///
/// # Arguments
/// * `design` — Taguchi (or any) design matrix produced by [`taguchi_array`].
/// * `sn_values` — SN ratios computed by [`signal_to_noise`], one per run.
///
/// # Returns
/// One [`SnFactorEffect`] per factor, in column order.
///
/// # Errors
/// [`DoeError::InsufficientResponses`] if `sn_values.len() != design.run_count()`.
///
/// [`taguchi_array`]: crate::design::taguchi::taguchi_array
///
/// # Examples
///
/// ```
/// use u_doe::design::taguchi::taguchi_array;
/// use u_doe::analysis::taguchi_sn::{signal_to_noise, sn_factor_effects, SnGoal};
///
/// let design = taguchi_array("L9", 4).unwrap();
/// let responses: Vec<Vec<f64>> = (0..9).map(|i| vec![10.0 + i as f64]).collect();
/// let sn = signal_to_noise(&responses, SnGoal::LargerIsBetter).unwrap();
/// let effects = sn_factor_effects(&design, &sn).unwrap();
/// assert_eq!(effects.len(), 4);
/// ```
pub fn sn_factor_effects(
    design: &DesignMatrix,
    sn_values: &[f64],
) -> Result<Vec<SnFactorEffect>, DoeError> {
    let n_runs = design.run_count();
    if sn_values.len() != n_runs {
        return Err(DoeError::InsufficientResponses {
            expected: n_runs,
            got: sn_values.len(),
        });
    }

    let n_factors = design.factor_count();
    let mut effects = Vec::with_capacity(n_factors);

    for factor in 0..n_factors {
        // Group runs by discretised level value.
        // key = round(coded_value * 10) as i64:
        //   -1.0 → -10,  0.0 → 0,  +1.0 → +10
        let mut groups: BTreeMap<i64, Vec<f64>> = BTreeMap::new();
        for (run, &sn) in sn_values.iter().enumerate() {
            let coded = design.get(run, factor);
            let key = (coded * 10.0).round() as i64;
            groups.entry(key).or_default().push(sn);
        }

        // Compute mean SN per level (BTreeMap iteration is key-sorted).
        let level_sn_means: Vec<f64> = groups
            .values()
            .map(|vals| vals.iter().sum::<f64>() / vals.len() as f64)
            .collect();

        let max_sn = level_sn_means
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let min_sn = level_sn_means
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);

        let optimal_level = level_sn_means
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).expect("SN values must be finite"))
            .map(|(i, _)| i)
            .expect("level_sn_means must be non-empty");

        effects.push(SnFactorEffect {
            factor_index: factor,
            level_sn_means,
            optimal_level,
            delta: max_sn - min_sn,
        });
    }

    Ok(effects)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::design::taguchi::taguchi_array;

    #[test]
    fn lb_higher_response_higher_sn() {
        let r = vec![vec![10.0, 10.0], vec![100.0, 100.0]];
        let sn = signal_to_noise(&r, SnGoal::LargerIsBetter).unwrap();
        assert!(sn[1] > sn[0]);
    }

    #[test]
    fn sb_lower_response_higher_sn() {
        let r = vec![vec![100.0, 100.0], vec![10.0, 10.0]];
        let sn = signal_to_noise(&r, SnGoal::SmallerIsBetter).unwrap();
        assert!(sn[1] > sn[0]);
    }

    #[test]
    fn nb_low_variance_high_sn() {
        let r1 = vec![100.0, 100.1, 99.9]; // high mean, low var
        let r2 = vec![100.0, 90.0, 110.0]; // high mean, high var
        let sn = signal_to_noise(&[r1, r2], SnGoal::NominalIsBest).unwrap();
        assert!(sn[0] > sn[1]);
    }

    #[test]
    fn lb_known_value_20db() {
        // y=[10,10]: SN = -10*log10(mean(1/100)) = -10*log10(0.01) = 20 dB
        let r = vec![vec![10.0, 10.0]];
        let sn = signal_to_noise(&r, SnGoal::LargerIsBetter).unwrap();
        assert!((sn[0] - 20.0).abs() < 1e-9, "SN={}", sn[0]);
    }

    #[test]
    fn sb_known_value_minus20db() {
        // y=[10,10]: SN = -10*log10(100) = -20 dB
        let r = vec![vec![10.0, 10.0]];
        let sn = signal_to_noise(&r, SnGoal::SmallerIsBetter).unwrap();
        assert!((sn[0] - (-20.0)).abs() < 1e-9, "SN={}", sn[0]);
    }

    #[test]
    fn lb_non_positive_error() {
        let r = vec![vec![10.0, -1.0]];
        assert!(signal_to_noise(&r, SnGoal::LargerIsBetter).is_err());
    }

    #[test]
    fn nb_insufficient_replicates_error() {
        let r = vec![vec![10.0]]; // need >= 2 for NB
        assert!(signal_to_noise(&r, SnGoal::NominalIsBest).is_err());
    }

    #[test]
    fn factor_effects_l9_four_factors() {
        let design = taguchi_array("L9", 4).unwrap();
        let responses: Vec<Vec<f64>> = (0..9)
            .map(|i| vec![10.0 + i as f64 * 2.0, 11.0 + i as f64 * 2.0])
            .collect();
        let sn = signal_to_noise(&responses, SnGoal::LargerIsBetter).unwrap();
        let effects = sn_factor_effects(&design, &sn).unwrap();
        assert_eq!(effects.len(), 4);
        for eff in &effects {
            assert_eq!(eff.level_sn_means.len(), 3); // L9 is 3-level
            assert!(eff.delta >= 0.0);
        }
    }

    #[test]
    fn sb_zero_response_yields_neg_infinity() {
        // y=[0.0]: mean(y²) = 0 → log10(0) = -inf → SN = +inf
        // This is a legitimate SB result for y=0 (perfectly small).
        let r = vec![vec![0.0]];
        let sn = signal_to_noise(&r, SnGoal::SmallerIsBetter).unwrap();
        assert!(sn[0].is_infinite() && sn[0] > 0.0);
    }

    #[test]
    fn nb_zero_variance_returns_zero() {
        // All identical replicates → var = 0 → degenerate → SN = 0.0
        let r = vec![vec![5.0, 5.0, 5.0]];
        let sn = signal_to_noise(&r, SnGoal::NominalIsBest).unwrap();
        assert_eq!(sn[0], 0.0);
    }

    #[test]
    fn factor_effects_mismatch_error() {
        let design = taguchi_array("L9", 4).unwrap();
        // Supply 5 SN values but design has 9 runs
        let sn = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!(sn_factor_effects(&design, &sn).is_err());
    }

    #[test]
    fn optimal_level_is_argmax() {
        let design = taguchi_array("L9", 1).unwrap();
        // For a single factor with 3 levels, assign ascending SN values to later runs.
        // L9 col-0: runs 0..2 = level -1, runs 3..5 = level 0, runs 6..8 = level +1
        // Assign SN: [1,1,1, 2,2,2, 3,3,3] → level +1 (index 2) should be optimal.
        let sn: Vec<f64> = (0..9).map(|i| (i / 3 + 1) as f64).collect();
        let effects = sn_factor_effects(&design, &sn).unwrap();
        assert_eq!(effects[0].optimal_level, 2);
    }
}
