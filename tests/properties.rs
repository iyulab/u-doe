//! Property tests: the algebra every two-level design and its analysis must
//! satisfy on any input, not only on the textbook examples the unit tests use.
//!
//! - A two-level design generator returns ±1 columns that each sum to zero and
//!   are pairwise orthogonal, with the run count the design promises.
//! - On a full factorial, `estimate_effects` agrees with least squares: each
//!   effect is twice the regression coefficient of its contrast, and the
//!   saturated model's sums of squares add up to the total.
//! - `doe_anova` decomposes the total sum of squares exactly and fits what it
//!   says it fits.

use proptest::prelude::*;
use u_doe::analysis::anova::doe_anova;
use u_doe::analysis::effects::estimate_effects;
use u_doe::design::factorial::{fractional_factorial, full_factorial};
use u_doe::design::plackett_burman::plackett_burman;
use u_doe::design::DesignMatrix;

const EPS: f64 = 1e-9;

fn column(design: &DesignMatrix, j: usize) -> Vec<f64> {
    (0..design.run_count()).map(|i| design.get(i, j)).collect()
}

/// Every entry ±1, every column balanced, every pair of columns orthogonal.
fn check_two_level_orthogonal(design: &DesignMatrix) -> Result<(), TestCaseError> {
    let (n, k) = (design.run_count(), design.factor_count());
    prop_assert_eq!(design.factor_names.len(), k);
    for row in &design.data {
        prop_assert_eq!(row.len(), k, "ragged design matrix");
        for &v in row {
            prop_assert!(v == 1.0 || v == -1.0, "entry {v} is not ±1");
        }
    }
    for j in 0..k {
        let sum: f64 = column(design, j).iter().sum();
        prop_assert!(sum.abs() < EPS, "column {j} sums to {sum}, not 0");
        for l in (j + 1)..k {
            let dot: f64 = column(design, j)
                .iter()
                .zip(column(design, l))
                .map(|(a, b)| a * b)
                .sum();
            prop_assert!(
                dot.abs() < EPS,
                "columns {j} and {l} are not orthogonal (dot {dot})"
            );
        }
    }
    prop_assert!(n > 0);
    Ok(())
}

/// Least-squares coefficient of a ±1 contrast on an orthogonal design:
/// `(x'x)^-1 x'y = (1/n) Σ x_i y_i`.
fn ls_coefficient(contrast: &[f64], responses: &[f64]) -> f64 {
    contrast
        .iter()
        .zip(responses)
        .map(|(x, y)| x * y)
        .sum::<f64>()
        / contrast.len() as f64
}

fn responses_for(runs: usize) -> impl Strategy<Value = Vec<f64>> {
    proptest::collection::vec(-100.0_f64..100.0, runs)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    // ── Design generators ────────────────────────────────────────────────

    #[test]
    fn full_factorial_is_a_balanced_orthogonal_array_of_distinct_runs(k in 2..=7usize) {
        let design = full_factorial(k).unwrap();
        prop_assert_eq!(design.run_count(), 1 << k);
        prop_assert_eq!(design.factor_count(), k);
        check_two_level_orthogonal(&design)?;
        let mut rows = design.data.clone();
        rows.sort_by(|a, b| a.partial_cmp(b).unwrap());
        rows.dedup();
        prop_assert_eq!(rows.len(), 1 << k, "a full factorial repeats a run");
    }

    #[test]
    fn fractional_factorial_is_balanced_and_orthogonal_in_its_main_effects(
        k in 3..=7usize,
        p in 1..=4usize,
    ) {
        // Not every (k, p) is a catalogued fraction; those are refused, which
        // is not this property's concern.
        if let Ok(design) = fractional_factorial(k, p) {
            prop_assert_eq!(design.run_count(), 1 << (k - p));
            prop_assert_eq!(design.factor_count(), k);
            check_two_level_orthogonal(&design)?;
        }
    }

    #[test]
    fn plackett_burman_is_balanced_and_orthogonal(k in 1..=19usize) {
        let design = plackett_burman(k).unwrap();
        prop_assert_eq!(design.factor_count(), k);
        prop_assert!(design.run_count() > k, "a screening design needs more runs than factors");
        prop_assert_eq!(design.run_count() % 4, 0, "Plackett-Burman run counts are multiples of 4");
        check_two_level_orthogonal(&design)?;
    }

    // ── Analysis ─────────────────────────────────────────────────────────

    #[test]
    fn effects_on_a_full_factorial_are_twice_the_least_squares_coefficients(
        (k, responses) in (2..=5usize).prop_flat_map(|k| (Just(k), responses_for(1 << k))),
    ) {
        let design = full_factorial(k).unwrap();
        let n = design.run_count() as f64;
        let effects = estimate_effects(&design, &responses, k).unwrap();
        // Saturated: every non-empty subset of factors is a term.
        prop_assert_eq!(effects.len(), (1usize << k) - 1);

        let total_ss: f64 = {
            let mean = responses.iter().sum::<f64>() / n;
            responses.iter().map(|y| (y - mean).powi(2)).sum()
        };
        let mut ss_sum = 0.0;
        for effect in &effects {
            let contrast: Vec<f64> = (0..design.run_count())
                .map(|i| effect.columns.iter().map(|&j| design.get(i, j)).product())
                .collect();
            let expected = 2.0 * ls_coefficient(&contrast, &responses);
            prop_assert!(
                (effect.estimate - expected).abs() < 1e-7,
                "{}: effect {} but 2 × least-squares coefficient is {expected}",
                effect.name, effect.estimate
            );
            // SS of a ±1 contrast: n · effect² / 4.
            let expected_ss = n * effect.estimate * effect.estimate / 4.0;
            prop_assert!(
                (effect.sum_of_squares - expected_ss).abs() < 1e-6 * (1.0 + expected_ss),
                "{}: SS {} but n·effect²/4 is {expected_ss}",
                effect.name, effect.sum_of_squares
            );
            ss_sum += effect.sum_of_squares;
        }
        // The saturated model explains everything: the sums of squares add up
        // to the total, and the percentages to 100.
        prop_assert!(
            (ss_sum - total_ss).abs() < 1e-6 * (1.0 + total_ss),
            "effect SS sum {ss_sum} != total SS {total_ss}"
        );
        if total_ss > 1e-9 {
            let pct: f64 = effects.iter().map(|e| e.percent_contribution).sum();
            prop_assert!((pct - 100.0).abs() < 1e-6, "percent contributions sum to {pct}");
        }
    }

    #[test]
    fn anova_decomposes_the_total_sum_of_squares_exactly(
        (k, responses) in (2..=5usize).prop_flat_map(|k| (Just(k), responses_for(1 << k))),
    ) {
        let design = full_factorial(k).unwrap();
        let names: Vec<&str> = design.factor_names.iter().map(String::as_str).collect();
        let result = doe_anova(&design, &responses, &names).unwrap();

        let effect_ss: f64 = result.effects.iter().map(|e| e.sum_of_squares).sum();
        prop_assert!(
            (effect_ss + result.residual_ss - result.total_ss).abs() < 1e-6 * (1.0 + result.total_ss),
            "model SS {effect_ss} + residual SS {} != total SS {}",
            result.residual_ss, result.total_ss
        );
        prop_assert!(result.r_squared >= -EPS && result.r_squared <= 1.0 + EPS);
        prop_assert!(result.r_squared_adj <= result.r_squared + EPS);
        prop_assert_eq!(result.residual_df, design.run_count() - 1 - k);
        prop_assert_eq!(result.fitted.len(), responses.len());
        for ((fitted, residual), y) in result.fitted.iter().zip(&result.residuals).zip(&responses) {
            prop_assert!((fitted + residual - y).abs() < 1e-7, "fitted + residual != response");
        }
        // The residuals of a least-squares fit are orthogonal to every
        // fitted contrast and sum to zero.
        prop_assert!(result.residuals.iter().sum::<f64>().abs() < 1e-6);
        for j in 0..k {
            let dot: f64 = column(&design, j).iter().zip(&result.residuals).map(|(x, r)| x * r).sum();
            prop_assert!(dot.abs() < 1e-6, "residuals are not orthogonal to factor {j}");
        }
    }
}
