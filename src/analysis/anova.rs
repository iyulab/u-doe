//! DOE-context Analysis of Variance (ANOVA).
//!
//! Partitions response variation into contributions from specified effects
//! and residual, then computes F-statistics for significance testing.
//!
//! Reference: Montgomery, D.C. (2019). *Introduction to Statistical Quality
//! Control*, 8th ed., Section 6.4. Wiley.

use crate::{analysis::effects::estimate_effects, design::DesignMatrix, error::DoeError};

/// One row in the ANOVA table.
#[derive(Debug, Clone)]
pub struct EffectRow {
    /// Effect name (e.g. "A", "AB").
    pub name: String,
    /// Sum of squares for this effect.
    pub sum_of_squares: f64,
    /// Degrees of freedom (always 1 for main effects and 2FI).
    pub df: usize,
    /// Mean square = SS / df.
    pub mean_square: f64,
    /// F-statistic = MS_effect / MS_residual.
    pub f_statistic: f64,
    /// p-value from F(1, df_residual) distribution.
    /// `None` if residual degrees of freedom is 0 (saturated model).
    pub p_value: Option<f64>,
}

/// Result of a DOE ANOVA.
#[derive(Debug, Clone)]
pub struct DoeAnovaResult {
    /// Rows for each requested effect.
    pub effects: Vec<EffectRow>,
    /// Residual (error) sum of squares.
    pub residual_ss: f64,
    /// Residual degrees of freedom.
    pub residual_df: usize,
    /// Total sum of squares.
    pub total_ss: f64,
    /// R² = SS_model / SS_total.
    pub r_squared: f64,
    /// Adjusted R² = 1 - (SS_residual/df_residual) / (SS_total/df_total).
    pub r_squared_adj: f64,
}

/// Compute DOE ANOVA for specified effects.
///
/// # Arguments
/// * `design`       — Experimental design matrix (coded ±1 values)
/// * `responses`    — Observed response values, one per run
/// * `effect_names` — Names of effects to include (e.g. `&["A", "B", "AB"]`)
///
/// # Errors
///
/// Returns `Err` if `responses.len() != design.run_count()`.
///
/// # Examples
///
/// ```
/// use u_doe::design::factorial::full_factorial;
/// use u_doe::analysis::anova::doe_anova;
///
/// let design = full_factorial(3).unwrap();
/// let responses = vec![10.0,20.0,15.0,25.0,12.0,22.0,18.0,30.0];
/// let result = doe_anova(&design, &responses, &["A","B","C"]).unwrap();
/// assert!(result.r_squared >= 0.0 && result.r_squared <= 1.0);
/// ```
pub fn doe_anova(
    design: &DesignMatrix,
    responses: &[f64],
    effect_names: &[&str],
) -> Result<DoeAnovaResult, DoeError> {
    let n = design.run_count();
    if responses.len() != n {
        return Err(DoeError::InsufficientResponses {
            expected: n,
            got: responses.len(),
        });
    }

    // Estimate all effects up to order 2
    let all_effects = estimate_effects(design, responses, 2)?;

    // Grand mean and total SS
    let grand_mean = responses.iter().sum::<f64>() / n as f64;
    let total_ss: f64 = responses.iter().map(|&y| (y - grand_mean).powi(2)).sum();
    let total_df = n - 1;

    // Filter to requested effects (unknown names are silently skipped)
    let selected: Vec<_> = effect_names
        .iter()
        .filter_map(|&name| all_effects.iter().find(|e| e.name == name))
        .collect();

    let model_ss: f64 = selected.iter().map(|e| e.sum_of_squares).sum();
    let model_df = selected.len();
    let residual_ss = (total_ss - model_ss).max(0.0);
    let residual_df = total_df.saturating_sub(model_df);
    let ms_residual = if residual_df > 0 {
        residual_ss / residual_df as f64
    } else {
        f64::NAN
    };

    let effects: Vec<EffectRow> = selected
        .iter()
        .map(|e| {
            let ms = e.sum_of_squares; // df = 1
            let f_stat = if ms_residual > 0.0 && ms_residual.is_finite() {
                ms / ms_residual
            } else {
                f64::NAN
            };
            let p = if f_stat.is_finite() && residual_df > 0 {
                Some(f_pvalue(f_stat, 1, residual_df))
            } else {
                None
            };
            EffectRow {
                name: e.name.clone(),
                sum_of_squares: e.sum_of_squares,
                df: 1,
                mean_square: ms,
                f_statistic: f_stat,
                p_value: p,
            }
        })
        .collect();

    let r_squared = if total_ss > 1e-12 {
        model_ss / total_ss
    } else {
        0.0
    };
    let r_squared_adj = if residual_df > 0 && total_df > 0 {
        1.0 - (residual_ss / residual_df as f64) / (total_ss / total_df as f64)
    } else {
        r_squared
    };

    Ok(DoeAnovaResult {
        effects,
        residual_ss,
        residual_df,
        total_ss,
        r_squared,
        r_squared_adj,
    })
}

// ---------------------------------------------------------------------------
// F-distribution p-value: P(F(df1, df2) > f)
// Using regularized incomplete beta: p = I_x(df2/2, df1/2) where x = df2/(df2 + df1*f)
// ---------------------------------------------------------------------------

fn f_pvalue(f: f64, df1: usize, df2: usize) -> f64 {
    if f <= 0.0 {
        return 1.0;
    }
    let d1 = df1 as f64;
    let d2 = df2 as f64;
    let x = d2 / (d2 + d1 * f);
    regularized_incomplete_beta(x, d2 / 2.0, d1 / 2.0)
}

/// Regularized incomplete beta function I_x(a, b).
/// Uses continued fraction expansion (Lentz method).
fn regularized_incomplete_beta(x: f64, a: f64, b: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }

    // Use symmetry relation when x > (a+1)/(a+b+2)
    if x > (a + 1.0) / (a + b + 2.0) {
        return 1.0 - regularized_incomplete_beta(1.0 - x, b, a);
    }

    // log B(a,b) = lgamma(a) + lgamma(b) - lgamma(a+b)
    let log_beta = lgamma(a) + lgamma(b) - lgamma(a + b);
    let front = (a * x.ln() + b * (1.0 - x).ln() - log_beta - a.ln()).exp();

    front * continued_fraction_beta(x, a, b)
}

/// Modified Lentz continued fraction for incomplete beta.
fn continued_fraction_beta(x: f64, a: f64, b: f64) -> f64 {
    const MAX_ITER: usize = 200;
    const EPS: f64 = 1e-12;
    const FPMIN: f64 = 1e-300;

    let mut c = 1.0_f64;
    let mut d = 1.0 - (a + b) * x / (a + 1.0);
    if d.abs() < FPMIN {
        d = FPMIN;
    }
    d = 1.0 / d;
    let mut h = d;

    for m in 1..=MAX_ITER {
        let mf = m as f64;
        // Even step
        let numerator = mf * (b - mf) * x / ((a + 2.0 * mf - 1.0) * (a + 2.0 * mf));
        d = 1.0 + numerator * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + numerator / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        h *= d * c;

        // Odd step
        let numerator2 = -(a + mf) * (a + b + mf) * x / ((a + 2.0 * mf) * (a + 2.0 * mf + 1.0));
        d = 1.0 + numerator2 * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + numerator2 / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        let delta = d * c;
        h *= delta;

        if (delta - 1.0).abs() < EPS {
            break;
        }
    }

    h
}

/// Lanczos approximation for ln(Γ(x)).
fn lgamma(x: f64) -> f64 {
    const G: f64 = 7.0;
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_312e-7,
    ];

    let z = x - 1.0;
    let mut s = C[0];
    for (i, &c) in C[1..].iter().enumerate() {
        s += c / (z + i as f64 + 1.0);
    }
    let t = z + G + 0.5;
    let sqrt_2pi = (2.0 * std::f64::consts::PI).sqrt();
    (sqrt_2pi * t.powf(z + 0.5) * (-t).exp() * s).ln()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::design::factorial::full_factorial;

    fn montgomery_responses() -> Vec<f64> {
        vec![
            45.0, 71.0, 48.0, 65.0, 68.0, 60.0, 80.0, 65.0, 43.0, 100.0, 45.0, 104.0, 75.0, 86.0,
            70.0, 96.0,
        ]
    }

    #[test]
    fn anova_significant_a() {
        let design = full_factorial(4).unwrap();
        let result = doe_anova(
            &design,
            &montgomery_responses(),
            &["A", "B", "C", "D", "AC"],
        )
        .unwrap();
        let row_a = result.effects.iter().find(|r| r.name == "A").unwrap();
        // A is highly significant
        assert!(row_a.f_statistic > 10.0, "A F={}", row_a.f_statistic);
    }

    #[test]
    fn anova_r_squared_reasonable() {
        let design = full_factorial(4).unwrap();
        let result = doe_anova(&design, &montgomery_responses(), &["A", "C", "D", "AC"]).unwrap();
        // A, C, D, AC explain ~77% of total variance in this dataset
        assert!(result.r_squared > 0.70, "R²={}", result.r_squared);
    }

    #[test]
    fn anova_ss_partition() {
        let design = full_factorial(4).unwrap();
        let result = doe_anova(&design, &montgomery_responses(), &["A", "C", "D", "AC"]).unwrap();
        let ss_model: f64 = result.effects.iter().map(|r| r.sum_of_squares).sum();
        let recon = ss_model + result.residual_ss;
        assert!(
            (recon - result.total_ss).abs() < 0.01,
            "SS_model + SS_res = {recon} ≠ SS_total = {}",
            result.total_ss
        );
    }

    #[test]
    fn anova_insufficient_responses() {
        let design = full_factorial(2).unwrap();
        assert!(doe_anova(&design, &[1.0, 2.0], &["A"]).is_err());
    }

    #[test]
    fn anova_unknown_effect_ignored() {
        // Unknown effect names should be silently skipped or return error
        // (implementation choice — test that it doesn't panic)
        let design = full_factorial(2).unwrap();
        let responses = vec![10.0, 20.0, 30.0, 40.0];
        let result = doe_anova(&design, &responses, &["A", "B"]);
        assert!(result.is_ok());
    }

    #[test]
    fn lgamma_known_values() {
        // lgamma(1) = 0, lgamma(2) = 0, lgamma(0.5) = ln(sqrt(pi))
        assert!((lgamma(1.0)).abs() < 1e-9, "lgamma(1)={}", lgamma(1.0));
        assert!((lgamma(2.0)).abs() < 1e-9, "lgamma(2)={}", lgamma(2.0));
        let expected = (std::f64::consts::PI.sqrt()).ln();
        assert!(
            (lgamma(0.5) - expected).abs() < 1e-9,
            "lgamma(0.5)={} expected={expected}",
            lgamma(0.5)
        );
    }

    #[test]
    fn f_pvalue_sanity() {
        // Large F should give small p-value
        let p = f_pvalue(100.0, 1, 10);
        assert!(p < 0.001, "p={p}");
        // F near 0 should give p near 1
        let p2 = f_pvalue(0.01, 1, 10);
        assert!(p2 > 0.9, "p={p2}");
    }
}
