//! Main effects and interaction estimation for two-level factorial designs.
//!
//! # Algorithm
//!
//! For a 2^k design with n runs, the effect of any term is estimated as:
//! `effect = (2/n) × Σᵢ contrast_i × yᵢ`
//!
//! where `contrast_i` is the product of the factor coded values (±1) for
//! that term in run i. For interactions, the contrast column is the
//! element-wise product of the involved factor columns.
//!
//! Sum of squares: `SS = n × effect² / 4`
//!
//! Reference: Montgomery, D.C. (2019). *Introduction to Statistical Quality
//! Control*, 8th ed., Section 6.3. Wiley.
//!
//! # Half-Normal Plot
//!
//! Daniel's half-normal probability plot identifies active effects in
//! unreplicated designs without a formal error estimate.
//!
//! Reference: Daniel, C. (1959). "Use of Half-Normal Plots in Interpreting
//! Factorial Two-Level Experiments". *Technometrics* 1(4), pp. 311–341.

use crate::{design::DesignMatrix, error::DoeError};

/// Estimated effect for a factorial design term.
#[derive(Debug, Clone)]
pub struct EffectEstimate {
    /// Term name: "A", "B", "AB", "ACD", etc.
    pub name: String,
    /// Effect estimate (in response units).
    pub estimate: f64,
    /// Sum of squares attributed to this effect.
    pub sum_of_squares: f64,
    /// Percentage contribution to total sum of squares.
    pub percent_contribution: f64,
}

/// Estimate main effects and interactions for a 2-level factorial design.
///
/// # Arguments
/// * `design`    — The experimental design matrix (coded ±1 values)
/// * `responses` — Observed response values, one per run
/// * `max_order` — Maximum interaction order: 1 = main effects, 2 = + 2FI
///
/// # Errors
///
/// Returns `Err` if `responses.len() != design.run_count()`.
///
/// # Examples
///
/// ```
/// use u_doe::design::factorial::full_factorial;
/// use u_doe::analysis::effects::estimate_effects;
///
/// let design = full_factorial(2).unwrap();
/// let responses = vec![28.0, 36.0, 18.0, 31.0];
/// let effects = estimate_effects(&design, &responses, 2).unwrap();
/// assert_eq!(effects.len(), 3); // A, B, AB
/// ```
pub fn estimate_effects(
    design: &DesignMatrix,
    responses: &[f64],
    max_order: usize,
) -> Result<Vec<EffectEstimate>, DoeError> {
    let n = design.run_count();
    let k = design.factor_count();

    if responses.len() != n {
        return Err(DoeError::InsufficientResponses {
            expected: n,
            got: responses.len(),
        });
    }

    // Collect terms to estimate
    let terms = build_terms(k, max_order);

    // Grand mean for SS_total
    let grand_mean = responses.iter().sum::<f64>() / n as f64;
    let ss_total: f64 = responses.iter().map(|&y| (y - grand_mean).powi(2)).sum();

    let mut estimates = Vec::with_capacity(terms.len());

    for term in &terms {
        // Contrast column = product of factor columns for this term
        let contrast: Vec<f64> = (0..n)
            .map(|run| term.iter().map(|&col| design.get(run, col)).product())
            .collect();

        // Effect = (2/n) * dot(contrast, responses)
        let dot: f64 = contrast.iter().zip(responses.iter()).map(|(&c, &y)| c * y).sum();
        let effect = (2.0 / n as f64) * dot;

        // SS = n * effect² / 4
        let ss = (n as f64) * effect * effect / 4.0;

        let pct = if ss_total > 1e-12 { ss / ss_total * 100.0 } else { 0.0 };

        estimates.push(EffectEstimate {
            name: term_name(term, &design.factor_names),
            estimate: effect,
            sum_of_squares: ss,
            percent_contribution: pct,
        });
    }

    Ok(estimates)
}

/// Generate (|effect|, half-normal quantile) data for a Half-Normal Plot.
///
/// Sorts effects by absolute value and assigns half-normal quantiles using
/// the formula: `Φ⁻¹((i + 0.5) / m)` for the i-th sorted absolute effect.
///
/// Reference: Daniel (1959), Technometrics 1(4).
pub fn half_normal_plot_data(effects: &[EffectEstimate]) -> Vec<(f64, f64)> {
    let m = effects.len();
    if m == 0 {
        return Vec::new();
    }

    let mut abs_effects: Vec<f64> = effects.iter().map(|e| e.estimate.abs()).collect();
    abs_effects.sort_by(|a, b| a.partial_cmp(b).expect("effect estimate must be finite"));

    abs_effects
        .iter()
        .enumerate()
        .map(|(i, &abs_e)| {
            // Half-normal plotting position: Φ⁻¹((1 + p_i) / 2) where p_i = (i + 0.5) / m.
            // This maps the uniform plotting position to the upper half of the standard
            // normal, giving a non-negative quantile for the half-normal distribution.
            let p = (i as f64 + 0.5) / m as f64;
            let quantile = normal_quantile((1.0 + p) / 2.0);
            (abs_e, quantile)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Build list of terms (as column index sets) up to `max_order`.
fn build_terms(k: usize, max_order: usize) -> Vec<Vec<usize>> {
    let mut terms = Vec::new();

    // Main effects: single-column terms
    for i in 0..k {
        terms.push(vec![i]);
    }

    // Two-factor interactions
    if max_order >= 2 {
        for i in 0..k {
            for j in (i + 1)..k {
                terms.push(vec![i, j]);
            }
        }
    }

    // Three-factor interactions (if requested)
    if max_order >= 3 {
        for i in 0..k {
            for j in (i + 1)..k {
                for l in (j + 1)..k {
                    terms.push(vec![i, j, l]);
                }
            }
        }
    }

    terms
}

/// Build term name from column indices and factor names.
fn term_name(cols: &[usize], factor_names: &[String]) -> String {
    cols.iter()
        .map(|&c| factor_names[c].as_str())
        .collect::<Vec<_>>()
        .join("")
}

/// Approximate inverse normal CDF (Acklam rational approximation).
///
/// Accurate to ~1e-9 for p in (0, 1).
/// Reference: Abramowitz & Stegun 26.2.17; P.J. Acklam's refinement.
fn normal_quantile(p: f64) -> f64 {
    const A: [f64; 6] = [
        -3.969_683_028_665_376e1,
         2.209_460_984_245_205e2,
        -2.759_285_104_469_687e2,
         1.383_577_518_672_69e2,
        -3.066_479_806_614_716e1,
         2.506_628_277_459_239,
    ];
    const B: [f64; 5] = [
        -5.447_609_879_822_406e1,
         1.615_858_368_580_409e2,
        -1.556_989_798_598_866e2,
         6.680_131_188_771_972e1,
        -1.328_068_155_288_572e1,
    ];
    const C: [f64; 6] = [
        -7.784_894_002_430_293e-3,
        -3.223_964_580_411_365e-1,
        -2.400_758_277_161_838,
        -2.549_732_539_343_734,
         4.374_664_141_464_968,
         2.938_163_982_698_783,
    ];
    const D: [f64; 4] = [
         7.784_695_709_041_462e-3,
         3.224_671_290_700_398e-1,
         2.445_134_137_142_996,
         3.754_408_661_907_416,
    ];

    let p_low = 0.02425;
    let p_high = 1.0 - p_low;

    if p < p_low {
        let q = (-2.0 * p.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= p_high {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::design::factorial::full_factorial;

    // Montgomery (2019) Example 6.2 — Filtration rate experiment (2^4 design)
    // p.235, Table 6.8
    fn montgomery_ex_6_2_responses() -> Vec<f64> {
        vec![
            45.0, 71.0, 48.0, 65.0, 68.0, 60.0, 80.0, 65.0,
            43.0, 100.0, 45.0, 104.0, 75.0, 86.0, 70.0, 96.0,
        ]
    }

    #[test]
    fn montgomery_main_effect_a() {
        let design = full_factorial(4).unwrap();
        let responses = montgomery_ex_6_2_responses();
        let effects = estimate_effects(&design, &responses, 1).unwrap();
        let a = effects.iter().find(|e| e.name == "A").unwrap();
        assert!((a.estimate - 21.625).abs() < 0.01, "A={}", a.estimate);
    }

    #[test]
    fn montgomery_main_effect_c() {
        let design = full_factorial(4).unwrap();
        let responses = montgomery_ex_6_2_responses();
        let effects = estimate_effects(&design, &responses, 1).unwrap();
        let c = effects.iter().find(|e| e.name == "C").unwrap();
        assert!((c.estimate - 9.875).abs() < 0.01, "C={}", c.estimate);
    }

    #[test]
    fn montgomery_main_effect_d() {
        let design = full_factorial(4).unwrap();
        let responses = montgomery_ex_6_2_responses();
        let effects = estimate_effects(&design, &responses, 1).unwrap();
        let d = effects.iter().find(|e| e.name == "D").unwrap();
        assert!((d.estimate - 14.625).abs() < 0.01, "D={}", d.estimate);
    }

    #[test]
    fn montgomery_interaction_ac() {
        let design = full_factorial(4).unwrap();
        let responses = montgomery_ex_6_2_responses();
        let effects = estimate_effects(&design, &responses, 2).unwrap();
        let ac = effects.iter().find(|e| e.name == "AC").unwrap();
        // Montgomery Table 6.8: AC = -18.125
        assert!((ac.estimate - (-18.125)).abs() < 0.01, "AC={}", ac.estimate);
    }

    #[test]
    fn effects_count_main_only() {
        let design = full_factorial(4).unwrap();
        let responses = montgomery_ex_6_2_responses();
        let effects = estimate_effects(&design, &responses, 1).unwrap();
        // max_order=1: only 4 main effects
        assert_eq!(effects.len(), 4);
    }

    #[test]
    fn effects_count_with_2fi() {
        let design = full_factorial(4).unwrap();
        let responses = montgomery_ex_6_2_responses();
        let effects = estimate_effects(&design, &responses, 2).unwrap();
        // max_order=2: 4 main + C(4,2)=6 two-factor = 10 effects
        assert_eq!(effects.len(), 10);
    }

    #[test]
    fn ss_percent_sum_near_100() {
        let design = full_factorial(4).unwrap();
        let responses = montgomery_ex_6_2_responses();
        let effects = estimate_effects(&design, &responses, 2).unwrap();
        let total_pct: f64 = effects.iter().map(|e| e.percent_contribution).sum();
        // For unreplicated 2^4 with max_order=2 (10 of 15 effects), the 3FI and 4FI
        // account for ~2.2% of total SS. Tolerance of 5.0 accommodates this.
        assert!((total_pct - 100.0).abs() < 5.0, "sum%={total_pct}");
    }

    #[test]
    fn insufficient_responses_error() {
        let design = full_factorial(2).unwrap(); // 4 runs
        let responses = vec![1.0, 2.0]; // only 2
        assert!(estimate_effects(&design, &responses, 1).is_err());
    }

    #[test]
    fn half_normal_plot_sorted() {
        let effects = vec![
            EffectEstimate { name: "A".into(), estimate: -5.0, sum_of_squares: 25.0, percent_contribution: 50.0 },
            EffectEstimate { name: "B".into(), estimate: 2.0, sum_of_squares: 4.0, percent_contribution: 8.0 },
            EffectEstimate { name: "C".into(), estimate: -1.0, sum_of_squares: 1.0, percent_contribution: 2.0 },
        ];
        let data = half_normal_plot_data(&effects);
        assert_eq!(data.len(), 3);
        // Sorted by |effect| ascending
        assert!(data[0].0 <= data[1].0);
        assert!(data[1].0 <= data[2].0);
        // Quantiles are positive (half-normal)
        for (_, q) in &data {
            assert!(*q >= 0.0, "quantile should be non-negative");
        }
    }
}
