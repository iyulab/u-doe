//! Lenth's method for judging effects in unreplicated two-level designs.
//!
//! An unreplicated or saturated design leaves no residual degrees of freedom,
//! so effects cannot be tested against an error mean square. Lenth's method
//! estimates the standard error of an effect from the effects themselves,
//! assuming most are inactive (effect sparsity):
//!
//! ```text
//! s0  = 1.5 · median |c_j|
//! PSE = 1.5 · median { |c_j| : |c_j| < 2.5 · s0 }
//! ME  = t(0.975, m / 3) · PSE
//! ```
//!
//! over the `m` contrasts. An effect whose magnitude exceeds the margin of
//! error `ME` is judged active; on a half-normal plot `ME` is the line that
//! separates them.
//!
//! `m` counts **distinct contrasts**, not model terms. In a fractional design
//! several terms share one contrast column -- in a 2^(5-2), `A`, `B:D` and
//! `C:E` are one column -- and counting each would repeat the same magnitude,
//! pulling the medians and inflating the degrees of freedom.
//!
//! Only the individual margin of error is reported. The simultaneous margin
//! (SME) is quoted with different quantiles by different sources.
//!
//! Reference: Lenth, R.V. (1989). "Quick and Easy Analysis of Unreplicated
//! Factorials". *Technometrics* 31(4), pp. 469–473.

use crate::{analysis::effects::EffectEstimate, design::DesignMatrix, error::DoeError};

/// Lenth's pseudo standard error and margin of error for a set of effects.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LenthResult {
    /// Pseudo standard error of an effect.
    pub pse: f64,
    /// Margin of error at 95%: `t(0.975, df) · pse`. Effects larger than this
    /// in magnitude are judged active.
    pub margin_of_error: f64,
    /// Degrees of freedom of the reference t distribution, `m / 3`.
    pub df: f64,
    /// Number of distinct contrasts `m` the estimate was based on.
    pub distinct_contrasts: usize,
}

/// Apply Lenth's method to effects estimated from `design`.
///
/// `effects` is usually the output of
/// [`estimate_effects`](crate::analysis::effects::estimate_effects) for the
/// same design. Terms whose contrast columns coincide (or are negatives of each
/// other) are counted once.
///
/// # Errors
///
/// [`DoeError::NotTwoLevelCoded`] if the design is not two-level coded;
/// [`DoeError::UnsupportedDesign`] if an effect names a column the design does
/// not have, or if there are fewer than 3 distinct contrasts -- the reference
/// distribution then has less than one degree of freedom.
///
/// # Examples
///
/// ```
/// use u_doe::design::factorial::full_factorial;
/// use u_doe::analysis::effects::estimate_effects;
/// use u_doe::analysis::lenth::lenth;
///
/// let design = full_factorial(3).unwrap();
/// let responses = [10.0, 20.0, 15.0, 25.0, 12.0, 22.0, 18.0, 30.0];
/// let effects = estimate_effects(&design, &responses, 3).unwrap();
/// let result = lenth(&design, &effects).unwrap();
/// assert_eq!(result.distinct_contrasts, 7);
/// ```
pub fn lenth(design: &DesignMatrix, effects: &[EffectEstimate]) -> Result<LenthResult, DoeError> {
    if let Some((run, factor, value)) = design.two_level_violation() {
        return Err(DoeError::NotTwoLevelCoded { run, factor, value });
    }
    let n = design.run_count();
    let k = design.factor_count();
    if let Some(e) = effects.iter().find(|e| e.columns.iter().any(|&c| c >= k)) {
        return Err(DoeError::UnsupportedDesign(format!(
            "effect {} names a column the design does not have ({k} factors)",
            e.name
        )));
    }

    // One magnitude per distinct contrast column. Aliased terms share a column
    // (or its negative) and therefore the same magnitude.
    let mut columns: Vec<Vec<f64>> = Vec::new();
    let mut magnitudes: Vec<f64> = Vec::new();
    for e in effects {
        let contrast: Vec<f64> = (0..n)
            .map(|run| e.columns.iter().map(|&c| design.get(run, c)).product())
            .collect();
        let seen = columns.iter().any(|other| {
            let dot: f64 = other.iter().zip(&contrast).map(|(a, b)| a * b).sum();
            (dot.abs() - n as f64).abs() < 0.5
        });
        if !seen {
            columns.push(contrast);
            magnitudes.push(e.estimate.abs());
        }
    }

    let m = magnitudes.len();
    if m < 3 {
        return Err(DoeError::UnsupportedDesign(format!(
            "Lenth's method needs at least 3 distinct contrasts, got {m}"
        )));
    }

    let s0 = 1.5 * median(&magnitudes);
    let trimmed: Vec<f64> = magnitudes
        .iter()
        .copied()
        .filter(|&c| c < 2.5 * s0)
        .collect();
    // With most effects exactly zero, s0 is zero and nothing is below 2.5·s0:
    // there is no spread among the inactive effects to estimate from.
    let pse = if trimmed.is_empty() {
        0.0
    } else {
        1.5 * median(&trimmed)
    };
    let df = m as f64 / 3.0;
    let t = u_numflow::special::t_distribution_quantile(0.975, df);

    Ok(LenthResult {
        pse,
        margin_of_error: t * pse,
        df,
        distinct_contrasts: m,
    })
}

/// Median of a non-empty slice.
fn median(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("effect magnitudes are finite"));
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[mid]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::effects::estimate_effects;
    use crate::design::factorial::{fractional_factorial, full_factorial};

    /// Montgomery, *Design and Analysis of Experiments*, Example 6.2: the
    /// filtration-rate 2^4, in standard (Yates) order. Montgomery reports
    /// PSE = 2.625 and ME = 6.75 for this experiment.
    #[test]
    fn reproduces_the_published_filtration_rate_example() {
        let design = full_factorial(4).expect("2^4");
        let y = [
            45.0, 71.0, 48.0, 65.0, 68.0, 60.0, 80.0, 65.0, 43.0, 100.0, 45.0, 104.0, 75.0, 86.0,
            70.0, 96.0,
        ];
        let effects = estimate_effects(&design, &y, 4).expect("all 15 effects");
        assert_eq!(effects.len(), 15);
        let a = effects.iter().find(|e| e.name == "A").expect("A");
        assert!((a.estimate - 21.625).abs() < 1e-9, "A = {}", a.estimate);

        let r = lenth(&design, &effects).expect("valid");
        assert_eq!(r.distinct_contrasts, 15);
        assert!((r.pse - 2.625).abs() < 1e-9, "PSE = {}", r.pse);
        assert!((r.df - 5.0).abs() < 1e-12);
        // t(0.975, 5) = 2.5706
        assert!(
            (r.margin_of_error - 6.748).abs() < 0.01,
            "ME = {}",
            r.margin_of_error
        );
    }

    /// In a 2^(5-2) every two-factor interaction shares a column with a main
    /// effect or another interaction: 8 runs carry only 7 distinct contrasts,
    /// however many terms are estimated.
    #[test]
    fn aliased_terms_count_once() {
        let design = fractional_factorial(5, 2).expect("2^(5-2)");
        let y = [12.0, 15.0, 11.0, 18.0, 14.0, 13.0, 17.0, 16.0];
        let effects = estimate_effects(&design, &y, 2).expect("mains and 2FIs");
        assert!(effects.len() > 7, "the model terms outnumber the contrasts");
        let r = lenth(&design, &effects).expect("valid");
        assert_eq!(r.distinct_contrasts, 7);
    }

    #[test]
    fn fewer_than_three_contrasts_is_an_error() {
        let design = full_factorial(2).expect("2^2");
        let effects = estimate_effects(&design, &[1.0, 2.0, 3.0, 5.0], 1).expect("A, B");
        assert!(matches!(
            lenth(&design, &effects),
            Err(DoeError::UnsupportedDesign(_))
        ));
    }
}
