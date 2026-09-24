//! Statistical power analysis for two-level factorial designs.
//!
//! The power of the t test of one effect in a 2^(k-p) design, computed the
//! way Montgomery (*Design and Analysis of Experiments*, §6.6) and the
//! standard DOE packages do: from the noncentral t distribution with the
//! error degrees of freedom of the fitted model.
//!
//! The error degrees of freedom are what a normal approximation leaves out,
//! and they matter most exactly where power is asked about -- small
//! screening designs with few replicates, where the approximation is always
//! optimistic. A 2^3 run twice (16 runs, full model, 8 error df) has power
//! 0.648 to detect an effect of 2 with σ = 1.5 at α = 0.05; the normal
//! approximation says 0.760.

use u_numflow::special::{noncentral_t_cdf, t_distribution_quantile};

use crate::error::DoeError;

/// Compute the power of the two-sided t test of one effect in a 2^(k-p)
/// factorial design.
///
/// With `N = n_replicates · 2^(k-p)` runs an effect (the difference between
/// the mean response at its high and low level) has standard error
/// `SE = 2σ/√N`. Its t statistic follows a noncentral t with noncentrality
/// `δ = effect_size / SE` and `df_E` error degrees of freedom, and
///
/// ```text
/// power = P(|T| > t(1 − α/2; df_E)) = 1 − F(c; df_E, δ) + F(−c; df_E, δ)
/// ```
///
/// `df_E = N − 1 − model_terms`, where `model_terms` counts the terms of the
/// fitted model other than the mean. `None` is the full model -- every main
/// effect and interaction the fraction can estimate, `2^(k-p) − 1` terms --
/// which leaves `df_E = 2^(k-p)·(n_replicates − 1)`. Pass a smaller count for
/// a reduced model, e.g. `Some(k)` for main effects only.
///
/// # Arguments
/// * `k` — total number of factors
/// * `p` — number of generators (the fraction is 2^(-p); p = 0 is a full factorial)
/// * `n_replicates` — number of replicates of the design
/// * `effect_size` — the detectable effect δ (in response units)
/// * `sigma` — process standard deviation σ
/// * `alpha` — type-I error rate (e.g., 0.05)
/// * `model_terms` — terms in the fitted model besides the mean (`None` = full model)
///
/// # Errors
/// [`DoeError::ParameterOutOfRange`] if `k` is 0, `p >= k`, `n_replicates` is
/// 0, or `model_terms` is 0 or more than the `2^(k-p) − 1` effects the
/// fraction can estimate; [`DoeError::ValueOutOfDomain`] unless `effect_size`
/// and `sigma` are finite and positive and `alpha` lies strictly between 0 and
/// 1; [`DoeError::NoErrorDegreesOfFreedom`] if the model leaves no error
/// degrees of freedom (an unreplicated full model), where no effect can be
/// tested and so there is no power to report.
///
/// # Examples
///
/// ```
/// use u_doe::power::two_level_factorial_power;
/// // 2^3, two replicates, full model: 8 error df.
/// let power = two_level_factorial_power(3, 0, 2, 2.0, 1.5, 0.05, None).unwrap();
/// assert!((power - 0.648).abs() < 1e-3);
/// // Unreplicated, full model: nothing left to estimate the error from.
/// assert!(two_level_factorial_power(3, 0, 1, 2.0, 1.5, 0.05, None).is_err());
/// // Unreplicated, main effects only: 8 − 1 − 3 = 4 error df.
/// assert!(two_level_factorial_power(3, 0, 1, 2.0, 1.5, 0.05, Some(3)).is_ok());
/// ```
pub fn two_level_factorial_power(
    k: usize,
    p: usize,
    n_replicates: usize,
    effect_size: f64,
    sigma: f64,
    alpha: f64,
    model_terms: Option<usize>,
) -> Result<f64, DoeError> {
    check_design(k, p)?;
    if n_replicates == 0 {
        return Err(DoeError::ParameterOutOfRange {
            parameter: "n_replicates",
            min: 1,
            max: None,
            got: 0,
        });
    }
    check_effect(effect_size, sigma, alpha)?;
    let design = Design::new(k, p, model_terms)?;
    let df = design
        .error_df(n_replicates)
        .ok_or(DoeError::NoErrorDegreesOfFreedom {
            terms: design.terms,
            runs: design.runs(n_replicates),
        })?;
    Ok(power_unchecked(
        &design,
        n_replicates,
        df,
        effect_size,
        sigma,
        alpha,
    ))
}

/// Find the minimum number of replicates to achieve `target_power`.
///
/// Searches `1..=max_replicates` and returns the first replicate count whose
/// power -- computed as in [`two_level_factorial_power`] -- reaches
/// `target_power`, or `None` if none in the range does. A replicate count
/// that leaves no error degrees of freedom (one replicate of the full model)
/// has no power and is passed over.
///
/// # Arguments
/// * `k`, `p`, `effect_size`, `sigma`, `alpha`, `model_terms` — as in [`two_level_factorial_power`]
/// * `target_power` — desired power (e.g., 0.80)
/// * `max_replicates` — upper bound for the search
///
/// # Errors
/// As [`two_level_factorial_power`] (except `NoErrorDegreesOfFreedom`), and
/// [`DoeError::ParameterOutOfRange`] if `max_replicates` is 0 or
/// [`DoeError::ValueOutOfDomain`] unless `target_power` lies strictly between
/// 0 and 1.
///
/// # Examples
///
/// ```
/// use u_doe::power::required_replicates;
/// assert_eq!(required_replicates(1, 0, 1.0, 1.0, 0.05, 0.80, 50, None).unwrap(), Some(17));
/// assert_eq!(required_replicates(1, 0, 1.0, 1.0, 0.05, 0.80, 10, None).unwrap(), None);
/// ```
#[allow(clippy::too_many_arguments)]
pub fn required_replicates(
    k: usize,
    p: usize,
    effect_size: f64,
    sigma: f64,
    alpha: f64,
    target_power: f64,
    max_replicates: usize,
    model_terms: Option<usize>,
) -> Result<Option<usize>, DoeError> {
    check_design(k, p)?;
    check_max_replicates(max_replicates)?;
    check_effect(effect_size, sigma, alpha)?;
    if !(target_power > 0.0 && target_power < 1.0) {
        return Err(DoeError::ValueOutOfDomain {
            index: None,
            parameter: "target_power",
            value: target_power,
            domain: OPEN_UNIT_INTERVAL,
        });
    }
    let design = Design::new(k, p, model_terms)?;
    Ok((1..=max_replicates).find(|&n| {
        design.error_df(n).is_some_and(|df| {
            power_unchecked(&design, n, df, effect_size, sigma, alpha) >= target_power
        })
    }))
}

/// Power curve: `(replicates, power)` pairs for `1..=max_replicates`.
///
/// Useful for visualising how power grows with the number of replicates. A
/// replicate count that leaves no error degrees of freedom (one replicate of
/// the full model) has no power and is left out, so the curve then starts at
/// two replicates.
///
/// # Errors
/// As [`two_level_factorial_power`] (except `NoErrorDegreesOfFreedom`), and
/// [`DoeError::ParameterOutOfRange`] if `max_replicates` is 0.
pub fn power_curve(
    k: usize,
    p: usize,
    effect_size: f64,
    sigma: f64,
    alpha: f64,
    max_replicates: usize,
    model_terms: Option<usize>,
) -> Result<Vec<(usize, f64)>, DoeError> {
    check_design(k, p)?;
    check_max_replicates(max_replicates)?;
    check_effect(effect_size, sigma, alpha)?;
    let design = Design::new(k, p, model_terms)?;
    Ok((1..=max_replicates)
        .filter_map(|n| {
            design.error_df(n).map(|df| {
                (
                    n,
                    power_unchecked(&design, n, df, effect_size, sigma, alpha),
                )
            })
        })
        .collect())
}

const POSITIVE: &str = "finite and greater than 0";
const OPEN_UNIT_INTERVAL: &str = "strictly between 0 and 1";

/// A checked 2^(k-p) design and the number of model terms fitted to it.
struct Design {
    k: usize,
    p: usize,
    terms: usize,
}

impl Design {
    /// `model_terms`, checked against the `2^(k-p) − 1` effects the fraction
    /// can estimate; `None` is all of them (the full model). `k` and `p` must
    /// already be checked.
    fn new(k: usize, p: usize, model_terms: Option<usize>) -> Result<Self, DoeError> {
        // Saturates for a fraction too large to count, which no term count reaches.
        let estimable = u32::try_from(k - p)
            .ok()
            .and_then(|e| 1usize.checked_shl(e))
            .map_or(usize::MAX, |runs| runs - 1);
        let terms = match model_terms {
            None => estimable,
            Some(terms) if (1..=estimable).contains(&terms) => terms,
            Some(terms) => {
                return Err(DoeError::ParameterOutOfRange {
                    parameter: "model_terms",
                    min: 1,
                    max: Some(estimable),
                    got: terms,
                })
            }
        };
        Ok(Self { k, p, terms })
    }

    /// Runs per replicate, as a float: an integer shift overflows once `k - p`
    /// reaches the width of `usize`, which is 32 bits in WebAssembly.
    fn runs_per_replicate(&self) -> f64 {
        2f64.powf((self.k - self.p) as f64)
    }

    fn runs(&self, n_replicates: usize) -> usize {
        (n_replicates as f64 * self.runs_per_replicate()).min(usize::MAX as f64) as usize
    }

    /// `N − 1 − terms`, or `None` when no error degree of freedom is left.
    fn error_df(&self, n_replicates: usize) -> Option<f64> {
        let df = n_replicates as f64 * self.runs_per_replicate() - 1.0 - self.terms as f64;
        (df >= 1.0).then_some(df)
    }
}

fn check_design(k: usize, p: usize) -> Result<(), DoeError> {
    if k == 0 {
        return Err(DoeError::ParameterOutOfRange {
            parameter: "k",
            min: 1,
            max: None,
            got: 0,
        });
    }
    if p >= k {
        return Err(DoeError::ParameterOutOfRange {
            parameter: "p",
            min: 0,
            max: Some(k - 1),
            got: p,
        });
    }
    Ok(())
}

fn check_max_replicates(max_replicates: usize) -> Result<(), DoeError> {
    if max_replicates == 0 {
        return Err(DoeError::ParameterOutOfRange {
            parameter: "max_replicates",
            min: 1,
            max: None,
            got: 0,
        });
    }
    Ok(())
}

fn check_effect(effect_size: f64, sigma: f64, alpha: f64) -> Result<(), DoeError> {
    for (parameter, value) in [("effect_size", effect_size), ("sigma", sigma)] {
        if !(value.is_finite() && value > 0.0) {
            return Err(DoeError::ValueOutOfDomain {
                index: None,
                parameter,
                value,
                domain: POSITIVE,
            });
        }
    }
    if !(alpha > 0.0 && alpha < 1.0) {
        return Err(DoeError::ValueOutOfDomain {
            index: None,
            parameter: "alpha",
            value: alpha,
            domain: OPEN_UNIT_INTERVAL,
        });
    }
    Ok(())
}

/// Power for inputs already checked, at `df` (≥ 1) error degrees of freedom.
fn power_unchecked(
    design: &Design,
    n_replicates: usize,
    df: f64,
    effect_size: f64,
    sigma: f64,
    alpha: f64,
) -> f64 {
    let total_runs = n_replicates as f64 * design.runs_per_replicate();
    let se_effect = sigma * 2.0 / total_runs.sqrt();
    let ncp = effect_size / se_effect;
    // Beyond the range the noncentral t series is accurate for, the test
    // rejects with certainty to double precision.
    if ncp > 37.0 {
        return 1.0;
    }
    let crit = t_distribution_quantile(1.0 - alpha / 2.0, df);
    let power = 1.0 - noncentral_t_cdf(crit, df, ncp) + noncentral_t_cdf(-crit, df, ncp);
    power.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn power(k: usize, p: usize, n: usize, effect: f64, sigma: f64, alpha: f64) -> f64 {
        two_level_factorial_power(k, p, n, effect, sigma, alpha, None).expect("valid power inputs")
    }

    fn power_terms(k: usize, p: usize, n: usize, effect: f64, sigma: f64, terms: usize) -> f64 {
        two_level_factorial_power(k, p, n, effect, sigma, 0.05, Some(terms))
            .expect("valid power inputs")
    }

    // -----------------------------------------------------------------------
    // Reference values. Oracle: the two-sided noncentral t power
    // 1 − F(c; df, δ) + F(−c; df, δ), each F computed by composite Simpson
    // integration of Φ(t·√(u/df) − δ) over the χ²(df) density, independently
    // of the series this crate uses.
    // -----------------------------------------------------------------------

    /// 2^3 × 2 replicates, full model: 16 runs, 8 error df, δ = 2/(2·1.5/4).
    /// The normal approximation this replaced reported 0.7601.
    #[test]
    fn replicated_2_3_full_model_matches_noncentral_t() {
        let pw = power(3, 0, 2, 2.0, 1.5, 0.05);
        assert!((pw - 0.648_009_710_6).abs() < 1e-6, "power = {pw}");
    }

    #[test]
    fn replicated_2_3_unit_sigma_matches_noncentral_t() {
        let pw = power(3, 0, 2, 2.0, 1.0, 0.05);
        assert!((pw - 0.936_742_872_9).abs() < 1e-6, "power = {pw}");
    }

    /// Unreplicated 2^3 fitted with main effects only: 8 − 1 − 3 = 4 error df.
    #[test]
    fn reduced_model_uses_its_own_error_df() {
        let pw = power_terms(3, 0, 1, 2.0, 1.5, 3);
        assert!((pw - 0.305_684_951_3).abs() < 1e-6, "power = {pw}");
    }

    /// k = 1, δ = σ = 1: 2n runs, 2n − 2 error df. Power first reaches 0.80
    /// at n = 17 (0.7814 at 16, 0.8070 at 17). The normal approximation said
    /// 16 -- one replicate short, the direction it always errs in.
    #[test]
    fn required_replicates_reference_n17() {
        assert_eq!(
            required_replicates(1, 0, 1.0, 1.0, 0.05, 0.80, 50, None).unwrap(),
            Some(17)
        );
        let at16 = power(1, 0, 16, 1.0, 1.0, 0.05);
        let at17 = power(1, 0, 17, 1.0, 1.0, 0.05);
        assert!((at16 - 0.781_397_792_5).abs() < 1e-6, "n=16: {at16}");
        assert!((at17 - 0.807_036_715_1).abs() < 1e-6, "n=17: {at17}");
    }

    /// Reaching the target exactly at `max_replicates` and not reaching it at
    /// all must stay distinguishable.
    #[test]
    fn required_replicates_tells_reached_at_the_bound_from_not_reached() {
        assert_eq!(
            required_replicates(1, 0, 1.0, 1.0, 0.05, 0.80, 17, None).unwrap(),
            Some(17)
        );
        assert_eq!(
            required_replicates(1, 0, 1.0, 1.0, 0.05, 0.80, 16, None).unwrap(),
            None
        );
    }

    /// With many error df the t test approaches the z test, so the exact power
    /// approaches the normal approximation.
    #[test]
    fn large_error_df_approaches_the_normal_approximation() {
        use u_numflow::special::{inverse_normal_cdf, standard_normal_cdf};
        let n = 400;
        let exact = power(2, 0, n, 0.3, 1.0, 0.05);
        let ncp = 0.3 / (2.0 / (4.0 * n as f64).sqrt());
        let z = inverse_normal_cdf(0.975);
        let normal = standard_normal_cdf(ncp - z) + standard_normal_cdf(-ncp - z);
        assert!(
            (exact - normal).abs() < 1e-3,
            "exact={exact} normal={normal}"
        );
        assert!(exact < normal, "the t test must be the less powerful one");
    }

    // -----------------------------------------------------------------------
    // Error degrees of freedom
    // -----------------------------------------------------------------------

    /// An unreplicated full model has no error df: there is no t test, so no
    /// power. It used to be answered with the normal approximation.
    #[test]
    fn unreplicated_full_model_is_refused() {
        assert_eq!(
            two_level_factorial_power(3, 0, 1, 2.0, 1.5, 0.05, None),
            Err(DoeError::NoErrorDegreesOfFreedom { terms: 7, runs: 8 })
        );
        // The same with the full model named explicitly.
        assert_eq!(
            two_level_factorial_power(3, 0, 1, 2.0, 1.5, 0.05, Some(7)),
            Err(DoeError::NoErrorDegreesOfFreedom { terms: 7, runs: 8 })
        );
    }

    #[test]
    fn model_terms_outside_what_the_fraction_estimates_are_refused() {
        for terms in [0, 8] {
            assert_eq!(
                two_level_factorial_power(3, 0, 2, 2.0, 1.5, 0.05, Some(terms)),
                Err(DoeError::ParameterOutOfRange {
                    parameter: "model_terms",
                    min: 1,
                    max: Some(7),
                    got: terms,
                })
            );
        }
        // A 2^(5-2) estimates 7 effects, not 31.
        assert!(matches!(
            two_level_factorial_power(5, 2, 2, 1.0, 1.0, 0.05, Some(8)),
            Err(DoeError::ParameterOutOfRange {
                parameter: "model_terms",
                max: Some(7),
                ..
            })
        ));
    }

    #[test]
    fn curve_leaves_out_replicate_counts_without_error_df() {
        let full = power_curve(3, 0, 2.0, 1.0, 0.05, 8, None).unwrap();
        assert_eq!(full.len(), 7);
        assert_eq!(full[0].0, 2);
        let reduced = power_curve(3, 0, 2.0, 1.0, 0.05, 8, Some(3)).unwrap();
        assert_eq!(reduced.len(), 8);
        assert_eq!(reduced[0].0, 1);
    }

    #[test]
    fn required_replicates_passes_over_counts_without_error_df() {
        // Power at n = 2 is 0.937 (above), n = 1 has none.
        assert_eq!(
            required_replicates(3, 0, 2.0, 1.0, 0.05, 0.50, 10, None).unwrap(),
            Some(2)
        );
    }

    // -----------------------------------------------------------------------
    // Shape
    // -----------------------------------------------------------------------

    #[test]
    fn power_increases_with_replicates() {
        let p2 = power(3, 0, 2, 2.0, 1.0, 0.05);
        let p3 = power(3, 0, 3, 2.0, 1.0, 0.05);
        let p4 = power(3, 0, 4, 2.0, 1.0, 0.05);
        assert!(p3 > p2 && p4 > p3, "p2={p2} p3={p3} p4={p4}");
    }

    #[test]
    fn power_increases_with_effect_size() {
        let small = power_terms(4, 0, 1, 0.5, 1.0, 4);
        let large = power_terms(4, 0, 1, 3.0, 1.0, 4);
        assert!(large > small, "small={small} large={large}");
    }

    #[test]
    fn fewer_model_terms_leave_more_error_df_and_more_power() {
        let main = power_terms(4, 0, 1, 1.5, 1.0, 4);
        let two_fi = power_terms(4, 0, 1, 1.5, 1.0, 10);
        assert!(main > two_fi, "main={main} with 2FI={two_fi}");
    }

    #[test]
    fn power_in_range() {
        let p = power(4, 0, 2, 2.0, 1.0, 0.05);
        assert!((0.0..=1.0).contains(&p), "p={p}");
    }

    #[test]
    fn required_replicates_monotone() {
        let r80 = required_replicates(3, 0, 2.0, 1.0, 0.05, 0.80, 10, None).unwrap();
        let r90 = required_replicates(3, 0, 2.0, 1.0, 0.05, 0.90, 10, None).unwrap();
        assert!(r90 >= r80, "r80={r80:?} r90={r90:?}");
    }

    #[test]
    fn power_curve_non_decreasing() {
        let curve = power_curve(3, 0, 2.0, 1.0, 0.05, 8, None).unwrap();
        for w in curve.windows(2) {
            assert!(
                w[1].1 >= w[0].1 - 1e-9,
                "non-monotone: ({}, {}) -> ({}, {})",
                w[0].0,
                w[0].1,
                w[1].0,
                w[1].1
            );
        }
    }

    #[test]
    fn power_fractional_vs_full() {
        // Fewer runs (and so fewer error df) at the same model: less power.
        let full = power_terms(4, 0, 1, 1.5, 1.0, 4);
        let frac = power_terms(4, 1, 1, 1.5, 1.0, 4);
        assert!(full > frac, "full={full} frac={frac}");
    }

    /// `2^(k-p)` used to be an integer shift, which overflows at the width of
    /// `usize`. A design that large has power 1 for any real effect.
    #[test]
    fn a_very_large_design_does_not_overflow() {
        assert_eq!(power(70, 0, 1, 1.0, 1.0, 0.05), 1.0);
    }

    /// Inputs with no power to report used to return 0.0, the same number a
    /// valid design with a tiny effect returns.
    #[test]
    fn invalid_inputs_are_refused() {
        let range = |parameter: &'static str, min, max, got| DoeError::ParameterOutOfRange {
            parameter,
            min,
            max,
            got,
        };
        let pw = |k, p, n, e, s, a| two_level_factorial_power(k, p, n, e, s, a, None);
        assert_eq!(pw(0, 0, 1, 1.0, 1.0, 0.05), Err(range("k", 1, None, 0)));
        assert_eq!(pw(3, 3, 1, 1.0, 1.0, 0.05), Err(range("p", 0, Some(2), 3)));
        assert_eq!(
            pw(3, 0, 0, 1.0, 1.0, 0.05),
            Err(range("n_replicates", 1, None, 0))
        );
        let domain = |f: Result<f64, DoeError>| match f {
            Err(DoeError::ValueOutOfDomain { parameter, .. }) => parameter,
            other => panic!("{other:?}"),
        };
        assert_eq!(domain(pw(3, 0, 2, 0.0, 1.0, 0.05)), "effect_size");
        assert_eq!(domain(pw(3, 0, 2, 1.0, -1.0, 0.05)), "sigma");
        assert_eq!(domain(pw(3, 0, 2, f64::NAN, 1.0, 0.05)), "effect_size");
        for alpha in [0.0, 1.0, 2.0, f64::NAN] {
            assert_eq!(domain(pw(3, 0, 2, 1.0, 1.0, alpha)), "alpha");
        }
        assert!(matches!(
            required_replicates(3, 0, 1.0, 1.0, 0.05, 1.0, 10, None),
            Err(DoeError::ValueOutOfDomain {
                parameter: "target_power",
                ..
            })
        ));
        assert_eq!(
            power_curve(3, 0, 1.0, 1.0, 0.05, 0, None),
            Err(range("max_replicates", 1, None, 0))
        );
    }
}
