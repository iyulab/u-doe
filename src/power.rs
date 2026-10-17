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

use std::ops::Bound::{self, Excluded, Unbounded};

use u_numflow::special::{noncentral_t_cdf, t_distribution_quantile};

use crate::error::DoeError;

/// The terms of the model fitted to a 2^(k-p) design, besides the mean: a
/// count, or a named model whose count the design resolves.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PowerModel {
    /// This many terms (`1..=2^(k-p) − 1`).
    Terms(usize),
    /// The `k` main effects.
    MainEffects,
    /// The `k` main effects and the `k(k−1)/2` two-factor interactions, as
    /// many of them as the fraction can estimate.
    TwoFactor,
    /// Every effect the fraction can estimate, `2^(k-p) − 1` of them.
    Full,
}

/// A power calculation with the quantities it was made from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PowerAnalysis {
    /// The power of the two-sided t test of one effect, in `[0, 1]`.
    pub power: f64,
    /// The runs in all: `n_replicates · 2^(k-p)`.
    pub runs: usize,
    /// The terms of the model besides the mean, as resolved from the
    /// [`PowerModel`] -- a named model capped at what the fraction can estimate.
    pub model_terms: usize,
    /// The error degrees of freedom: `runs − 1 − model_terms`.
    pub error_df: usize,
}

/// The number of model terms `model` stands for in a 2^(k-p) design -- what
/// the power calculation fits.
///
/// A named model asks for more terms than a small fraction can estimate (a
/// resolution III fraction aliases its interactions, so `TwoFactor` is capped
/// at the `2^(k-p) − 1` effects there are); the count returned is the capped
/// one. A [`PowerModel::Terms`] count outside `1..=2^(k-p) − 1` is refused.
///
/// # Errors
/// [`DoeError::ParameterOutOfRange`] if `k` is 0, `p >= k`, or a `Terms` count
/// is outside the range.
///
/// # Examples
///
/// ```
/// use u_doe::power::{model_terms, PowerModel};
/// assert_eq!(model_terms(3, 0, PowerModel::MainEffects).unwrap(), 3);
/// assert_eq!(model_terms(3, 0, PowerModel::TwoFactor).unwrap(), 6);
/// // 2^(5-2) has 8 runs: 7 effects, though 5 + 10 terms were asked for.
/// assert_eq!(model_terms(5, 2, PowerModel::TwoFactor).unwrap(), 7);
/// ```
pub fn model_terms(k: usize, p: usize, model: PowerModel) -> Result<usize, DoeError> {
    check_design(k, p)?;
    Ok(Design::new(k, p, model)?.terms)
}

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
/// fraction can estimate, or unless `effect_size` and `sigma` are positive
/// (`exclusive_min` 0) and `alpha` lies strictly between 0 and 1
/// (`exclusive_min` 0, `exclusive_max` 1); [`DoeError::ValueNotFinite`] if one
/// of those three is NaN or infinite; [`DoeError::NoErrorDegreesOfFreedom`] if the model leaves no error
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
    two_level_factorial_analysis(
        k,
        p,
        n_replicates,
        effect_size,
        sigma,
        alpha,
        model_terms.map_or(PowerModel::Full, PowerModel::Terms),
    )
    .map(|a| a.power)
}

/// [`two_level_factorial_power`] with the quantities it is made from: the run
/// count, the model's term count and the error degrees of freedom, and the
/// model given by name as well as by count -- so a caller shows them beside
/// the power without applying the estimability cap and the `N − 1 − terms`
/// rule a second time.
///
/// # Errors
/// As [`two_level_factorial_power`].
///
/// # Examples
///
/// ```
/// use u_doe::power::{two_level_factorial_analysis, PowerModel};
/// let a = two_level_factorial_analysis(3, 0, 2, 2.0, 1.5, 0.05, PowerModel::Full).unwrap();
/// assert!((a.power - 0.648).abs() < 1e-3);
/// assert_eq!((a.runs, a.model_terms, a.error_df), (16, 7, 8));
/// // Main effects only: 16 − 1 − 3 = 12 error df.
/// let m = two_level_factorial_analysis(3, 0, 2, 2.0, 1.5, 0.05, PowerModel::MainEffects).unwrap();
/// assert_eq!((m.model_terms, m.error_df), (3, 12));
/// ```
pub fn two_level_factorial_analysis(
    k: usize,
    p: usize,
    n_replicates: usize,
    effect_size: f64,
    sigma: f64,
    alpha: f64,
    model: PowerModel,
) -> Result<PowerAnalysis, DoeError> {
    check_design(k, p)?;
    if n_replicates == 0 {
        return Err(DoeError::count_out_of_range("n_replicates", 1.., 0));
    }
    check_effect(effect_size, sigma, alpha)?;
    let design = Design::new(k, p, model)?;
    let df = design
        .error_df(n_replicates)
        .ok_or(DoeError::NoErrorDegreesOfFreedom {
            terms: design.terms,
            runs: design.runs(n_replicates),
        })?;
    Ok(PowerAnalysis {
        power: power_unchecked(&design, n_replicates, df, effect_size, sigma, alpha),
        runs: design.runs(n_replicates),
        model_terms: design.terms,
        error_df: df as usize,
    })
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
/// [`DoeError::ParameterOutOfRange`] if `max_replicates` is 0 or unless
/// `target_power` lies strictly between 0 and 1.
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
    DoeError::check_real("target_power", target_power, OPEN_UNIT_INTERVAL)?;
    let design = Design::new(
        k,
        p,
        model_terms.map_or(PowerModel::Full, PowerModel::Terms),
    )?;
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
    let design = Design::new(
        k,
        p,
        model_terms.map_or(PowerModel::Full, PowerModel::Terms),
    )?;
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

/// What an effect size or a standard deviation takes: greater than 0.
const POSITIVE: (Bound<f64>, Bound<f64>) = (Excluded(0.0), Unbounded);
/// What a probability -- `alpha`, a target power -- takes: strictly between 0
/// and 1.
const OPEN_UNIT_INTERVAL: (Bound<f64>, Bound<f64>) = (Excluded(0.0), Excluded(1.0));

/// A checked 2^(k-p) design and the number of model terms fitted to it.
struct Design {
    k: usize,
    p: usize,
    terms: usize,
}

impl Design {
    /// `model`'s terms: a count checked against the `2^(k-p) − 1` effects the
    /// fraction can estimate, or a named model capped at that number. `k` and
    /// `p` must already be checked.
    fn new(k: usize, p: usize, model: PowerModel) -> Result<Self, DoeError> {
        // Saturates for a fraction too large to count, which no term count reaches.
        let estimable = u32::try_from(k - p)
            .ok()
            .and_then(|e| 1usize.checked_shl(e))
            .map_or(usize::MAX, |runs| runs - 1);
        let named = |wanted: usize| wanted.min(estimable);
        let terms = match model {
            PowerModel::Full => estimable,
            PowerModel::MainEffects => named(k),
            PowerModel::TwoFactor => named(k.saturating_add(k.saturating_mul(k - 1) / 2)),
            PowerModel::Terms(terms) if (1..=estimable).contains(&terms) => terms,
            PowerModel::Terms(terms) => {
                return Err(DoeError::count_out_of_range(
                    "model_terms",
                    1..=estimable,
                    terms,
                ))
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
        return Err(DoeError::count_out_of_range("k", 1.., 0));
    }
    if p >= k {
        return Err(DoeError::count_out_of_range("p", 0..k, p));
    }
    Ok(())
}

fn check_max_replicates(max_replicates: usize) -> Result<(), DoeError> {
    if max_replicates == 0 {
        return Err(DoeError::count_out_of_range("max_replicates", 1.., 0));
    }
    Ok(())
}

fn check_effect(effect_size: f64, sigma: f64, alpha: f64) -> Result<(), DoeError> {
    DoeError::check_real("effect_size", effect_size, POSITIVE)?;
    DoeError::check_real("sigma", sigma, POSITIVE)?;
    DoeError::check_real("alpha", alpha, OPEN_UNIT_INTERVAL)?;
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
mod analysis_tests {
    use super::*;

    #[test]
    fn named_models_resolve_and_cap_at_what_the_fraction_estimates() {
        // 2^3: 7 effects. Main effects 3, with two-factor interactions 6.
        assert_eq!(model_terms(3, 0, PowerModel::MainEffects).unwrap(), 3);
        assert_eq!(model_terms(3, 0, PowerModel::TwoFactor).unwrap(), 6);
        assert_eq!(model_terms(3, 0, PowerModel::Full).unwrap(), 7);
        // 2^(5-2): 8 runs, 7 effects -- 5 + 10 terms asked for, 7 estimable.
        assert_eq!(model_terms(5, 2, PowerModel::TwoFactor).unwrap(), 7);
        // 2^(6-3): 8 runs, 7 effects, 6 main effects.
        assert_eq!(model_terms(6, 3, PowerModel::MainEffects).unwrap(), 6);
        assert_eq!(model_terms(6, 3, PowerModel::TwoFactor).unwrap(), 7);
        // A count is checked, not capped.
        assert!(model_terms(3, 0, PowerModel::Terms(8)).is_err());
        assert!(model_terms(3, 0, PowerModel::Terms(0)).is_err());
        assert!(model_terms(0, 0, PowerModel::Full).is_err());
    }

    #[test]
    fn the_analysis_carries_runs_terms_and_error_df() {
        for (k, p, n, model) in [
            (3, 0, 2, PowerModel::Full),
            (3, 0, 1, PowerModel::MainEffects),
            (4, 1, 3, PowerModel::TwoFactor),
            (5, 2, 4, PowerModel::TwoFactor),
        ] {
            let a = two_level_factorial_analysis(k, p, n, 2.0, 1.5, 0.05, model).unwrap();
            let runs = n * (1usize << (k - p));
            assert_eq!(a.runs, runs);
            assert_eq!(a.model_terms, model_terms(k, p, model).unwrap());
            assert_eq!(a.error_df, runs - 1 - a.model_terms);
            // The same number the bare function gives for the resolved count.
            let bare =
                two_level_factorial_power(k, p, n, 2.0, 1.5, 0.05, Some(a.model_terms)).unwrap();
            assert_eq!(a.power, bare);
        }
        // No error df left: refused with the figures, as before.
        assert!(matches!(
            two_level_factorial_analysis(3, 0, 1, 2.0, 1.5, 0.05, PowerModel::Full),
            Err(DoeError::NoErrorDegreesOfFreedom { terms: 7, runs: 8 })
        ));
    }
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
                    index: None,
                    min: Some(1.0),
                    max: Some(7.0),
                    exclusive_min: None,
                    exclusive_max: None,
                    got: terms as f64,
                })
            );
        }
        // A 2^(5-2) estimates 7 effects, not 31.
        assert!(matches!(
            two_level_factorial_power(5, 2, 2, 1.0, 1.0, 0.05, Some(8)),
            Err(DoeError::ParameterOutOfRange {
                parameter: "model_terms",
                max: Some(7.0),
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
        let range = |parameter: &'static str, min: usize, max: Option<usize>, got: usize| {
            DoeError::ParameterOutOfRange {
                parameter,
                index: None,
                min: Some(min as f64),
                max: max.map(|m| m as f64),
                exclusive_min: None,
                exclusive_max: None,
                got: got as f64,
            }
        };
        let pw = |k, p, n, e, s, a| two_level_factorial_power(k, p, n, e, s, a, None);
        assert_eq!(pw(0, 0, 1, 1.0, 1.0, 0.05), Err(range("k", 1, None, 0)));
        assert_eq!(pw(3, 3, 1, 1.0, 1.0, 0.05), Err(range("p", 0, Some(2), 3)));
        assert_eq!(
            pw(3, 0, 0, 1.0, 1.0, 0.05),
            Err(range("n_replicates", 1, None, 0))
        );
        let open =
            |parameter: &'static str, high: Option<f64>, got: f64| DoeError::ParameterOutOfRange {
                parameter,
                index: None,
                min: None,
                max: None,
                exclusive_min: Some(0.0),
                exclusive_max: high,
                got,
            };
        assert_eq!(
            pw(3, 0, 2, 0.0, 1.0, 0.05),
            Err(open("effect_size", None, 0.0))
        );
        assert_eq!(pw(3, 0, 2, 1.0, -1.0, 0.05), Err(open("sigma", None, -1.0)));
        assert_eq!(
            pw(3, 0, 2, f64::NAN, 1.0, 0.05),
            Err(DoeError::ValueNotFinite {
                parameter: "effect_size",
                index: None
            })
        );
        for alpha in [0.0, 1.0, 2.0] {
            assert_eq!(
                pw(3, 0, 2, 1.0, 1.0, alpha),
                Err(open("alpha", Some(1.0), alpha))
            );
        }
        assert_eq!(
            pw(3, 0, 2, 1.0, 1.0, f64::NAN),
            Err(DoeError::ValueNotFinite {
                parameter: "alpha",
                index: None
            })
        );
        assert_eq!(
            required_replicates(3, 0, 1.0, 1.0, 0.05, 1.0, 10, None),
            Err(open("target_power", Some(1.0), 1.0))
        );
        assert_eq!(
            power_curve(3, 0, 1.0, 1.0, 0.05, 0, None),
            Err(range("max_replicates", 1, None, 0))
        );
    }
}
