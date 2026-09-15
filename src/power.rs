//! Statistical power analysis for factorial designs.
//!
//! Computes the power of 2^(k-p) fractional factorial designs to detect
//! a specified effect size using the normal approximation.

use u_numflow::special::{inverse_normal_cdf, standard_normal_cdf};

use crate::error::DoeError;

/// Compute the statistical power of a 2^(k-p) factorial design.
///
/// Uses the normal approximation:
/// - SE_effect = σ × 2 / √(n_replicates × 2^(k-p))
/// - ncp = effect_size / SE_effect
/// - power ≈ Φ(ncp − z_{α/2}) + Φ(−ncp − z_{α/2})
///
/// # Arguments
/// * `k` — total number of factors
/// * `p` — number of generators (the fraction is 2^(-p); p = 0 is a full factorial)
/// * `n_replicates` — number of replicates of the design
/// * `effect_size` — the detectable effect δ (in response units)
/// * `sigma` — process standard deviation σ
/// * `alpha` — type-I error rate (e.g., 0.05)
///
/// # Errors
/// [`DoeError::ParameterOutOfRange`] if `k` is 0, `p >= k` or `n_replicates`
/// is 0; [`DoeError::ValueOutOfDomain`] unless `effect_size` and `sigma` are
/// finite and positive and `alpha` lies strictly between 0 and 1. Such inputs
/// have no power to report -- a number, even 0, would read as one.
///
/// # Examples
///
/// ```
/// use u_doe::power::two_level_factorial_power;
/// let power = two_level_factorial_power(1, 0, 16, 1.0, 1.0, 0.05).unwrap();
/// assert!(power > 0.80 && power < 0.81);
/// assert!(two_level_factorial_power(3, 0, 1, 2.0, 0.0, 0.05).is_err());
/// ```
pub fn two_level_factorial_power(
    k: usize,
    p: usize,
    n_replicates: usize,
    effect_size: f64,
    sigma: f64,
    alpha: f64,
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
    Ok(power_unchecked(
        k,
        p,
        n_replicates,
        effect_size,
        sigma,
        alpha,
    ))
}

/// Find the minimum number of replicates to achieve `target_power`.
///
/// Searches `1..=max_replicates` and returns the first replicate count whose
/// power reaches `target_power`, or `None` if none in the range does.
///
/// # Arguments
/// * `k`, `p`, `effect_size`, `sigma`, `alpha` — as in [`two_level_factorial_power`]
/// * `target_power` — desired power (e.g., 0.80)
/// * `max_replicates` — upper bound for the search
///
/// # Errors
/// As [`two_level_factorial_power`], and [`DoeError::ParameterOutOfRange`] if
/// `max_replicates` is 0 or [`DoeError::ValueOutOfDomain`] unless
/// `target_power` lies strictly between 0 and 1.
///
/// # Examples
///
/// ```
/// use u_doe::power::required_replicates;
/// assert_eq!(required_replicates(1, 0, 1.0, 1.0, 0.05, 0.80, 50).unwrap(), Some(16));
/// assert_eq!(required_replicates(1, 0, 1.0, 1.0, 0.05, 0.80, 10).unwrap(), None);
/// ```
pub fn required_replicates(
    k: usize,
    p: usize,
    effect_size: f64,
    sigma: f64,
    alpha: f64,
    target_power: f64,
    max_replicates: usize,
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
    Ok((1..=max_replicates)
        .find(|&n| power_unchecked(k, p, n, effect_size, sigma, alpha) >= target_power))
}

/// Power curve: `(replicates, power)` pairs for `1..=max_replicates`.
///
/// Useful for visualising how power grows with the number of replicates.
///
/// # Errors
/// As [`two_level_factorial_power`], and [`DoeError::ParameterOutOfRange`] if
/// `max_replicates` is 0.
pub fn power_curve(
    k: usize,
    p: usize,
    effect_size: f64,
    sigma: f64,
    alpha: f64,
    max_replicates: usize,
) -> Result<Vec<(usize, f64)>, DoeError> {
    check_design(k, p)?;
    check_max_replicates(max_replicates)?;
    check_effect(effect_size, sigma, alpha)?;
    Ok((1..=max_replicates)
        .map(|n| (n, power_unchecked(k, p, n, effect_size, sigma, alpha)))
        .collect())
}

const POSITIVE: &str = "finite and greater than 0";
const OPEN_UNIT_INTERVAL: &str = "strictly between 0 and 1";

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

/// Power for inputs already checked.
fn power_unchecked(
    k: usize,
    p: usize,
    n_replicates: usize,
    effect_size: f64,
    sigma: f64,
    alpha: f64,
) -> f64 {
    // 2^(k-p) as a float: an integer shift overflows once k - p reaches the
    // width of `usize`, which is 32 bits in WebAssembly.
    let total_runs = n_replicates as f64 * 2f64.powf((k - p) as f64);
    let se_effect = sigma * 2.0 / total_runs.sqrt();
    let ncp = effect_size / se_effect;
    let z_alpha_half = inverse_normal_cdf(1.0 - alpha / 2.0);
    let power = standard_normal_cdf(ncp - z_alpha_half) + standard_normal_cdf(-ncp - z_alpha_half);
    power.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn power(k: usize, p: usize, n: usize, effect: f64, sigma: f64, alpha: f64) -> f64 {
        two_level_factorial_power(k, p, n, effect, sigma, alpha).expect("valid power inputs")
    }

    #[test]
    fn power_increases_with_replicates() {
        let p1 = power(3, 0, 1, 2.0, 1.0, 0.05);
        let p2 = power(3, 0, 2, 2.0, 1.0, 0.05);
        let p3 = power(3, 0, 3, 2.0, 1.0, 0.05);
        assert!(p2 > p1 && p3 > p2, "p1={p1} p2={p2} p3={p3}");
    }

    #[test]
    fn power_increases_with_effect_size() {
        let small = power(4, 0, 1, 0.5, 1.0, 0.05);
        let large = power(4, 0, 1, 3.0, 1.0, 0.05);
        assert!(large > small, "small={small} large={large}");
    }

    #[test]
    fn power_in_range() {
        let p = power(4, 0, 2, 2.0, 1.0, 0.05);
        assert!((0.0..=1.0).contains(&p), "p={p}");
    }

    #[test]
    fn required_replicates_monotone() {
        let r80 = required_replicates(3, 0, 2.0, 1.0, 0.05, 0.80, 10).unwrap();
        let r90 = required_replicates(3, 0, 2.0, 1.0, 0.05, 0.90, 10).unwrap();
        assert!(r90 >= r80, "r80={r80:?} r90={r90:?}");
    }

    #[test]
    fn power_curve_non_decreasing() {
        let curve = power_curve(3, 0, 2.0, 1.0, 0.05, 8).unwrap();
        assert_eq!(curve.len(), 8);
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
        assert_eq!(
            two_level_factorial_power(0, 0, 1, 1.0, 1.0, 0.05),
            Err(range("k", 1, None, 0))
        );
        assert_eq!(
            two_level_factorial_power(3, 3, 1, 1.0, 1.0, 0.05),
            Err(range("p", 0, Some(2), 3))
        );
        assert_eq!(
            two_level_factorial_power(3, 0, 0, 1.0, 1.0, 0.05),
            Err(range("n_replicates", 1, None, 0))
        );
        let domain = |f: Result<f64, DoeError>| match f {
            Err(DoeError::ValueOutOfDomain { parameter, .. }) => parameter,
            other => panic!("{other:?}"),
        };
        assert_eq!(
            domain(two_level_factorial_power(3, 0, 1, 0.0, 1.0, 0.05)),
            "effect_size"
        );
        assert_eq!(
            domain(two_level_factorial_power(3, 0, 1, 1.0, -1.0, 0.05)),
            "sigma"
        );
        assert_eq!(
            domain(two_level_factorial_power(3, 0, 1, f64::NAN, 1.0, 0.05)),
            "effect_size"
        );
        for alpha in [0.0, 1.0, 2.0, f64::NAN] {
            assert_eq!(
                domain(two_level_factorial_power(3, 0, 1, 1.0, 1.0, alpha)),
                "alpha"
            );
        }
        assert!(matches!(
            required_replicates(3, 0, 1.0, 1.0, 0.05, 1.0, 10),
            Err(DoeError::ValueOutOfDomain {
                parameter: "target_power",
                ..
            })
        ));
        assert_eq!(
            power_curve(3, 0, 1.0, 1.0, 0.05, 0),
            Err(range("max_replicates", 1, None, 0))
        );
    }

    /// `2^(k-p)` used to be an integer shift, which overflows at the width of
    /// `usize`. A design that large has power 1 for any real effect.
    #[test]
    fn a_very_large_design_does_not_overflow() {
        assert_eq!(power(70, 0, 1, 1.0, 1.0, 0.05), 1.0);
    }

    #[test]
    fn power_fractional_vs_full() {
        // Fractional (fewer runs) should have lower power than full factorial
        // at the same replicate count.
        let full = power(4, 0, 1, 1.5, 1.0, 0.05);
        let frac = power(4, 1, 1, 1.5, 1.0, 0.05);
        assert!(full > frac, "full={full} frac={frac}");
    }

    // -----------------------------------------------------------------------
    // Reference-value tests (Task 12)
    // Cohen (1988) §2; Montgomery (2020) §3.7
    //
    // For k=1, p=0 (2^1 = 2 runs per replicate):
    //   SE_effect = σ·2 / √(n·2)  =  σ / √(n/2)
    //   ncp       = δ / SE_effect  = δ·√(n/2) / σ
    //
    // Setting ncp = z_{α/2} + z_β  with α=0.05, power=0.8:
    //   z_{α/2}  = 1.960,  z_β = 0.842
    //   ncp²    = (1.960 + 0.842)² = 7.851
    //   n/2     = 7.851  →  n = 15.7  →  ceil = 16 replicates
    // -----------------------------------------------------------------------

    /// required_replicates for k=1, α=0.05, power=0.8, δ=σ=1 must equal 16.
    /// Reference: n = 2·(z_{α/2}+z_β)²·(σ/δ)² = 2·(1.960+0.842)² ≈ 15.7 → 16
    #[test]
    fn required_replicates_reference_n16() {
        let n = required_replicates(1, 0, 1.0, 1.0, 0.05, 0.80, 50).unwrap();
        assert_eq!(n, Some(16), "expected 16 replicates, got {n:?}");
    }

    /// Reaching the target exactly at `max_replicates` and not reaching it at
    /// all used to return the same number.
    #[test]
    fn required_replicates_tells_reached_at_the_bound_from_not_reached() {
        assert_eq!(
            required_replicates(1, 0, 1.0, 1.0, 0.05, 0.80, 16).unwrap(),
            Some(16)
        );
        assert_eq!(
            required_replicates(1, 0, 1.0, 1.0, 0.05, 0.80, 15).unwrap(),
            None
        );
    }

    /// Power at n=16 replicates (k=1, p=0, δ=σ=1, α=0.05) must be ≥ 0.80.
    /// ncp = 1·√(16/2)/1 = √8 ≈ 2.828; Power = Φ(2.828−1.960) ≈ 0.807
    #[test]
    fn power_at_n16_above_080() {
        let pw = power(1, 0, 16, 1.0, 1.0, 0.05);
        assert!(pw >= 0.80, "expected power ≥ 0.80 at n=16, got {pw}");
        // And not too far above (sanity check: below 0.90)
        assert!(pw < 0.90, "power at n=16 suspiciously high: {pw}");
    }

    /// Power at n=15 replicates must be below 0.80 (n=16 is the minimum).
    #[test]
    fn power_at_n15_below_080() {
        let pw = power(1, 0, 15, 1.0, 1.0, 0.05);
        assert!(pw < 0.80, "expected power < 0.80 at n=15, got {pw}");
    }
}
