//! Statistical power analysis for factorial designs.
//!
//! Computes the power of 2^(k-p) fractional factorial designs to detect
//! a specified effect size using the normal approximation.

use u_numflow::special::{inverse_normal_cdf, standard_normal_cdf};

/// Compute the statistical power of a 2^(k-p) factorial design.
///
/// Uses the normal approximation:
/// - SE_effect = σ × 2 / √(n_replicates × 2^(k-p))
/// - ncp = effect_size / SE_effect
/// - power ≈ Φ(ncp − z_{α/2}) + Φ(−ncp − z_{α/2})
///
/// # Arguments
/// * `k` — total number of factors
/// * `p` — number of generators (k-p = resolution; p=0 means full factorial)
/// * `n_replicates` — number of replicates of the design
/// * `effect_size` — the detectable effect δ (in response units)
/// * `sigma` — process standard deviation σ
/// * `alpha` — type-I error rate (e.g., 0.05)
///
/// # Returns
/// Power in [0, 1], or 0.0 for invalid inputs (k=0, p≥k, n_replicates=0,
/// effect_size≤0, sigma≤0).
pub fn two_level_factorial_power(
    k: usize,
    p: usize,
    n_replicates: usize,
    effect_size: f64,
    sigma: f64,
    alpha: f64,
) -> f64 {
    if k == 0 || p >= k || n_replicates == 0 || effect_size <= 0.0 || sigma <= 0.0 {
        return 0.0;
    }
    let runs_per_rep = 1usize << (k - p); // 2^(k-p)
    let total_runs = (n_replicates * runs_per_rep) as f64;
    let se_effect = sigma * 2.0 / total_runs.sqrt();
    let ncp = effect_size / se_effect;
    let z_alpha_half = inverse_normal_cdf(1.0 - alpha / 2.0);
    let power = standard_normal_cdf(ncp - z_alpha_half) + standard_normal_cdf(-ncp - z_alpha_half);
    power.clamp(0.0, 1.0)
}

/// Find the minimum number of replicates to achieve `target_power`.
///
/// Searches `1..=max_replicates`. Returns `max_replicates` if the target
/// is not reached within the search range (or 0 for invalid inputs).
///
/// # Arguments
/// * `k`, `p`, `effect_size`, `sigma`, `alpha` — as in [`two_level_factorial_power`]
/// * `target_power` — desired power (e.g., 0.80)
/// * `max_replicates` — upper bound for the search
pub fn required_replicates(
    k: usize,
    p: usize,
    effect_size: f64,
    sigma: f64,
    alpha: f64,
    target_power: f64,
    max_replicates: usize,
) -> usize {
    if k == 0 || p >= k || max_replicates == 0 || effect_size <= 0.0 || sigma <= 0.0 {
        return 0;
    }
    for n in 1..=max_replicates {
        let pw = two_level_factorial_power(k, p, n, effect_size, sigma, alpha);
        if pw >= target_power {
            return n;
        }
    }
    max_replicates
}

/// Power curve: `(replicates, power)` pairs for `1..=max_replicates`.
///
/// Useful for visualising how power grows with the number of replicates.
pub fn power_curve(
    k: usize,
    p: usize,
    effect_size: f64,
    sigma: f64,
    alpha: f64,
    max_replicates: usize,
) -> Vec<(usize, f64)> {
    (1..=max_replicates)
        .map(|n| {
            let pw = two_level_factorial_power(k, p, n, effect_size, sigma, alpha);
            (n, pw)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn power_increases_with_replicates() {
        let p1 = two_level_factorial_power(3, 0, 1, 2.0, 1.0, 0.05);
        let p2 = two_level_factorial_power(3, 0, 2, 2.0, 1.0, 0.05);
        let p3 = two_level_factorial_power(3, 0, 3, 2.0, 1.0, 0.05);
        assert!(p2 > p1 && p3 > p2, "p1={p1} p2={p2} p3={p3}");
    }

    #[test]
    fn power_increases_with_effect_size() {
        let small = two_level_factorial_power(4, 0, 1, 0.5, 1.0, 0.05);
        let large = two_level_factorial_power(4, 0, 1, 3.0, 1.0, 0.05);
        assert!(large > small, "small={small} large={large}");
    }

    #[test]
    fn power_in_range() {
        let p = two_level_factorial_power(4, 0, 2, 2.0, 1.0, 0.05);
        assert!((0.0..=1.0).contains(&p), "p={p}");
    }

    #[test]
    fn required_replicates_monotone() {
        let r80 = required_replicates(3, 0, 2.0, 1.0, 0.05, 0.80, 10);
        let r90 = required_replicates(3, 0, 2.0, 1.0, 0.05, 0.90, 10);
        assert!(r90 >= r80, "r80={r80} r90={r90}");
    }

    #[test]
    fn power_curve_non_decreasing() {
        let curve = power_curve(3, 0, 2.0, 1.0, 0.05, 8);
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

    #[test]
    fn invalid_inputs_return_zero() {
        assert_eq!(two_level_factorial_power(0, 0, 1, 1.0, 1.0, 0.05), 0.0);
        assert_eq!(two_level_factorial_power(3, 3, 1, 1.0, 1.0, 0.05), 0.0); // p >= k
        assert_eq!(two_level_factorial_power(3, 0, 0, 1.0, 1.0, 0.05), 0.0);
    }

    #[test]
    fn power_fractional_vs_full() {
        // Fractional (fewer runs) should have lower power than full factorial
        // at the same replicate count.
        let full = two_level_factorial_power(4, 0, 1, 1.5, 1.0, 0.05);
        let frac = two_level_factorial_power(4, 1, 1, 1.5, 1.0, 0.05);
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
        let n = required_replicates(1, 0, 1.0, 1.0, 0.05, 0.80, 50);
        assert_eq!(n, 16, "expected 16 replicates, got {n}");
    }

    /// Power at n=16 replicates (k=1, p=0, δ=σ=1, α=0.05) must be ≥ 0.80.
    /// ncp = 1·√(16/2)/1 = √8 ≈ 2.828; Power = Φ(2.828−1.960) ≈ 0.807
    #[test]
    fn power_at_n16_above_080() {
        let pw = two_level_factorial_power(1, 0, 16, 1.0, 1.0, 0.05);
        assert!(pw >= 0.80, "expected power ≥ 0.80 at n=16, got {pw}");
        // And not too far above (sanity check: below 0.90)
        assert!(pw < 0.90, "power at n=16 suspiciously high: {pw}");
    }

    /// Power at n=15 replicates must be below 0.80 (n=16 is the minimum).
    #[test]
    fn power_at_n15_below_080() {
        let pw = two_level_factorial_power(1, 0, 15, 1.0, 1.0, 0.05);
        assert!(pw < 0.80, "expected power < 0.80 at n=15, got {pw}");
    }
}
