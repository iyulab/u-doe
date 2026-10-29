//! A response written in another unit gives the same analysis: every ratio an
//! analysis reports -- F, p, t, the share of each effect, Lenth's margin of
//! error in units of the pseudo standard error -- is the same whether the
//! responses are in picounits or teraunits.

use u_doe::analysis::anova::doe_anova;
use u_doe::analysis::effects::estimate_effects;
use u_doe::analysis::least_squares::fit_least_squares;
use u_doe::analysis::lenth::lenth;
use u_doe::design::factorial::full_factorial;

const SCALES: [f64; 7] = [1e-12, 1e-9, 1e-6, 1e-3, 1e3, 1e9, 1e12];

fn responses() -> Vec<f64> {
    // 2^3 with one clear main effect, one weaker, and noise.
    let mut state = 20261029u64;
    let design = full_factorial(3).expect("2^3");
    (0..design.run_count())
        .map(|i| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let noise = ((state >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 0.6;
            20.0 + 3.0 * design.get(i, 0) + 0.8 * design.get(i, 1) + noise
        })
        .collect()
}

fn close(a: f64, b: f64, what: &str, k: f64) {
    assert!(
        (a - b).abs() <= 1e-7 * b.abs().max(1e-3),
        "{what} at scale {k:e}: {a} vs {b}"
    );
}

#[test]
fn the_analyses_give_the_same_ratios_in_any_unit() {
    let design = full_factorial(3).expect("2^3");
    let y = responses();
    let names = ["A", "B", "C", "A:B"];
    let base_anova = doe_anova(&design, &y, &names).expect("anova");
    let base_ls = fit_least_squares(&design, &y, &names).expect("fit");
    let base_effects = estimate_effects(&design, &y, 3).expect("effects");
    let base_lenth = lenth(&design, &base_effects).expect("lenth");

    for k in SCALES {
        let yk: Vec<f64> = y.iter().map(|v| v * k).collect();
        let anova = doe_anova(&design, &yk, &names).expect("anova");
        close(anova.r_squared, base_anova.r_squared, "R2", k);
        for (a, b) in anova.effects.iter().zip(&base_anova.effects) {
            assert_eq!(
                a.f_statistic.is_some(),
                b.f_statistic.is_some(),
                "{} F present at {k:e}",
                a.name
            );
            if let (Some(x), Some(z)) = (a.f_statistic, b.f_statistic) {
                close(x, z, &format!("{} F", a.name), k);
            }
            if let (Some(x), Some(z)) = (a.p_value, b.p_value) {
                close(x, z, &format!("{} p", a.name), k);
            }
        }
        let ls = fit_least_squares(&design, &yk, &names).expect("fit");
        for (a, b) in ls.terms.iter().zip(&base_ls.terms) {
            if let (Some(x), Some(z)) = (a.t_statistic, b.t_statistic) {
                close(x, z, &format!("{} t", a.name), k);
            }
        }
        let effects = estimate_effects(&design, &yk, 3).expect("effects");
        for (a, b) in effects.iter().zip(&base_effects) {
            close(
                a.percent_contribution,
                b.percent_contribution,
                &format!("{} share", a.name),
                k,
            );
        }
        let l = lenth(&design, &effects).expect("lenth");
        close(
            l.margin_of_error / l.pse,
            base_lenth.margin_of_error / base_lenth.pse,
            "ME/PSE",
            k,
        );
        for (a, b) in effects.iter().zip(&base_effects) {
            close(
                a.estimate.abs() / l.pse,
                b.estimate.abs() / base_lenth.pse,
                &format!("{} |effect|/PSE", a.name),
                k,
            );
        }
    }
}
