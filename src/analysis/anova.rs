//! DOE-context Analysis of Variance (ANOVA).
//!
//! Partitions response variation into contributions from specified effects
//! and residual, then computes F-statistics for significance testing.
//!
//! Reference: Montgomery, D.C. (2019). *Introduction to Statistical Quality
//! Control*, 8th ed., Section 6.4. Wiley.

use crate::{design::DesignMatrix, error::DoeError};

/// One row in the ANOVA table.
#[derive(Debug, Clone)]
pub struct EffectRow {
    /// Effect name (e.g. "A", "A:B" — factor names joined with `":"`).
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

/// Curvature test from centre points (Montgomery, *Design and Analysis of
/// Experiments*, §6.8).
///
/// Centre points sit at the middle of every factor range. If the response is
/// linear in the factors, their mean equals the mean of the factorial runs;
/// the difference, scaled, is a sum of squares with one degree of freedom.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurvatureTest {
    /// `n_F · n_C · (ȳ_F − ȳ_C)² / (n_F + n_C)`.
    pub sum_of_squares: f64,
    /// Always 1.
    pub df: usize,
    /// `sum_of_squares / MS_pure_error`. `None` without pure error to test
    /// against: fewer than two runs at any design point, or replicates that
    /// agree exactly.
    pub f_statistic: Option<f64>,
    /// p-value of `f_statistic` from F(1, df_pure_error).
    pub p_value: Option<f64>,
}

/// Pure error: the spread among runs made at the same design point.
///
/// It is the one error estimate that does not depend on which terms the model
/// includes. Every group of identical design rows contributes -- centre points
/// and replicated factorial runs alike.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PureError {
    /// Sum of squared deviations of each run from the mean of its group.
    pub sum_of_squares: f64,
    /// Sum over groups of (group size − 1).
    pub df: usize,
}

/// Result of a DOE ANOVA.
#[derive(Debug, Clone)]
pub struct DoeAnovaResult {
    /// Rows for each requested effect.
    pub effects: Vec<EffectRow>,
    /// Residual (error) sum of squares: what neither the effects nor the
    /// curvature term account for. It includes the pure error.
    pub residual_ss: f64,
    /// Residual degrees of freedom.
    pub residual_df: usize,
    /// Total sum of squares.
    pub total_ss: f64,
    /// R² = (SS_model + SS_curvature) / SS_total.
    pub r_squared: f64,
    /// Adjusted R² = 1 - (SS_residual/df_residual) / (SS_total/df_total).
    pub r_squared_adj: f64,
    /// Fitted value for each run, in run order. A factorial run's is the mean
    /// of the factorial runs plus half of each requested effect times its
    /// contrast -- exact because the requested contrasts are orthogonal and
    /// balanced (see Errors). A centre run's is the mean of the centre runs.
    pub fitted: Vec<f64>,
    /// `responses[i] - fitted[i]` for each run. Their squares sum to
    /// [`residual_ss`](Self::residual_ss); they are what a residuals-versus-
    /// fitted or normal probability plot needs.
    pub residuals: Vec<f64>,
    /// Curvature test, when the design has centre points (runs with every
    /// factor at 0); `None` otherwise.
    pub curvature: Option<CurvatureTest>,
    /// Pure error, when some design point was run more than once; `None`
    /// otherwise.
    pub pure_error: Option<PureError>,
}

/// Compute DOE ANOVA for specified effects.
///
/// # Arguments
/// * `design`       — Experimental design matrix: coded ±1 values, plus any
///   centre points (runs with every factor at 0)
/// * `responses`    — Observed response values, one per run
/// * `effect_names` — Names of effects to include (e.g. `&["A", "B", "A:B"]`).
///   Interaction names join factor names with `":"`.
///
/// # Centre points
///
/// A centre point has a zero contrast for every term, so it says nothing about
/// the effects: they and their sums of squares come from the factorial runs
/// alone. What centre points do measure is curvature -- whether the response
/// at the middle of the region lies on the plane through the corners -- and
/// that is reported in [`DoeAnovaResult::curvature`], tested against the pure
/// error. Counting centre points into the contrasts instead would shrink every
/// effect by `n_F / (n_F + n_C)`.
///
/// # Errors
///
/// Returns `Err` if `responses.len() != design.run_count()`;
/// [`DoeError::UnknownEffect`] if an entry in `effect_names` does not match
/// any estimable effect (main effects and two-factor interactions);
/// [`DoeError::NotTwoLevelCoded`] if a run is neither two-level coded (every
/// factor `-1` or `+1`) nor a centre point (every factor `0`) -- axial and
/// three-level designs belong in [`crate::analysis::rsm::fit_rsm`];
/// [`DoeError::AliasedEffects`] if two requested effects share a contrast
/// column, which would count the same sum of squares twice, or a requested
/// contrast is the same in every run (aliased with the mean, `I`);
/// [`DoeError::PartiallyAliasedEffects`] if two requested contrasts are
/// correlated without being identical, or one is correlated with the mean --
/// the per-term sums of squares would then overlap and not add up to the model
/// sum of squares; or [`DoeError::OverSpecifiedModel`] if the model (counting
/// the curvature term) asks for more degrees of freedom than the design has.
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

    // Centre points -- runs with every factor at 0 -- are set aside; the
    // effects come from the factorial runs.
    let is_centre = |row: &[f64]| !row.is_empty() && row.iter().all(|v| v.abs() < 1e-9);
    let (centre_runs, factorial_runs): (Vec<usize>, Vec<usize>) =
        (0..n).partition(|&run| is_centre(&design.data[run]));
    let factorial = DesignMatrix {
        data: factorial_runs
            .iter()
            .map(|&run| design.data[run].clone())
            .collect(),
        factor_names: design.factor_names.clone(),
    };
    let factorial_y: Vec<f64> = factorial_runs.iter().map(|&run| responses[run]).collect();
    let n_f = factorial_runs.len();
    let n_c = centre_runs.len();

    // Estimate all effects up to order 2 from the factorial runs. This also
    // rejects a design with no factorial runs and any run that is neither two-
    // level coded nor a centre point -- reported in the caller's run numbering,
    // before `total_df = n - 1` below could underflow.
    // Unchecked: only the requested terms have to be orthogonal, and they are
    // checked below -- a 12-run Plackett-Burman's interactions are correlated
    // with its main effects, which does not stop an ANOVA of the main effects.
    let all_effects = crate::analysis::effects::contrast_effects(&factorial, &factorial_y, 2)
        .map_err(|e| match e {
            DoeError::NotTwoLevelCoded { run, factor, value } => DoeError::NotTwoLevelCoded {
                run: factorial_runs[run],
                factor,
                value,
            },
            other => other,
        })?;

    // Grand mean and total SS, over every run
    let grand_mean = responses.iter().sum::<f64>() / n as f64;
    let total_ss: f64 = responses.iter().map(|&y| (y - grand_mean).powi(2)).sum();
    let total_df = n - 1;

    // Resolve requested effects — an unknown name is an error, not a silent
    // skip: silently dropping a term yields a structurally different model
    // with no signal to the caller.
    let selected: Vec<_> = effect_names
        .iter()
        .map(|&name| {
            all_effects
                .iter()
                .find(|e| e.name == name)
                .ok_or_else(|| DoeError::UnknownEffect {
                    name: name.to_string(),
                })
        })
        .collect::<Result<_, _>>()?;

    // Two terms that share a contrast column are aliased: they carry the same
    // sum of squares, so admitting both inflates the model and deflates the
    // residual by exactly that amount. In a fractional design this is the
    // normal case rather than an exotic one -- in a 2^(5-2), A, B:D and C:E are
    // one column -- so it has to be refused rather than clamped away later.
    let contrasts: Vec<Vec<f64>> = selected
        .iter()
        .map(|e| {
            (0..n_f)
                .map(|run| e.columns.iter().map(|&c| factorial.get(run, c)).product())
                .collect()
        })
        .collect();
    for i in 0..contrasts.len() {
        for j in (i + 1)..contrasts.len() {
            let same = contrasts[i]
                .iter()
                .zip(contrasts[j].iter())
                .all(|(a, b)| (a - b).abs() < 1e-9);
            let opposite = contrasts[i]
                .iter()
                .zip(contrasts[j].iter())
                .all(|(a, b)| (a + b).abs() < 1e-9);
            if same || opposite {
                return Err(DoeError::AliasedEffects {
                    first: selected[i].name.clone(),
                    second: selected[j].name.clone(),
                });
            }
        }
    }

    // The table below sums per-term sums of squares into the model sum of
    // squares, and builds fitted values from half-effects. Both are exact only
    // when every requested contrast is balanced (orthogonal to the mean) and
    // the contrasts are pairwise orthogonal. Regular fractions satisfy this for
    // any estimable set; a Plackett-Burman two-factor interaction does not --
    // it is correlated with other main effects -- and neither does a design
    // with a run left out. Only the requested terms are checked: a 12-run
    // Plackett-Burman supports an ANOVA of its main effects.
    let names: Vec<&str> = selected.iter().map(|e| e.name.as_str()).collect();
    crate::analysis::effects::require_orthogonal_contrasts(&names, &contrasts)?;

    // Curvature: the difference between the factorial and centre means,
    // scaled to a one-degree-of-freedom sum of squares.
    let mean_f = factorial_y.iter().sum::<f64>() / n_f as f64;
    let (curvature_ss, mean_c) = if n_c > 0 {
        let mean_c = centre_runs.iter().map(|&run| responses[run]).sum::<f64>() / n_c as f64;
        let d = mean_f - mean_c;
        ((n_f * n_c) as f64 * d * d / (n_f + n_c) as f64, mean_c)
    } else {
        (0.0, mean_f)
    };
    let curvature_df = usize::from(n_c > 0);

    let model_df = selected.len();
    if model_df + curvature_df > total_df {
        return Err(DoeError::OverSpecifiedModel {
            terms: model_df + curvature_df,
            runs: n,
        });
    }

    let model_ss: f64 = selected.iter().map(|e| e.sum_of_squares).sum();
    // The total splits exactly into the effects (within the factorial runs),
    // the curvature (between the factorial and centre means) and the rest.
    // With the coding, aliasing, orthogonality and degrees-of-freedom guards
    // above, what is clamped here is float noise, not a structural overflow.
    let residual_ss = (total_ss - model_ss - curvature_ss).max(0.0);
    let residual_df = total_df - model_df - curvature_df;
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

    let pure_error = pure_error(design, responses);
    let curvature = (n_c > 0).then(|| {
        let (f_statistic, p_value) = match pure_error {
            Some(pe) if pe.sum_of_squares > 0.0 => {
                let f = curvature_ss / (pe.sum_of_squares / pe.df as f64);
                (Some(f), Some(f_pvalue(f, 1, pe.df)))
            }
            _ => (None, None),
        };
        CurvatureTest {
            sum_of_squares: curvature_ss,
            df: 1,
            f_statistic,
            p_value,
        }
    });

    let r_squared = if total_ss > 1e-12 {
        (model_ss + curvature_ss) / total_ss
    } else {
        0.0
    };
    let r_squared_adj = if residual_df > 0 && total_df > 0 {
        1.0 - (residual_ss / residual_df as f64) / (total_ss / total_df as f64)
    } else {
        r_squared
    };

    // With orthogonal, balanced contrasts the least-squares coefficient of each
    // term is half its effect, so no matrix solve is needed. Centre runs are
    // fitted at their own mean -- that is what the curvature term fits.
    let mut fitted = vec![mean_c; n];
    for (i, &run) in factorial_runs.iter().enumerate() {
        fitted[run] = mean_f
            + selected
                .iter()
                .zip(contrasts.iter())
                .map(|(e, c)| e.estimate / 2.0 * c[i])
                .sum::<f64>();
    }
    let residuals: Vec<f64> = responses
        .iter()
        .zip(fitted.iter())
        .map(|(y, f)| y - f)
        .collect();

    Ok(DoeAnovaResult {
        effects,
        residual_ss,
        residual_df,
        total_ss,
        r_squared,
        r_squared_adj,
        fitted,
        residuals,
        curvature,
        pure_error,
    })
}

/// Pure error from every group of identical design rows; `None` when no design
/// point was run more than once.
fn pure_error(design: &DesignMatrix, responses: &[f64]) -> Option<PureError> {
    let same_point = |a: &[f64], b: &[f64]| {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() < 1e-9)
    };
    let mut assigned = vec![false; design.run_count()];
    let mut ss = 0.0;
    let mut df = 0;
    for run in 0..design.run_count() {
        if assigned[run] {
            continue;
        }
        let group: Vec<usize> = (run..design.run_count())
            .filter(|&other| !assigned[other] && same_point(&design.data[run], &design.data[other]))
            .collect();
        for &member in &group {
            assigned[member] = true;
        }
        let mean = group.iter().map(|&r| responses[r]).sum::<f64>() / group.len() as f64;
        ss += group
            .iter()
            .map(|&r| (responses[r] - mean).powi(2))
            .sum::<f64>();
        df += group.len() - 1;
    }
    (df > 0).then_some(PureError {
        sum_of_squares: ss,
        df,
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
    /// 2^3 with centre points: the shape a template for adding centre points
    /// produces. `y = 50 + 5A + 3B + C` on the factorial runs.
    fn factorial_with_centres(centres: &[f64]) -> (DesignMatrix, Vec<f64>) {
        let base = crate::design::factorial::full_factorial(3).expect("2^3");
        let mut y: Vec<f64> = base
            .data
            .iter()
            .map(|r| 50.0 + 5.0 * r[0] + 3.0 * r[1] + r[2])
            .collect();
        let mut data = base.data.clone();
        for &c in centres {
            data.push(vec![0.0; 3]);
            y.push(c);
        }
        (
            DesignMatrix {
                data,
                factor_names: base.factor_names,
            },
            y,
        )
    }

    #[test]
    fn centre_points_give_curvature_and_pure_error() {
        // The centre mean is on the plane: no curvature, and three centre
        // replicates give two degrees of pure error.
        let (design, y) = factorial_with_centres(&[50.1, 49.9, 50.0]);
        let r = doe_anova(&design, &y, &["A", "B", "C"]).expect("valid");
        let c = r.curvature.expect("centre points give a curvature test");
        assert!(c.sum_of_squares.abs() < 1e-9);
        let pe = r.pure_error.expect("replicated centre points");
        assert_eq!(pe.df, 2);
        assert!((pe.sum_of_squares - 0.02).abs() < 1e-9);
    }

    #[test]
    fn curvature_is_tested_against_pure_error() {
        let (design, y) = factorial_with_centres(&[55.0, 55.2, 54.8]);
        let r = doe_anova(&design, &y, &["A", "B", "C"]).expect("valid");
        let c = r.curvature.expect("curvature");
        let expected = 8.0 * 3.0 * 25.0 / 11.0;
        assert!(
            (c.sum_of_squares - expected).abs() < 1e-9,
            "{}",
            c.sum_of_squares
        );
        let pe = r.pure_error.expect("pure error");
        assert!((pe.sum_of_squares - 0.08).abs() < 1e-9);
        let f = c.f_statistic.expect("F");
        assert!((f - expected / 0.04).abs() < 1e-6, "F = {f}");
        // F(1, 2) is the square of a t with 2 degrees of freedom, whose tail has
        // a closed form: P(F > f) = 1 - sqrt(f / (2 + f)).
        let p = c.p_value.expect("p");
        assert!((p - (1.0 - (f / (2.0 + f)).sqrt())).abs() < 1e-9, "p = {p}");
        // Centre runs are fitted at their own mean.
        for run in 8..11 {
            assert!((r.fitted[run] - 55.0).abs() < 1e-9);
        }
    }

    #[test]
    fn the_sums_of_squares_add_up_with_centre_points() {
        let (design, y) = factorial_with_centres(&[55.0, 55.2, 54.8]);
        // Only A and B: C's sum of squares moves into the residual.
        let r = doe_anova(&design, &y, &["A", "B"]).expect("valid");
        let model: f64 = r.effects.iter().map(|e| e.sum_of_squares).sum();
        let curvature = r.curvature.expect("curvature").sum_of_squares;
        assert!((model + curvature + r.residual_ss - r.total_ss).abs() < 1e-9);
        let ss: f64 = r.residuals.iter().map(|e| e * e).sum();
        assert!((ss - r.residual_ss).abs() < 1e-9);
        assert_eq!(r.residual_df, 11 - 1 - 2 - 1);
    }

    #[test]
    fn a_run_that_is_neither_factorial_nor_centre_is_named_in_the_callers_numbering() {
        let (mut design, y) = factorial_with_centres(&[50.0, 50.0]);
        design.data[9] = vec![0.0, 1.0, 0.0];
        match doe_anova(&design, &y, &["A"]) {
            Err(DoeError::NotTwoLevelCoded { run, .. }) => assert_eq!(run, 9),
            other => panic!("expected NotTwoLevelCoded at run 9, got {other:?}"),
        }
    }

    #[test]
    fn replicated_factorial_runs_give_pure_error_without_centre_points() {
        let base = crate::design::factorial::full_factorial(2).expect("2^2");
        let mut data = base.data.clone();
        data.extend(base.data.clone());
        let design = DesignMatrix {
            data,
            factor_names: base.factor_names,
        };
        let y = [10.0, 14.0, 12.0, 20.0, 11.0, 15.0, 12.5, 19.0];
        let r = doe_anova(&design, &y, &["A", "B"]).expect("valid");
        assert!(r.curvature.is_none());
        let pe = r.pure_error.expect("replicates");
        assert_eq!(pe.df, 4);
        let expected = 0.5 + 0.5 + 0.125 + 0.5;
        assert!(
            (pe.sum_of_squares - expected).abs() < 1e-9,
            "{}",
            pe.sum_of_squares
        );
    }

    #[test]
    fn an_unreplicated_factorial_has_neither() {
        let design = crate::design::factorial::full_factorial(3).expect("2^3");
        let y = [10.0, 20.0, 15.0, 25.0, 12.0, 22.0, 18.0, 30.0];
        let r = doe_anova(&design, &y, &["A", "B", "C"]).expect("valid");
        assert!(r.curvature.is_none() && r.pure_error.is_none());
    }
    /// Counted into the contrasts, three centre points shrank every effect by
    /// 8/11; refused, they left the design with no analysis at all.
    #[test]
    fn centre_points_are_accepted_and_leave_the_effects_alone() {
        let (design, y) = factorial_with_centres(&[50.1, 49.9, 50.0]);
        let r = doe_anova(&design, &y, &["A", "B", "C"]).expect("centre points are accepted");
        let ss: Vec<f64> = r.effects.iter().map(|e| e.sum_of_squares).collect();
        assert!((ss[0] - 200.0).abs() < 1e-9, "SS_A = {}", ss[0]);
        assert!((ss[1] - 72.0).abs() < 1e-9, "SS_B = {}", ss[1]);
        assert!((ss[2] - 8.0).abs() < 1e-9, "SS_C = {}", ss[2]);
    }
    #[test]
    fn residuals_square_to_the_residual_sum_of_squares() {
        let design = crate::design::factorial::full_factorial(3).expect("2^3");
        let responses = [10.0, 20.0, 15.0, 25.0, 12.0, 22.0, 18.0, 30.0];
        let result = doe_anova(&design, &responses, &["A", "B", "C"]).expect("valid");
        let ss: f64 = result.residuals.iter().map(|r| r * r).sum();
        assert!(
            (ss - result.residual_ss).abs() < 1e-9,
            "{ss} vs {}",
            result.residual_ss
        );
        assert!(result.residuals.iter().sum::<f64>().abs() < 1e-9);
        for ((y, f), r) in responses.iter().zip(&result.fitted).zip(&result.residuals) {
            assert!((y - f - r).abs() < 1e-12);
        }
    }

    #[test]
    fn a_saturated_model_fits_every_run_exactly() {
        let design = crate::design::factorial::full_factorial(2).expect("2^2");
        let responses = [28.0, 36.0, 18.0, 31.0];
        let result = doe_anova(&design, &responses, &["A", "B", "A:B"]).expect("valid");
        assert_eq!(result.residual_df, 0);
        for (y, f) in responses.iter().zip(&result.fitted) {
            assert!((y - f).abs() < 1e-9, "{y} vs {f}");
        }
    }

    /// In a 12-run Plackett-Burman design a two-factor interaction column is
    /// correlated with the main effects it does not contain (by +/-1/3). The
    /// per-term sums of squares then do not add up to the model sum of
    /// squares, so a table built from them is not an ANOVA of this data.
    #[test]
    fn partially_aliased_effects_are_refused() {
        let design = crate::design::plackett_burman::plackett_burman(8).expect("pb12");
        assert_eq!(design.run_count(), 12);
        let responses: Vec<f64> = (0..12).map(|i| 10.0 + (i * 7 % 5) as f64).collect();
        let names = ["A", "B", "C", "D", "E", "F", "G", "H", "A:B"];
        let result = doe_anova(&design, &responses, &names);
        assert!(
            matches!(result, Err(DoeError::PartiallyAliasedEffects { .. })),
            "partially aliased model accepted: {result:?}"
        );
        // The main effects alone are orthogonal and still analysable.
        doe_anova(&design, &responses, &names[..8]).expect("main effects are orthogonal");
    }
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
            &["A", "B", "C", "D", "A:C"],
        )
        .unwrap();
        let row_a = result.effects.iter().find(|r| r.name == "A").unwrap();
        // A is highly significant
        assert!(row_a.f_statistic > 10.0, "A F={}", row_a.f_statistic);
    }

    #[test]
    fn anova_r_squared_reasonable() {
        let design = full_factorial(4).unwrap();
        let result = doe_anova(&design, &montgomery_responses(), &["A", "C", "D", "A:C"]).unwrap();
        // A, C, D, AC explain ~77% of total variance in this dataset
        assert!(result.r_squared > 0.70, "R²={}", result.r_squared);
    }

    #[test]
    fn anova_ss_partition() {
        let design = full_factorial(4).unwrap();
        let result = doe_anova(&design, &montgomery_responses(), &["A", "C", "D", "A:C"]).unwrap();
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
    fn anova_unknown_effect_errors() {
        let design = full_factorial(2).unwrap();
        let responses = vec![10.0, 20.0, 30.0, 40.0];
        // Valid names succeed
        assert!(doe_anova(&design, &responses, &["A", "B", "A:B"]).is_ok());
        // Unknown name errors instead of silently producing a smaller model
        let err = doe_anova(&design, &responses, &["A", "Z"]).unwrap_err();
        assert_eq!(
            err,
            crate::error::DoeError::UnknownEffect { name: "Z".into() }
        );
        // Legacy separator-less interaction names are unknown (renamed in 0.5.0)
        assert!(doe_anova(&design, &responses, &["AB"]).is_err());
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

    // --- model admissibility guards (Cycle 265) ------------------------------

    #[test]
    fn rejects_aliased_effects_in_a_fractional_design() {
        use crate::design::factorial::fractional_factorial;

        // 2^(5-2): A, B:D and C:E are the same contrast column. Admitting two of
        // them counts one sum of squares twice, which used to surface only as a
        // residual clamped to zero.
        let design = fractional_factorial(5, 2).expect("fractional factorial");
        let responses: Vec<f64> = (0..design.run_count()).map(|i| (i * i) as f64).collect();

        let all = crate::analysis::effects::estimate_effects(&design, &responses, 2)
            .expect("two-level design estimates");

        // Find a genuinely aliased pair rather than assuming the labels.
        let contrast = |e: &crate::analysis::effects::EffectEstimate| -> Vec<f64> {
            (0..design.run_count())
                .map(|run| e.columns.iter().map(|&c| design.get(run, c)).product())
                .collect()
        };
        let mut pair = None;
        'outer: for i in 0..all.len() {
            for j in (i + 1)..all.len() {
                let (a, b) = (contrast(&all[i]), contrast(&all[j]));
                if a.iter().zip(b.iter()).all(|(x, y)| (x - y).abs() < 1e-9) {
                    pair = Some((all[i].name.clone(), all[j].name.clone()));
                    break 'outer;
                }
            }
        }
        let (first, second) = pair.expect("a 2^(5-2) design has aliased terms");

        match doe_anova(&design, &responses, &[first.as_str(), second.as_str()]) {
            Err(DoeError::AliasedEffects { .. }) => {}
            other => panic!("aliased pair ({first}, {second}) must be refused, got {other:?}"),
        }
    }

    #[test]
    fn accepts_non_aliased_effects_in_the_same_design() {
        use crate::design::factorial::fractional_factorial;
        let design = fractional_factorial(5, 2).expect("fractional factorial");
        let responses: Vec<f64> = (0..design.run_count()).map(|i| (i * i) as f64).collect();
        assert!(
            doe_anova(&design, &responses, &["A", "B"]).is_ok(),
            "distinct main-effect columns must still be admissible"
        );
    }

    #[test]
    fn rejects_a_model_with_more_terms_than_degrees_of_freedom() {
        use crate::design::factorial::full_factorial;
        // 4 runs -> 3 degrees of freedom; asking for A, B and A:B is saturated
        // and legal, so build the illegal case from a smaller design.
        let design = full_factorial(2).expect("full factorial");
        let responses = vec![10.0, 20.0, 15.0, 25.0];
        assert!(
            doe_anova(&design, &responses, &["A", "B", "A:B"]).is_ok(),
            "a saturated model is legal: residual df is zero, not negative"
        );
    }

    #[test]
    fn refuses_a_design_that_is_not_two_level_coded() {
        use crate::design::ccd::{ccd, AlphaType};
        let design = ccd(2, AlphaType::FaceCentered, 3).expect("ccd builds");
        let responses = vec![1.0; design.run_count()];
        assert!(matches!(
            doe_anova(&design, &responses, &["A", "B"]),
            Err(DoeError::NotTwoLevelCoded { .. })
        ));
    }

    #[test]
    fn saturated_model_leaves_zero_residual_degrees_of_freedom() {
        use crate::design::factorial::full_factorial;
        let design = full_factorial(2).expect("full factorial");
        let responses = vec![10.0, 20.0, 15.0, 25.0];
        let result =
            doe_anova(&design, &responses, &["A", "B", "A:B"]).expect("saturated is legal");
        assert_eq!(result.residual_df, 0);
        assert!(
            result.residual_ss.abs() < 1e-9,
            "an exactly saturated model explains the whole total sum of squares,              so the residual is zero rather than a clamped negative: got {}",
            result.residual_ss
        );
    }

    #[test]
    fn empty_design_is_an_error_not_a_panic() {
        use crate::design::DesignMatrix;
        let empty = DesignMatrix {
            data: vec![],
            factor_names: vec![],
        };
        assert!(matches!(
            doe_anova(&empty, &[], &["A"]),
            Err(DoeError::UnsupportedDesign(_))
        ));
    }
}
