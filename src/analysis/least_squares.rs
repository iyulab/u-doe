//! Least-squares analysis of two-level data whose runs are not balanced.
//!
//! [`crate::analysis::anova::doe_anova`] and
//! [`crate::analysis::effects::estimate_effects`] use the contrast formula,
//! which is exact only when every contrast is balanced and the contrasts are
//! pairwise orthogonal -- and they refuse any input that is not. A run left
//! out, unequal replication of the corners, or a design whose terms are
//! correlated all fail that test, and the experiment would otherwise lose its
//! analysis. This module fits the same model by ordinary least squares and
//! reports each term's **Type III (partial) sum of squares**: the increase in
//! the residual sum of squares when that term alone is dropped from the model.
//!
//! On balanced, orthogonal data the two agree exactly: each least-squares
//! coefficient is half the contrast effect, and the Type III sum of squares
//! is `n · effect² / 4`.
//!
//! Reference: Montgomery, D.C. (2017). *Design and Analysis of Experiments*,
//! 9th ed., §6.5 (regression model of a 2^k) and §10.6 (unbalanced data);
//! Searle, S.R. (1987). *Linear Models for Unbalanced Data*, Ch. 4.

use u_numflow::matrix::Matrix;

use crate::analysis::anova::f_pvalue;
use crate::{design::DesignMatrix, error::DoeError};

/// One term of a least-squares fit.
#[derive(Debug, Clone)]
pub struct LeastSquaresTerm {
    /// Term name, e.g. `"A"` or `"A:B"`.
    pub name: String,
    /// Least-squares coefficient on the ±1 coded contrast.
    pub coefficient: f64,
    /// Effect on the ±1 scale: twice the coefficient.
    pub effect: f64,
    /// Standard error of the coefficient, `sqrt(MS_residual · c_jj)`.
    pub std_error: f64,
    /// `coefficient / std_error`.
    pub t_statistic: f64,
    /// Type III (partial) sum of squares: the residual sum of squares of the
    /// model without this term, minus that of the full model.
    pub sum_of_squares: f64,
    /// Degrees of freedom (always 1 for a single contrast).
    pub df: usize,
    /// `sum_of_squares / df`.
    pub mean_square: f64,
    /// `mean_square / MS_residual`; equals `t_statistic²`.
    pub f_statistic: f64,
    /// p-value of the F test; `None` when there are no residual degrees of
    /// freedom.
    pub p_value: Option<f64>,
}

/// Result of a least-squares fit with Type III sums of squares.
#[derive(Debug, Clone)]
pub struct LeastSquaresFit {
    /// Intercept: the fitted response with every contrast at 0.
    pub intercept: f64,
    /// One entry per requested term, in the order requested.
    pub terms: Vec<LeastSquaresTerm>,
    /// Residual sum of squares.
    pub residual_ss: f64,
    /// Residual degrees of freedom: `n − terms − 1`.
    pub residual_df: usize,
    /// Total sum of squares about the mean.
    pub total_ss: f64,
    /// `1 − residual_ss / total_ss`.
    pub r_squared: f64,
    /// Adjusted R².
    pub r_squared_adj: f64,
    /// Fitted values, one per run.
    pub fitted: Vec<f64>,
    /// Residuals `y − fitted`, one per run.
    pub residuals: Vec<f64>,
}

/// Fit the requested terms to two-level coded data by ordinary least squares
/// and report Type III sums of squares.
///
/// # Arguments
/// * `design`       — Experimental design matrix, every entry coded ±1. The
///   runs need not be balanced: a run may be left out and corners may be
///   replicated unequally.
/// * `responses`    — Observed response values, one per run
/// * `effect_names` — Terms to fit, e.g. `&["A", "B", "A:B"]`. Interaction
///   names join factor names with `":"`; any order is accepted.
///
/// # Balanced data
///
/// When every requested contrast sums to zero and the contrasts are pairwise
/// orthogonal, this returns the same effects, sums of squares, residual, R²
/// and fitted values as [`crate::analysis::anova::doe_anova`].
///
/// # Errors
///
/// [`DoeError::InsufficientResponses`] if `responses.len() != design.run_count()`;
/// [`DoeError::UnknownEffect`] if a name is not a factor or a `:`-joined
/// combination of distinct factors; [`DoeError::NotTwoLevelCoded`] if any
/// entry is not ±1 -- centre points belong in `doe_anova`, which tests them
/// for curvature, and axial points in [`crate::analysis::rsm::fit_rsm`];
/// [`DoeError::OverSpecifiedModel`] if the terms and intercept outnumber the
/// runs; [`DoeError::AliasedEffects`] if the requested columns are linearly
/// dependent (two aliased terms of a regular fraction, or a contrast that is
/// constant over the runs kept), so the coefficients are not identifiable.
///
/// # Examples
///
/// ```
/// use u_doe::design::factorial::full_factorial;
/// use u_doe::analysis::least_squares::fit_least_squares;
///
/// // 2^2 with one run lost: the contrast formula cannot analyse it.
/// let mut design = full_factorial(2).unwrap();
/// design.data.pop();
/// let responses = vec![28.0, 36.0, 18.0];
/// let fit = fit_least_squares(&design, &responses, &["A", "B"]).unwrap();
/// assert_eq!(fit.terms.len(), 2);
/// assert_eq!(fit.residual_df, 0);
/// ```
pub fn fit_least_squares(
    design: &DesignMatrix,
    responses: &[f64],
    effect_names: &[&str],
) -> Result<LeastSquaresFit, DoeError> {
    let n = design.run_count();
    let k = design.factor_count();
    if responses.len() != n {
        return Err(DoeError::InsufficientResponses {
            expected: n,
            got: responses.len(),
        });
    }
    if n == 0 || k == 0 {
        return Err(DoeError::UnsupportedDesign(format!(
            "design has {n} runs and {k} factors; both must be non-zero"
        )));
    }
    if let Some((run, factor, value)) = design.two_level_violation() {
        return Err(DoeError::NotTwoLevelCoded { run, factor, value });
    }

    let terms: Vec<Vec<usize>> = effect_names
        .iter()
        .map(|name| resolve_term(name, &design.factor_names))
        .collect::<Result<_, _>>()?;

    let p = terms.len() + 1;
    if p > n {
        return Err(DoeError::OverSpecifiedModel { terms: p, runs: n });
    }

    // Model matrix: intercept, then one product column per term.
    let rows: Vec<Vec<f64>> = (0..n)
        .map(|run| {
            let mut row = Vec::with_capacity(p);
            row.push(1.0);
            row.extend(
                terms
                    .iter()
                    .map(|cols| cols.iter().map(|&c| design.get(run, c)).product::<f64>()),
            );
            row
        })
        .collect();
    let refs: Vec<&[f64]> = rows.iter().map(Vec::as_slice).collect();
    let x = Matrix::from_rows(&refs);
    let xt = x.transpose();
    let xtx = xt
        .mul_mat(&x)
        .map_err(|e| DoeError::MatrixError(e.to_string()))?;

    // Linear dependence among the columns means the coefficients are not
    // identifiable; name the pair (or the constant column) rather than
    // returning a matrix error.
    if let Some((i, j)) = dependent_columns(&xtx, n) {
        let name = |c: usize| {
            if c == 0 {
                "I".to_string()
            } else {
                effect_names[c - 1].to_string()
            }
        };
        return Err(DoeError::AliasedEffects {
            first: name(i),
            second: name(j),
        });
    }

    let xtx_inv = xtx
        .inverse()
        .map_err(|e| DoeError::MatrixError(e.to_string()))?;
    let xty = xt
        .mul_vec(responses)
        .map_err(|e| DoeError::MatrixError(e.to_string()))?;
    let beta = xtx_inv
        .mul_vec(&xty)
        .map_err(|e| DoeError::MatrixError(e.to_string()))?;

    let fitted: Vec<f64> = rows
        .iter()
        .map(|row| row.iter().zip(&beta).map(|(x, b)| x * b).sum())
        .collect();
    let residuals: Vec<f64> = responses.iter().zip(&fitted).map(|(y, f)| y - f).collect();
    let residual_ss: f64 = residuals.iter().map(|r| r * r).sum();
    let residual_df = n - p;
    let grand_mean = responses.iter().sum::<f64>() / n as f64;
    let total_ss: f64 = responses.iter().map(|&y| (y - grand_mean).powi(2)).sum();
    let ms_residual = if residual_df > 0 {
        residual_ss / residual_df as f64
    } else {
        f64::NAN
    };

    // Type III sum of squares of a single-column term: dropping column j from
    // the model raises the residual sum of squares by b_j² / c_jj, where c_jj
    // is the j-th diagonal entry of (X'X)⁻¹ (Searle 1987, eq. 4.24).
    let terms = effect_names
        .iter()
        .enumerate()
        .map(|(idx, name)| {
            let j = idx + 1;
            let b = beta[j];
            let c = xtx_inv.get(j, j);
            let ss = b * b / c;
            let std_error = (ms_residual * c).sqrt();
            let t = b / std_error;
            let f = if ms_residual > 0.0 && ms_residual.is_finite() {
                ss / ms_residual
            } else {
                f64::NAN
            };
            let p_value = (f.is_finite() && residual_df > 0).then(|| f_pvalue(f, 1, residual_df));
            LeastSquaresTerm {
                name: name.to_string(),
                coefficient: b,
                effect: 2.0 * b,
                std_error,
                t_statistic: t,
                sum_of_squares: ss,
                df: 1,
                mean_square: ss,
                f_statistic: f,
                p_value,
            }
        })
        .collect();

    let r_squared = if total_ss > 1e-12 {
        1.0 - residual_ss / total_ss
    } else {
        0.0
    };
    let total_df = n - 1;
    let r_squared_adj = if residual_df > 0 && total_df > 0 && total_ss > 1e-12 {
        1.0 - (residual_ss / residual_df as f64) / (total_ss / total_df as f64)
    } else {
        r_squared
    };

    Ok(LeastSquaresFit {
        intercept: beta[0],
        terms,
        residual_ss,
        residual_df,
        total_ss,
        r_squared,
        r_squared_adj,
        fitted,
        residuals,
    })
}

/// Resolve `"A:B"` to the factor columns it names, refusing an unknown factor
/// or a factor named twice in one term.
fn resolve_term(name: &str, factor_names: &[String]) -> Result<Vec<usize>, DoeError> {
    let unknown = || DoeError::UnknownEffect {
        name: name.to_string(),
    };
    let mut cols = Vec::new();
    for part in name.split(':') {
        let c = factor_names
            .iter()
            .position(|f| f == part)
            .ok_or_else(unknown)?;
        if cols.contains(&c) {
            return Err(unknown());
        }
        cols.push(c);
    }
    if cols.is_empty() {
        return Err(unknown());
    }
    Ok(cols)
}

/// A pair of model columns that are the same or opposite up to sign over the
/// runs, or a column that is constant (aliased with the intercept). Entries of
/// `X'X` are integers, so `|x_i · x_j| = n` is exact.
///
/// This catches the aliasing a two-level design produces -- coinciding
/// contrasts, and a contrast that the runs kept never vary. Rank deficiency
/// of any other shape is left to the inverse, which reports it as a
/// [`DoeError::MatrixError`].
fn dependent_columns(xtx: &Matrix, n: usize) -> Option<(usize, usize)> {
    let p = xtx.rows();
    for i in 0..p {
        for j in (i + 1)..p {
            if (xtx.get(i, j).abs() - n as f64).abs() < 0.5 {
                return Some((i, j));
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::anova::doe_anova;
    use crate::analysis::effects::estimate_effects;
    use crate::design::factorial::{fractional_factorial, full_factorial};

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9 * (1.0 + a.abs().max(b.abs()))
    }

    /// 2³ × 2 replicates, y = 50 + 5A + 3B + C + 2AB with a fixed noise
    /// pattern, as in the report that motivated this module.
    fn replicated_cube() -> (DesignMatrix, Vec<f64>) {
        let base = full_factorial(3).unwrap();
        let mut data = Vec::new();
        let mut y = Vec::new();
        let noise = [0.3, -0.2, 0.1, 0.4, -0.3, 0.2, -0.1, 0.0];
        for rep in 0..2 {
            for (i, row) in base.data.iter().enumerate() {
                let (a, b, c) = (row[0], row[1], row[2]);
                data.push(row.clone());
                y.push(50.0 + 5.0 * a + 3.0 * b + c + 2.0 * a * b + noise[i] * (rep as f64 + 1.0));
            }
        }
        (
            DesignMatrix {
                data,
                factor_names: base.factor_names.clone(),
            },
            y,
        )
    }

    #[test]
    fn matches_doe_anova_on_balanced_data() {
        let (design, y) = replicated_cube();
        let names = ["A", "B", "C", "A:B"];
        let fit = fit_least_squares(&design, &y, &names).unwrap();
        let anova = doe_anova(&design, &y, &names).unwrap();

        for (t, e) in fit.terms.iter().zip(&anova.effects) {
            assert_eq!(t.name, e.name);
            assert!(
                close(t.sum_of_squares, e.sum_of_squares),
                "{}: {} vs {}",
                t.name,
                t.sum_of_squares,
                e.sum_of_squares
            );
            assert!(close(t.f_statistic, e.f_statistic));
            assert!(close(t.p_value.unwrap(), e.p_value.unwrap()));
        }
        assert!(close(fit.residual_ss, anova.residual_ss));
        assert_eq!(fit.residual_df, anova.residual_df);
        assert!(close(fit.total_ss, anova.total_ss));
        assert!(close(fit.r_squared, anova.r_squared));
        assert!(close(fit.r_squared_adj, anova.r_squared_adj));
        for (a, b) in fit.fitted.iter().zip(&anova.fitted) {
            assert!(close(*a, *b));
        }
        // Effects are the contrast effects, and the coefficients half of them.
        let effects = estimate_effects(&design, &y, 2).unwrap();
        for t in &fit.terms {
            let e = effects.iter().find(|e| e.name == t.name).unwrap();
            assert!(
                close(t.effect, e.estimate),
                "{}: {} vs {}",
                t.name,
                t.effect,
                e.estimate
            );
            assert!(close(t.coefficient, e.estimate / 2.0));
        }
        let grand = y.iter().sum::<f64>() / y.len() as f64;
        assert!(close(fit.intercept, grand));
    }

    #[test]
    fn analyses_the_run_the_contrast_formula_refuses() {
        let (mut design, mut y) = replicated_cube();
        design.data.remove(5);
        y.remove(5);
        let names = ["A", "B", "C", "A:B"];
        assert!(matches!(
            doe_anova(&design, &y, &names),
            Err(DoeError::PartiallyAliasedEffects { .. })
        ));

        let fit = fit_least_squares(&design, &y, &names).unwrap();
        assert_eq!(fit.residual_df, 15 - 5);
        // The model still recovers the generating coefficients to within the
        // noise: |b_A − 5| small, |b_AB − 2| small.
        assert!(
            (fit.terms[0].coefficient - 5.0).abs() < 0.3,
            "{}",
            fit.terms[0].coefficient
        );
        assert!(
            (fit.terms[3].coefficient - 2.0).abs() < 0.3,
            "{}",
            fit.terms[3].coefficient
        );

        // Type III SS is the residual increase from dropping the term,
        // computed independently by refitting without it.
        for (idx, name) in names.iter().enumerate() {
            let reduced: Vec<&str> = names.iter().copied().filter(|n| n != name).collect();
            let without = fit_least_squares(&design, &y, &reduced).unwrap();
            let expected = without.residual_ss - fit.residual_ss;
            assert!(
                close(fit.terms[idx].sum_of_squares, expected),
                "{name}: Type III {} vs refit {expected}",
                fit.terms[idx].sum_of_squares
            );
        }
        // Partial sums of squares of correlated terms do not add up to the
        // model sum of squares -- that is what makes them Type III.
        let model_ss = fit.total_ss - fit.residual_ss;
        let sum: f64 = fit.terms.iter().map(|t| t.sum_of_squares).sum();
        assert!(
            !close(sum, model_ss),
            "unbalanced data should not partition exactly"
        );
        // t² = F.
        for t in &fit.terms {
            assert!(close(t.t_statistic * t.t_statistic, t.f_statistic));
        }
    }

    #[test]
    fn unequal_replication_is_fitted_and_effects_stay_twice_the_coefficients() {
        let base = full_factorial(2).unwrap();
        // Corner (+1,+1) run three times, the others once.
        let mut data = base.data.clone();
        data.push(vec![1.0, 1.0]);
        data.push(vec![1.0, 1.0]);
        let design = DesignMatrix {
            data,
            factor_names: base.factor_names.clone(),
        };
        let y = vec![10.0, 14.0, 12.0, 20.0, 21.0, 19.0];
        let fit = fit_least_squares(&design, &y, &["A", "B", "A:B"]).unwrap();
        assert_eq!(fit.residual_df, 2);
        // Saturated in the four cells: fitted equals the cell means.
        assert!(close(fit.fitted[0], 10.0));
        assert!(close(fit.fitted[3], 20.0));
        assert!(close(fit.fitted[4], 20.0));
        for t in &fit.terms {
            assert!(close(t.effect, 2.0 * t.coefficient));
        }
    }

    #[test]
    fn refuses_aliased_terms_and_names_them() {
        // In the 2^(4-1) with I = ABCD, the columns of D and A:B:C coincide.
        let design = fractional_factorial(4, 1).unwrap();
        let y = vec![1.0, 2.0, 3.0, 5.0, 4.0, 6.0, 7.0, 9.0];
        let err = fit_least_squares(&design, &y, &["D", "A:B:C"]).unwrap_err();
        assert!(
            matches!(&err, DoeError::AliasedEffects { first, second } if first == "D" && second == "A:B:C"),
            "{err:?}"
        );
        // Either alone is fine.
        assert!(fit_least_squares(&design, &y, &["A", "B", "C", "D"]).is_ok());
    }

    #[test]
    fn refuses_a_constant_contrast_over_the_runs_kept() {
        let mut design = full_factorial(2).unwrap();
        // Keep only the runs with A = +1.
        design.data.retain(|row| row[0] > 0.0);
        let err = fit_least_squares(&design, &[1.0, 2.0], &["A"]).unwrap_err();
        assert!(
            matches!(&err, DoeError::AliasedEffects { first, second } if first == "I" && second == "A"),
            "{err:?}"
        );
    }

    #[test]
    fn rejects_unknown_terms_centre_points_and_overspecified_models() {
        let design = full_factorial(2).unwrap();
        let y = vec![1.0, 2.0, 3.0, 4.0];
        assert!(matches!(
            fit_least_squares(&design, &y, &["A", "Z"]),
            Err(DoeError::UnknownEffect { .. })
        ));
        assert!(matches!(
            fit_least_squares(&design, &y, &["A:A"]),
            Err(DoeError::UnknownEffect { .. })
        ));
        assert!(matches!(
            fit_least_squares(&design, &y, &["A", "B", "A:B", "A:B"]),
            Err(DoeError::OverSpecifiedModel { terms: 5, runs: 4 })
        ));
        let mut centred = design.clone();
        centred.data.push(vec![0.0, 0.0]);
        assert!(matches!(
            fit_least_squares(&centred, &[1.0, 2.0, 3.0, 4.0, 2.5], &["A"]),
            Err(DoeError::NotTwoLevelCoded { run: 4, .. })
        ));
        assert!(matches!(
            fit_least_squares(&design, &[1.0, 2.0], &["A"]),
            Err(DoeError::InsufficientResponses {
                expected: 4,
                got: 2
            })
        ));
    }

    #[test]
    fn saturated_fit_has_no_residual_and_no_p_values() {
        let design = full_factorial(2).unwrap();
        let y = vec![1.0, 2.0, 3.0, 7.0];
        let fit = fit_least_squares(&design, &y, &["A", "B", "A:B"]).unwrap();
        assert_eq!(fit.residual_df, 0);
        assert!(fit.residual_ss.abs() < 1e-9);
        assert!(fit.terms.iter().all(|t| t.p_value.is_none()));
        assert!(fit.terms.iter().all(|t| t.f_statistic.is_nan()));
        assert!(close(fit.r_squared, 1.0));
    }
}
