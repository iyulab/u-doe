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

use crate::analysis::anova::{f_pvalue, BlockTerm};
use crate::analysis::has_no_variation;
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
    /// Standard error of the coefficient, `sqrt(MS_residual · c_jj)`. `None`
    /// with no residual degrees of freedom.
    pub std_error: Option<f64>,
    /// `coefficient / std_error`. `None` when `std_error` is `None` or 0.
    pub t_statistic: Option<f64>,
    /// Type III (partial) sum of squares: the residual sum of squares of the
    /// model without this term, minus that of the full model.
    pub sum_of_squares: f64,
    /// Degrees of freedom (always 1 for a single contrast).
    pub df: usize,
    /// `sum_of_squares / df`.
    pub mean_square: f64,
    /// `mean_square / MS_residual`; equals `t_statistic²`. `None` with no
    /// residual degrees of freedom or a residual of exactly zero.
    pub f_statistic: Option<f64>,
    /// p-value of the F test; `None` whenever `f_statistic` is.
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
    /// Residual degrees of freedom: `n − terms − 1`, less the block degrees of
    /// freedom when there are blocks.
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
    /// The block term, when `blocks` were given: its Type III sum of squares
    /// (the rise in the residual sum of squares when the block columns are
    /// dropped) on `blocks − 1` degrees of freedom; `None` otherwise. As in
    /// [`doe_anova`](crate::analysis::anova::doe_anova), no F test.
    pub blocks: Option<BlockTerm>,
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
/// * `blocks`       — The block of each run (any labels), or `None`.
///
/// # Blocks
///
/// A block is a fixed effect: one column per block but the first, coded
/// `+1` in that block, `−1` in the first and 0 elsewhere (sum-to-zero), so the
/// intercept is the average of the block means at the centre of the design.
/// Unlike [`doe_anova`](crate::analysis::anova::doe_anova), a term need not
/// be orthogonal to the blocks -- an unequal number of runs per block, or a
/// run lost from one, is what least squares is for; the Type III sums of
/// squares are each adjusted for the blocks. Only a term that is constant
/// within every block is refused: it is confounded with them.
///
/// # Balanced data
///
/// When every requested contrast sums to zero and the contrasts are pairwise
/// orthogonal, this returns the same effects, sums of squares, residual, R²
/// and fitted values as [`crate::analysis::anova::doe_anova`].
///
/// # Errors
///
/// [`DoeError::DimensionMismatch`] if `responses.len() != design.run_count()`;
/// [`DoeError::UnknownEffect`] if a name is not a factor or a `:`-joined
/// combination of distinct factors; [`DoeError::NotTwoLevelCoded`] if any
/// entry is not ±1 -- centre points belong in `doe_anova`, which tests them
/// for curvature, and axial points in [`crate::analysis::rsm::fit_rsm`];
/// [`DoeError::OverSpecifiedModel`] if the terms and intercept outnumber the
/// runs; [`DoeError::AliasedEffects`] if the requested columns are linearly
/// dependent (two aliased terms of a regular fraction, or a contrast that is
/// constant over the runs kept), so the coefficients are not identifiable --
/// with `Block` as `second` for a term constant within every block;
/// [`DoeError::DimensionMismatch`] (`parameter` = `blocks`) if `blocks` does
/// not hold one block per run.
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
/// let fit = fit_least_squares(&design, &responses, &["A", "B"], None).unwrap();
/// assert_eq!(fit.terms.len(), 2);
/// assert_eq!(fit.residual_df, 0);
/// ```
pub fn fit_least_squares(
    design: &DesignMatrix,
    responses: &[f64],
    effect_names: &[&str],
    blocks: Option<&[usize]>,
) -> Result<LeastSquaresFit, DoeError> {
    let n = design.run_count();
    if responses.len() != n {
        return Err(DoeError::DimensionMismatch {
            parameter: "responses",
            index: None,
            expected: n,
            got: responses.len(),
            counted: "one response per run",
        });
    }
    if let Some(blocks) = blocks {
        if blocks.len() != n {
            return Err(DoeError::DimensionMismatch {
                parameter: "blocks",
                index: None,
                expected: n,
                got: blocks.len(),
                counted: "one block per run",
            });
        }
    }
    design.check_not_empty()?;
    design.check_shape()?;
    if let Some((run, factor, value)) = design.two_level_violation() {
        return Err(DoeError::NotTwoLevelCoded { run, factor, value });
    }

    let terms: Vec<Vec<usize>> = effect_names
        .iter()
        .map(|name| resolve_term(name, &design.factor_names))
        .collect::<Result<_, _>>()?;

    // The block of each run, numbered 0.. in label order.
    let block_of: Vec<usize> = match blocks {
        Some(labels) => {
            let mut sorted: Vec<usize> = labels.to_vec();
            sorted.sort_unstable();
            sorted.dedup();
            labels
                .iter()
                .map(|l| sorted.binary_search(l).expect("every label is in the list"))
                .collect()
        }
        None => vec![0; n],
    };
    let n_blocks = block_of.iter().max().map_or(1, |&b| b + 1);
    let block_df = n_blocks - 1;

    // One degree of freedom goes to the mean, one to each term and
    // `blocks − 1` to the blocks.
    let k = terms.len() + 1;
    let p = k + block_df;
    if p > n {
        return Err(DoeError::OverSpecifiedModel {
            terms: terms.len() + block_df,
            runs: n,
        });
    }

    // Model matrix: intercept, one product column per term, then the block
    // columns (sum-to-zero against the first block) -- last, so a term keeps
    // its column index whether or not there are blocks.
    let rows: Vec<Vec<f64>> = (0..n)
        .map(|run| {
            let mut row = Vec::with_capacity(p);
            row.push(1.0);
            row.extend(
                terms
                    .iter()
                    .map(|cols| cols.iter().map(|&c| design.get(run, c)).product::<f64>()),
            );
            row.extend((1..n_blocks).map(|b| match block_of[run] {
                0 => -1.0,
                own if own == b => 1.0,
                _ => 0.0,
            }));
            row
        })
        .collect();

    // Linear dependence among the term columns means the coefficients are not
    // identifiable; name the pair (or the constant column) rather than
    // returning a matrix error. Checked on the unblocked columns: a two-block
    // column is ±1 like a contrast and would pass for one.
    let term_rows: Vec<&[f64]> = rows.iter().map(|r| &r[..k]).collect();
    let term_x = Matrix::from_rows(&term_rows);
    let term_xtx = term_x
        .transpose()
        .mul_mat(&term_x)
        .expect("X' is k x n and X is n x k by construction");
    if let Some((i, j)) = dependent_columns(&term_xtx, n) {
        // Column 0 is the mean. It goes second, as `I`, the place every other
        // path in the crate puts it -- a caller reads `second == "I"` as "this
        // term is constant over the runs".
        let term = |c: usize| effect_names[c - 1].to_string();
        let (first, second) = if i == 0 {
            (term(j), "I".to_string())
        } else {
            (term(i), term(j))
        };
        return Err(DoeError::AliasedEffects { first, second });
    }
    // A term constant within every block lies in the span of the blocks: it is
    // confounded with them.
    if n_blocks > 1 {
        for (idx, name) in effect_names.iter().enumerate() {
            let mut level = vec![None::<f64>; n_blocks];
            let constant = rows.iter().zip(&block_of).all(|(row, &b)| {
                let v = row[idx + 1];
                match level[b] {
                    None => {
                        level[b] = Some(v);
                        true
                    }
                    Some(seen) => (seen - v).abs() < 1e-9,
                }
            });
            if constant {
                return Err(DoeError::AliasedEffects {
                    first: name.to_string(),
                    second: "Block".to_string(),
                });
            }
        }
    }

    let refs: Vec<&[f64]> = rows.iter().map(Vec::as_slice).collect();
    let x = Matrix::from_rows(&refs);
    let xt = x.transpose();
    let xtx = xt
        .mul_mat(&x)
        .expect("X' is p x n and X is n x p by construction");
    let xtx_inv = xtx.inverse().map_err(|_| DoeError::SingularModel)?;
    let xty = xt
        .mul_vec(responses)
        .expect("X' has n columns and there are n responses, checked above");
    let beta = xtx_inv
        .mul_vec(&xty)
        .expect("(X'X)^-1 is p x p and X'y has p entries by construction");

    let fitted: Vec<f64> = rows
        .iter()
        .map(|row| row.iter().zip(&beta).map(|(x, b)| x * b).sum())
        .collect();
    let residuals: Vec<f64> = responses.iter().zip(&fitted).map(|(y, f)| y - f).collect();
    let residual_ss: f64 = residuals.iter().map(|r| r * r).sum();
    let residual_df = n - p;
    let grand_mean = responses.iter().sum::<f64>() / n as f64;
    let total_ss: f64 = responses.iter().map(|&y| (y - grand_mean).powi(2)).sum();
    let ms_residual = (residual_df > 0).then(|| residual_ss / residual_df as f64);

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
            let std_error = ms_residual.map(|ms| (ms * c).sqrt());
            let t = std_error.filter(|&se| se > 0.0).map(|se| b / se);
            let f = ms_residual.filter(|&ms| ms > 0.0).map(|ms| ss / ms);
            let p_value = f.map(|f| f_pvalue(f, 1, residual_df));
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

    // The blocks' Type III sum of squares: dropping their columns raises the
    // residual sum of squares by b_B' (C_BB)⁻¹ b_B, C_BB the block block of
    // (X'X)⁻¹ (Searle 1987, eq. 4.24 for several columns at once).
    let block_term = match blocks {
        None => None,
        Some(_) if block_df == 0 => Some(BlockTerm {
            sum_of_squares: 0.0,
            df: 0,
            mean_square: 0.0,
        }),
        Some(_) => {
            let c_bb = Matrix::from_rows(
                &(k..p)
                    .map(|i| (k..p).map(|j| xtx_inv.get(i, j)).collect::<Vec<f64>>())
                    .collect::<Vec<_>>()
                    .iter()
                    .map(Vec::as_slice)
                    .collect::<Vec<_>>(),
            );
            let c_bb_inv = c_bb.inverse().map_err(|_| DoeError::SingularModel)?;
            let b_b = &beta[k..p];
            let ss = c_bb_inv
                .mul_vec(b_b)
                .expect("(C_BB)^-1 is square in the block columns")
                .iter()
                .zip(b_b)
                .map(|(a, b)| a * b)
                .sum::<f64>()
                .max(0.0);
            Some(BlockTerm {
                sum_of_squares: ss,
                df: block_df,
                mean_square: ss / block_df as f64,
            })
        }
    };

    let r_squared = if !has_no_variation(total_ss, responses) {
        1.0 - residual_ss / total_ss
    } else {
        0.0
    };
    let total_df = n - 1;
    let r_squared_adj = if residual_df > 0 && total_df > 0 && !has_no_variation(total_ss, responses)
    {
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
        blocks: block_term,
    })
}

/// Resolve `"A:B"` to the factor columns it names, refusing an unknown factor
/// or a factor named twice in one term.
fn resolve_term(name: &str, factor_names: &[String]) -> Result<Vec<usize>, DoeError> {
    let unknown = || DoeError::UnknownEffect {
        effect: name.to_string(),
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

    /// A 2^3 with the four runs at C = -1 left out: C is constant over the runs
    /// kept. The mean goes second, so the message says what happened.
    #[test]
    fn a_factor_with_one_level_left_is_aliased_with_the_mean() {
        let full = full_factorial(3).unwrap();
        let kept: Vec<Vec<f64>> = full.data.into_iter().filter(|r| r[2] > 0.0).collect();
        let y: Vec<f64> = (0..kept.len()).map(|i| 50.0 + 0.3 * i as f64).collect();
        let design = DesignMatrix {
            data: kept,
            factor_names: full.factor_names,
        };
        let err = fit_least_squares(&design, &y, &["A", "B", "C"], None).map(|f| f.intercept);
        assert_eq!(
            err,
            Err(DoeError::AliasedEffects {
                first: "C".into(),
                second: "I".into(),
            })
        );
        let message = err.expect_err("refused").to_string();
        assert!(
            message.contains("aliased with the overall mean"),
            "{message}"
        );
    }

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
        let fit = fit_least_squares(&design, &y, &names, None).unwrap();
        let anova = doe_anova(&design, &y, &names, None).unwrap();

        for (t, e) in fit.terms.iter().zip(&anova.effects) {
            assert_eq!(t.name, e.name);
            assert!(
                close(t.sum_of_squares, e.sum_of_squares),
                "{}: {} vs {}",
                t.name,
                t.sum_of_squares,
                e.sum_of_squares
            );
            assert!(close(t.f_statistic.unwrap(), e.f_statistic.unwrap()));
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
            doe_anova(&design, &y, &names, None),
            Err(DoeError::PartiallyAliasedEffects { .. })
        ));

        let fit = fit_least_squares(&design, &y, &names, None).unwrap();
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
            let without = fit_least_squares(&design, &y, &reduced, None).unwrap();
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
            let (t_stat, f) = (t.t_statistic.unwrap(), t.f_statistic.unwrap());
            assert!(close(t_stat * t_stat, f));
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
        let fit = fit_least_squares(&design, &y, &["A", "B", "A:B"], None).unwrap();
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
        let err = fit_least_squares(&design, &y, &["D", "A:B:C"], None).unwrap_err();
        assert!(
            matches!(&err, DoeError::AliasedEffects { first, second } if first == "D" && second == "A:B:C"),
            "{err:?}"
        );
        // Either alone is fine.
        assert!(fit_least_squares(&design, &y, &["A", "B", "C", "D"], None).is_ok());
    }

    #[test]
    fn rejects_unknown_terms_centre_points_and_overspecified_models() {
        let design = full_factorial(2).unwrap();
        let y = vec![1.0, 2.0, 3.0, 4.0];
        assert!(matches!(
            fit_least_squares(&design, &y, &["A", "Z"], None),
            Err(DoeError::UnknownEffect { .. })
        ));
        assert!(matches!(
            fit_least_squares(&design, &y, &["A:A"], None),
            Err(DoeError::UnknownEffect { .. })
        ));
        assert!(matches!(
            fit_least_squares(&design, &y, &["A", "B", "A:B", "A:B"], None),
            // Four terms need four degrees of freedom; four runs leave three
            // after the mean. `terms` counts the terms, not the mean.
            Err(DoeError::OverSpecifiedModel { terms: 4, runs: 4 })
        ));
        let mut centred = design.clone();
        centred.data.push(vec![0.0, 0.0]);
        assert!(matches!(
            fit_least_squares(&centred, &[1.0, 2.0, 3.0, 4.0, 2.5], &["A"], None),
            Err(DoeError::NotTwoLevelCoded { run: 4, .. })
        ));
        assert!(matches!(
            fit_least_squares(&design, &[1.0, 2.0], &["A"], None),
            Err(DoeError::DimensionMismatch {
                parameter: "responses",
                expected: 4,
                got: 2,
                ..
            })
        ));
    }

    #[test]
    fn saturated_fit_has_no_residual_and_no_p_values() {
        let design = full_factorial(2).unwrap();
        let y = vec![1.0, 2.0, 3.0, 7.0];
        let fit = fit_least_squares(&design, &y, &["A", "B", "A:B"], None).unwrap();
        assert_eq!(fit.residual_df, 0);
        assert!(fit.residual_ss.abs() < 1e-9);
        assert!(fit.terms.iter().all(|t| t.p_value.is_none()));
        assert!(fit
            .terms
            .iter()
            .all(|t| t.f_statistic.is_none() && t.std_error.is_none() && t.t_statistic.is_none()));
        assert!(close(fit.r_squared, 1.0));
    }
}

/// Type III least squares against statsmodels' OLS and `anova_lm(typ=3)`, on
/// an unbalanced two-level design: eight corners less two, three of the rest
/// replicated, the model A, B, C and A:B. The responses are drawn from a known
/// model with a seeded NumPy generator.
#[cfg(test)]
mod statsmodels_reference {
    use super::*;

    fn close(got: f64, want: f64) -> bool {
        (got - want).abs() <= 1e-9 * want.abs().max(1.0)
    }

    #[test]
    fn an_unbalanced_fit_reproduces_statsmodels() {
        let design = DesignMatrix {
            data: vec![
                vec![-1.0, -1.0, -1.0],
                vec![1.0, -1.0, -1.0],
                vec![-1.0, 1.0, -1.0],
                vec![-1.0, -1.0, 1.0],
                vec![1.0, -1.0, 1.0],
                vec![1.0, 1.0, 1.0],
                vec![-1.0, -1.0, -1.0],
                vec![-1.0, -1.0, -1.0],
                vec![1.0, -1.0, 1.0],
            ],
            factor_names: vec!["A".into(), "B".into(), "C".into()],
        };
        let responses = [
            9.881, 11.9, 4.167, 10.474, 13.549, 12.328, 9.01, 8.853, 13.484,
        ];
        let fit = fit_least_squares(&design, &responses, &["A", "B", "C", "A:B"], None)
            .expect("a valid fit");
        assert!(
            close(fit.intercept, 9.786161764705883),
            "intercept {}",
            fit.intercept
        );
        // name, coefficient, standard error, t, Type III SS, F, p
        let terms: [(&str, f64, f64, f64, f64, f64, f64); 4] = [
            (
                "A",
                2.39674999999999,
                0.204525455053169,
                11.7185902330687,
                22.9776422499999,
                137.325357050572,
                0.000303287891821598,
            ),
            (
                "B",
                -1.53866176470588,
                0.164520207843679,
                -9.3524180699301,
                14.6353310708556,
                87.467723754755,
                0.00072787808654291,
            ),
            (
                "C",
                0.704882352941182,
                0.171835752832631,
                4.1020703859443,
                2.81553507843141,
                16.8269814512412,
                0.0148279126984472,
            ),
            (
                "A:B",
                0.978867647058818,
                0.171835752832631,
                5.69653073311374,
                5.42969726593133,
                32.4504623933093,
                0.00469210566782265,
            ),
        ];
        for (t, (name, coefficient, se, tstat, ss, f, p)) in fit.terms.iter().zip(terms) {
            assert_eq!(t.name, name);
            assert!(
                close(t.coefficient, coefficient),
                "{name} coefficient {}",
                t.coefficient
            );
            assert!(close(t.effect, 2.0 * coefficient), "{name} effect");
            assert!(
                close(t.std_error.expect("residual df"), se),
                "{name} std error"
            );
            assert!(
                close(t.t_statistic.expect("residual df"), tstat),
                "{name} t"
            );
            assert!(
                close(t.sum_of_squares, ss),
                "{name} Type III SS {}",
                t.sum_of_squares
            );
            assert!(close(t.f_statistic.expect("residual df"), f), "{name} F");
            assert!(
                close(t.p_value.expect("residual df"), p),
                "{name} p {}",
                t.p_value.unwrap()
            );
        }
        assert!(close(fit.residual_ss, 0.669290588235294));
        assert_eq!(fit.residual_df, 4);
        assert!(close(fit.total_ss, 68.8446408888889));
        assert!(close(fit.r_squared, 0.990278247085122));
        assert!(close(fit.r_squared_adj, 0.980556494170244));
        for (got, want) in fit.fitted.iter().zip([
            9.20205882352941,
            12.0378235294118,
            4.16700000000001,
            10.6118235294118,
            13.4475882352941,
            12.328,
            9.20205882352941,
            9.20205882352941,
            13.4475882352941,
        ]) {
            assert!(close(*got, want), "fitted {got} vs {want}");
        }
    }
}

/// Blocks in the least-squares fit: on balanced blocked data the same table
/// as `doe_anova`; on unbalanced blocks (unequal days, a run lost) the Type
/// III table statsmodels 0.15.0 computes (`ols("y ~ A + B + C + A:B +
/// C(block, Sum)")`, `anova_lm(typ=3)`, seeded NumPy noise).
#[cfg(test)]
mod blocks {
    use super::*;
    use crate::analysis::anova::doe_anova;
    use crate::design::factorial::blocked_factorial;

    fn close(got: f64, want: f64) -> bool {
        (got - want).abs() <= 1e-9 * want.abs().max(1.0)
    }

    #[test]
    fn balanced_blocks_give_the_doe_anova_table() {
        let blocked = blocked_factorial(4, 4).expect("2^4 in 4 blocks");
        let y: Vec<f64> = blocked
            .design
            .data
            .iter()
            .zip(&blocked.blocks)
            .enumerate()
            .map(|(run, (r, &b))| {
                20.0 + 3.0 * r[0] - 2.0 * r[1]
                    + r[0] * r[2]
                    + 1.5 * b as f64
                    + ((run * 7 % 5) as f64 - 2.0) * 0.2
            })
            .collect();
        let names = ["A", "B", "C", "A:C"];
        let a = doe_anova(&blocked.design, &y, &names, Some(&blocked.blocks)).expect("anova");
        let f = fit_least_squares(&blocked.design, &y, &names, Some(&blocked.blocks))
            .expect("least squares");
        for (row, term) in a.effects.iter().zip(&f.terms) {
            assert!(
                close(term.sum_of_squares, row.sum_of_squares),
                "{}",
                row.name
            );
            assert!(close(
                term.f_statistic.expect("F"),
                row.f_statistic.expect("F")
            ));
        }
        let (ab, fb) = (a.blocks.expect("blocks"), f.blocks.expect("blocks"));
        assert_eq!(ab.df, fb.df);
        assert!(
            close(fb.sum_of_squares, ab.sum_of_squares),
            "{} vs {}",
            fb.sum_of_squares,
            ab.sum_of_squares
        );
        assert!(close(f.residual_ss, a.residual_ss));
        assert_eq!(f.residual_df, a.residual_df);
        assert!(close(f.r_squared, a.r_squared));
        for (x, y) in f.fitted.iter().zip(&a.fitted) {
            assert!(close(*x, *y), "{x} vs {y}");
        }
    }

    const ROWS: [[f64; 4]; 15] = [
        [-1.0, -1.0, -1.0, -1.0],
        [1.0, -1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0, -1.0],
        [1.0, 1.0, -1.0, -1.0],
        [-1.0, -1.0, 1.0, -1.0],
        [-1.0, 1.0, 1.0, -1.0],
        [1.0, 1.0, 1.0, -1.0],
        [-1.0, -1.0, -1.0, 1.0],
        [1.0, -1.0, -1.0, 1.0],
        [-1.0, 1.0, -1.0, 1.0],
        [1.0, 1.0, -1.0, 1.0],
        [-1.0, -1.0, 1.0, 1.0],
        [1.0, -1.0, 1.0, 1.0],
        [-1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0],
    ];
    const BLOCKS: [usize; 15] = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2];
    const Y: [f64; 15] = [
        27.98022, 32.410836, 28.344161, 36.913045, 23.644812, 30.508867, 39.026795, 33.313916,
        36.377557, 33.276572, 33.908528, 21.427257, 25.535704, 23.796101, 31.486565,
    ];
    const FITTED: [f64; 15] = [
        27.584515475748216,
        31.686167768317855,
        28.708167917956676,
        36.967157960268324,
        24.347064877709,
        30.5356035510836,
        38.794593593395255,
        32.649401706914354,
        36.75105399948401,
        33.77305414912282,
        34.405316022187826,
        21.785222939628493,
        25.88687523219814,
        22.908875381836953,
        31.167865424148612,
    ];

    #[test]
    fn unbalanced_blocks_reproduce_statsmodels_type_iii() {
        let design = DesignMatrix {
            data: ROWS.iter().map(|r| r.to_vec()).collect(),
            factor_names: ["A", "B", "C", "D"].map(String::from).to_vec(),
        };
        // Labels other than 0..b: only their equality matters.
        let labels: Vec<usize> = BLOCKS.iter().map(|b| 10 * b + 7).collect();
        let fit = fit_least_squares(&design, &Y, &["A", "B", "C", "A:B"], Some(&labels))
            .expect("a valid fit");
        // name, coefficient, SS, F, p
        let terms: [(&str, f64, f64, f64, f64); 4] = [
            (
                "A",
                3.0901605837203254,
                131.48235265009174,
                293.8548367330318,
                1.3635848133201175e-07,
            ),
            (
                "B",
                1.601160658539733,
                36.667753141952346,
                81.95013547993433,
                1.7753229332507856e-05,
            ),
            (
                "C",
                -1.6187252990196095,
                28.133442374314136,
                62.876484554875894,
                4.653990554576451e-05,
            ),
            (
                "A:B",
                1.0393344374355027,
                15.364834856240028,
                34.33944515828725,
                0.00037847158463040774,
            ),
        ];
        for (term, (name, coef, ss, f, p)) in fit.terms.iter().zip(terms) {
            assert_eq!(term.name, name);
            assert!(
                close(term.coefficient, coef),
                "{name} b {}",
                term.coefficient
            );
            assert!(
                close(term.sum_of_squares, ss),
                "{name} SS {}",
                term.sum_of_squares
            );
            assert!(close(term.f_statistic.expect("F"), f), "{name} F");
            assert!(close(term.p_value.expect("p"), p), "{name} p");
        }
        assert!(
            close(fit.intercept, 30.452125079248376),
            "{}",
            fit.intercept
        );
        let block = fit.blocks.expect("blocks");
        assert_eq!(block.df, 2);
        assert!(
            close(block.sum_of_squares, 133.42524565527805),
            "{}",
            block.sum_of_squares
        );
        assert!(close(fit.residual_ss, 3.5795184891115186));
        assert_eq!(fit.residual_df, 8);
        assert!(close(fit.r_squared, 0.9908247451533986));
        assert!(close(fit.r_squared_adj, 0.9839433040184475));
        for (run, (got, want)) in fit.fitted.iter().zip(FITTED).enumerate() {
            assert!(close(*got, want), "fitted[{run}] {got} vs {want}");
        }
        // The block term by its definition: the rise in the residual when the
        // block columns are dropped.
        let without = fit_least_squares(&design, &Y, &["A", "B", "C", "A:B"], None).expect("fit");
        assert!(close(
            without.residual_ss - fit.residual_ss,
            block.sum_of_squares
        ));
    }

    #[test]
    fn a_term_constant_within_every_block_is_confounded() {
        let design = DesignMatrix {
            data: ROWS.iter().map(|r| r.to_vec()).collect(),
            factor_names: ["A", "B", "C", "D"].map(String::from).to_vec(),
        };
        // Two blocks by D: D is the block difference.
        let by_d: Vec<usize> = ROWS.iter().map(|r| usize::from(r[3] > 0.0)).collect();
        assert_eq!(
            fit_least_squares(&design, &Y, &["A", "D"], Some(&by_d)).map(|f| f.intercept),
            Err(DoeError::AliasedEffects {
                first: "D".into(),
                second: "Block".into()
            })
        );
        assert!(matches!(
            fit_least_squares(&design, &Y, &["A"], Some(&[0, 1])),
            Err(DoeError::DimensionMismatch {
                parameter: "blocks",
                ..
            })
        ));
        // 15 runs, 4 terms + mean + 14 blocks is too many.
        let each: Vec<usize> = (0..15).collect();
        assert!(matches!(
            fit_least_squares(&design, &Y, &["A", "B", "C", "A:B"], Some(&each)),
            Err(DoeError::OverSpecifiedModel {
                terms: 18,
                runs: 15
            })
        ));
    }
}
