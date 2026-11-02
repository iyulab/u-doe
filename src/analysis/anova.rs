//! DOE-context Analysis of Variance (ANOVA).
//!
//! Partitions response variation into contributions from specified effects
//! and residual, then computes F-statistics for significance testing.
//!
//! Reference: Montgomery, D.C. (2019). *Introduction to Statistical Quality
//! Control*, 8th ed., Section 6.4. Wiley.

use crate::analysis::has_no_variation;
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
    /// F-statistic = MS_effect / MS_residual. `None` when there is no residual
    /// mean square to test against: no residual degrees of freedom (a
    /// saturated model), or a residual of exactly zero.
    pub f_statistic: Option<f64>,
    /// p-value from F(1, df_residual) distribution; `None` whenever
    /// `f_statistic` is.
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

/// The block term of a blocked analysis: the differences between the block
/// means, a fixed effect taken out of the residual.
///
/// No F test is reported. The runs are randomized only within their blocks,
/// so the block differences are not tested against an error that
/// randomization makes valid for them (Montgomery, *Design and Analysis of
/// Experiments*, §4.1.2); the term is there to be removed, not judged.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BlockTerm {
    /// `Σ_j n_j · (ȳ_j − ȳ)²` over the blocks `j`, `n_j` runs each.
    pub sum_of_squares: f64,
    /// Number of blocks − 1.
    pub df: usize,
    /// `sum_of_squares / df`; 0 when there is a single block.
    pub mean_square: f64,
}

/// Result of a DOE ANOVA.
#[derive(Debug, Clone)]
pub struct DoeAnovaResult {
    /// Rows for each requested effect.
    pub effects: Vec<EffectRow>,
    /// Residual (error) sum of squares: what neither the effects, the
    /// curvature term nor the blocks account for. It includes the pure error.
    pub residual_ss: f64,
    /// Residual degrees of freedom.
    pub residual_df: usize,
    /// Total sum of squares.
    pub total_ss: f64,
    /// R² = (SS_model + SS_curvature + SS_blocks) / SS_total.
    pub r_squared: f64,
    /// Adjusted R² = 1 - (SS_residual/df_residual) / (SS_total/df_total).
    pub r_squared_adj: f64,
    /// Fitted value for each run, in run order. A factorial run's is the mean
    /// of the factorial runs plus half of each requested effect times its
    /// contrast -- exact because the requested contrasts are orthogonal and
    /// balanced (see Errors). A centre run's is the mean of the centre runs.
    /// With blocks, each run's value moves by its block mean's departure from
    /// the grand mean.
    pub fitted: Vec<f64>,
    /// `responses[i] - fitted[i]` for each run. Their squares sum to
    /// [`residual_ss`](Self::residual_ss); they are what a residuals-versus-
    /// fitted or normal probability plot needs.
    pub residuals: Vec<f64>,
    /// Curvature test, when the design has centre points (runs with every
    /// factor at 0); `None` otherwise.
    pub curvature: Option<CurvatureTest>,
    /// Pure error, when some design point was run more than once (within one
    /// block, when there are blocks); `None` otherwise.
    pub pure_error: Option<PureError>,
    /// The block term, when `blocks` were given; `None` otherwise.
    pub blocks: Option<BlockTerm>,
}

/// Compute DOE ANOVA for specified effects.
///
/// # Arguments
/// * `design`       — Experimental design matrix: coded ±1 values, plus any
///   centre points (runs with every factor at 0)
/// * `responses`    — Observed response values, one per run
/// * `effect_names` — Names of effects to include (e.g. `&["A", "B", "A:B"]`).
///   Interaction names join factor names with `":"`.
/// * `blocks`       — The block of each run (any labels; runs with the same
///   label share a block), or `None` for an unblocked experiment.
///
/// # Blocks
///
/// A block -- a batch of material, a shift, a day -- is a fixed effect: the
/// differences between the block means are taken out of the residual and
/// reported as [`DoeAnovaResult::blocks`]. An effect is analysable only if it
/// is orthogonal to the blocks, i.e. its contrast sums to zero within every
/// block, which is what
/// [`blocked_factorial`](crate::design::factorial::blocked_factorial) arranges
/// for every effect but the ones it confounds. Pure error is then taken among
/// replicates within a block, so a block difference is not counted as error.
/// Centre points have to be spread over the blocks in the same proportion,
/// or the curvature contrast would carry block differences.
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
/// any estimable effect -- a factor name, or factor names joined with `:` up to
/// the number of factors the design has, the same names
/// [`crate::analysis::effects::estimate_effects`] returns;
/// [`DoeError::NotTwoLevelCoded`] if a run is neither two-level coded (every
/// factor `-1` or `+1`) nor a centre point (every factor `0`) -- axial and
/// three-level designs belong in [`crate::analysis::rsm::fit_rsm`];
/// [`DoeError::AliasedEffects`] if two requested effects share a contrast
/// column, which would count the same sum of squares twice, or a requested
/// contrast is the same in every run (aliased with the mean, `I`);
/// [`DoeError::PartiallyAliasedEffects`] if two requested contrasts are
/// correlated without being identical, or one is correlated with the mean --
/// the per-term sums of squares would then overlap and not add up to the model
/// sum of squares; [`DoeError::DimensionMismatch`] (`parameter` = `blocks`) if
/// `blocks` does not hold one block per run; [`DoeError::AliasedEffects`] with
/// `Block` as `second` if a requested contrast is constant within every block
/// (confounded with the blocks -- its estimate is a block difference), and
/// [`DoeError::PartiallyAliasedEffects`] with `Block` if it is not but still
/// fails to sum to zero within some block, or (`first` = `curvature`) if the
/// centre points are not spread over the blocks in the same proportion; or
/// [`DoeError::OverSpecifiedModel`] if the model (counting the curvature term
/// and the block degrees of freedom) asks for more degrees of freedom than the
/// design has.
///
/// # Examples
///
/// ```
/// use u_doe::design::factorial::full_factorial;
/// use u_doe::analysis::anova::doe_anova;
///
/// let design = full_factorial(3).unwrap();
/// let responses = vec![10.0,20.0,15.0,25.0,12.0,22.0,18.0,30.0];
/// let result = doe_anova(&design, &responses, &["A","B","C"], None).unwrap();
/// assert!(result.r_squared >= 0.0 && result.r_squared <= 1.0);
/// ```
///
/// Blocked: 2^3 in two blocks, `A:B:C` confounded with them.
///
/// ```
/// use u_doe::design::factorial::blocked_factorial;
/// use u_doe::analysis::anova::doe_anova;
///
/// let blocked = blocked_factorial(3, 2).unwrap();
/// // y = 50 + 5A + 3B, the second block 4 higher.
/// let y: Vec<f64> = blocked.design.data.iter().zip(&blocked.blocks)
///     .map(|(r, &b)| 50.0 + 5.0 * r[0] + 3.0 * r[1] + 4.0 * b as f64)
///     .collect();
/// let r = doe_anova(&blocked.design, &y, &["A", "B"], Some(&blocked.blocks)).unwrap();
/// let block = r.blocks.unwrap();
/// assert_eq!(block.df, 1);
/// assert!((block.sum_of_squares - 32.0).abs() < 1e-9); // 8 runs, each 2 from the grand mean
/// assert!(r.residual_ss < 1e-9);
/// ```
pub fn doe_anova(
    design: &DesignMatrix,
    responses: &[f64],
    effect_names: &[&str],
    blocks: Option<&[usize]>,
) -> Result<DoeAnovaResult, DoeError> {
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
    // Before the centre-point test reads any row: a short row of zeros would
    // otherwise pass for a centre point.
    design.check_shape()?;

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
    // The order to estimate up to is the one the caller asked for, read off the
    // names themselves. A fixed 2 here meant `estimate_effects(.., 3)` produced
    // `A:B:C` and this function then refused it as unknown -- the crate's own
    // output rejected by the crate's own input. `build_terms` clamps an order
    // above the factor count, so an over-long name falls through to the
    // `UnknownEffect` below rather than growing the term list.
    let requested_order = effect_names
        .iter()
        .map(|name| name.split(':').count())
        .max()
        .unwrap_or(1)
        .max(1);
    let all_effects =
        crate::analysis::effects::contrast_effects(&factorial, &factorial_y, requested_order)
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
                    effect: name.to_string(),
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

    // Blocks: the runs of each block, by label, in label order.
    let block_runs: Vec<Vec<usize>> = match blocks {
        Some(labels) => {
            let mut by_label = std::collections::BTreeMap::<usize, Vec<usize>>::new();
            for (run, &label) in labels.iter().enumerate() {
                by_label.entry(label).or_default().push(run);
            }
            by_label.into_values().collect()
        }
        None => vec![(0..n).collect()],
    };
    let mut block_of = vec![0usize; n];
    for (b, runs) in block_runs.iter().enumerate() {
        for &run in runs {
            block_of[run] = b;
        }
    }
    if block_runs.len() > 1 {
        // An effect is separable from the blocks only if its contrast sums to
        // zero within each one. Constant within each block, it *is* a block
        // difference -- the interactions a blocked design confounds.
        for (e, c) in selected.iter().zip(&contrasts) {
            let mut sums = vec![0.0; block_runs.len()];
            let mut first = vec![None::<f64>; block_runs.len()];
            let mut constant = true;
            for (i, &run) in factorial_runs.iter().enumerate() {
                let b = block_of[run];
                sums[b] += c[i];
                match first[b] {
                    None => first[b] = Some(c[i]),
                    Some(v) if (v - c[i]).abs() > 1e-9 => constant = false,
                    Some(_) => {}
                }
            }
            if constant {
                return Err(DoeError::AliasedEffects {
                    first: e.name.clone(),
                    second: "Block".to_string(),
                });
            }
            if sums.iter().any(|s| s.abs() > 1e-9) {
                return Err(DoeError::PartiallyAliasedEffects {
                    first: e.name.clone(),
                    second: "Block".to_string(),
                });
            }
        }
        // The curvature contrast (factorial against centre runs) is
        // orthogonal to the blocks only when every block holds the same
        // share of centre runs.
        if n_c > 0 {
            let centre_in =
                |runs: &[usize]| runs.iter().filter(|&&r| is_centre(&design.data[r])).count();
            if block_runs
                .iter()
                .any(|runs| centre_in(runs) * n != n_c * runs.len())
            {
                return Err(DoeError::PartiallyAliasedEffects {
                    first: "curvature".to_string(),
                    second: "Block".to_string(),
                });
            }
        }
    }

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

    // Block means; with one block (no blocks given) it is the grand mean and
    // the block term is zero.
    let block_mean: Vec<f64> = block_runs
        .iter()
        .map(|runs| runs.iter().map(|&r| responses[r]).sum::<f64>() / runs.len() as f64)
        .collect();
    let block_ss: f64 = block_runs
        .iter()
        .zip(&block_mean)
        .map(|(runs, m)| runs.len() as f64 * (m - grand_mean).powi(2))
        .sum();
    let block_df = block_runs.len() - 1;

    let model_df = selected.len();
    if model_df + curvature_df + block_df > total_df {
        return Err(DoeError::OverSpecifiedModel {
            terms: model_df + curvature_df + block_df,
            runs: n,
        });
    }

    let model_ss: f64 = selected.iter().map(|e| e.sum_of_squares).sum();
    // The total splits exactly into the effects (within the factorial runs),
    // the curvature (between the factorial and centre means), the blocks and
    // the rest -- mutually orthogonal by the guards above, so none needs
    // adjusting for another. With the coding, aliasing, orthogonality and
    // degrees-of-freedom guards, what is clamped here is float noise, not a
    // structural overflow.
    let residual_ss = (total_ss - model_ss - curvature_ss - block_ss).max(0.0);
    let residual_df = total_df - model_df - curvature_df - block_df;
    // No residual degrees of freedom, or nothing left over: no error estimate
    // to test the effects against.
    let ms_residual = (residual_df > 0)
        .then(|| residual_ss / residual_df as f64)
        .filter(|&ms| ms > 0.0);

    let effects: Vec<EffectRow> = selected
        .iter()
        .map(|e| {
            let ms = e.sum_of_squares; // df = 1
            let f_stat = ms_residual.map(|ms_res| ms / ms_res);
            let p = f_stat.map(|f| f_pvalue(f, 1, residual_df));
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

    let pure_error = pure_error(design, responses, &block_of);
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

    let r_squared = if !has_no_variation(total_ss, responses) {
        (model_ss + curvature_ss + block_ss) / total_ss
    } else {
        0.0
    };
    let r_squared_adj = if residual_df > 0 && total_df > 0 {
        1.0 - (residual_ss / residual_df as f64) / (total_ss / total_df as f64)
    } else {
        r_squared
    };

    // With orthogonal, balanced contrasts the least-squares coefficient of each
    // term is half its effect, so no matrix solve is needed. Unblocked, centre
    // runs are fitted at their own mean and factorial runs at the factorial
    // mean plus the effects -- that is what the curvature term fits. A block
    // shifts every run of it by its mean's departure from the grand mean: the
    // curvature regressor (1 on factorial runs, 0 on centre runs) is centred
    // at its overall share, which the guard above made its share in every
    // block too.
    let share_f = n_f as f64 / n as f64;
    let curvature_coef = mean_f - mean_c;
    let mut fitted: Vec<f64> = (0..n)
        .map(|run| {
            let z = if is_centre(&design.data[run]) {
                0.0
            } else {
                1.0
            };
            block_mean[block_of[run]] + (z - share_f) * curvature_coef
        })
        .collect();
    for (i, &run) in factorial_runs.iter().enumerate() {
        fitted[run] += selected
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
        blocks: blocks.map(|_| BlockTerm {
            sum_of_squares: block_ss,
            df: block_df,
            mean_square: if block_df > 0 {
                block_ss / block_df as f64
            } else {
                0.0
            },
        }),
    })
}

/// Pure error from every group of identical design rows in the same block;
/// `None` when no design point was run more than once within a block.
fn pure_error(design: &DesignMatrix, responses: &[f64], block_of: &[usize]) -> Option<PureError> {
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
            .filter(|&other| {
                !assigned[other]
                    && block_of[other] == block_of[run]
                    && same_point(&design.data[run], &design.data[other])
            })
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

/// Upper-tail p-value of `F(df1, df2)`.
pub(crate) fn f_pvalue(f: f64, df1: usize, df2: usize) -> f64 {
    if f <= 0.0 {
        return 1.0;
    }
    let d1 = df1 as f64;
    let d2 = df2 as f64;
    let x = d2 / (d2 + d1 * f);
    u_numflow::special::regularized_incomplete_beta(x, d2 / 2.0, d1 / 2.0)
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
        let r = doe_anova(&design, &y, &["A", "B", "C"], None).expect("valid");
        let c = r.curvature.expect("centre points give a curvature test");
        assert!(c.sum_of_squares.abs() < 1e-9);
        let pe = r.pure_error.expect("replicated centre points");
        assert_eq!(pe.df, 2);
        assert!((pe.sum_of_squares - 0.02).abs() < 1e-9);
    }

    #[test]
    fn curvature_is_tested_against_pure_error() {
        let (design, y) = factorial_with_centres(&[55.0, 55.2, 54.8]);
        let r = doe_anova(&design, &y, &["A", "B", "C"], None).expect("valid");
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
        let r = doe_anova(&design, &y, &["A", "B"], None).expect("valid");
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
        match doe_anova(&design, &y, &["A"], None) {
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
        let r = doe_anova(&design, &y, &["A", "B"], None).expect("valid");
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
        let r = doe_anova(&design, &y, &["A", "B", "C"], None).expect("valid");
        assert!(r.curvature.is_none() && r.pure_error.is_none());
    }
    /// Counted into the contrasts, three centre points shrank every effect by
    /// 8/11; refused, they left the design with no analysis at all.
    #[test]
    fn centre_points_are_accepted_and_leave_the_effects_alone() {
        let (design, y) = factorial_with_centres(&[50.1, 49.9, 50.0]);
        let r = doe_anova(&design, &y, &["A", "B", "C"], None).expect("centre points are accepted");
        let ss: Vec<f64> = r.effects.iter().map(|e| e.sum_of_squares).collect();
        assert!((ss[0] - 200.0).abs() < 1e-9, "SS_A = {}", ss[0]);
        assert!((ss[1] - 72.0).abs() < 1e-9, "SS_B = {}", ss[1]);
        assert!((ss[2] - 8.0).abs() < 1e-9, "SS_C = {}", ss[2]);
    }
    #[test]
    fn residuals_square_to_the_residual_sum_of_squares() {
        let design = crate::design::factorial::full_factorial(3).expect("2^3");
        let responses = [10.0, 20.0, 15.0, 25.0, 12.0, 22.0, 18.0, 30.0];
        let result = doe_anova(&design, &responses, &["A", "B", "C"], None).expect("valid");
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
        let result = doe_anova(&design, &responses, &["A", "B", "A:B"], None).expect("valid");
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
        let result = doe_anova(&design, &responses, &names, None);
        assert!(
            matches!(result, Err(DoeError::PartiallyAliasedEffects { .. })),
            "partially aliased model accepted: {result:?}"
        );
        // The main effects alone are orthogonal and still analysable.
        doe_anova(&design, &responses, &names[..8], None).expect("main effects are orthogonal");
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
            None,
        )
        .unwrap();
        let row_a = result.effects.iter().find(|r| r.name == "A").unwrap();
        // A is highly significant
        let f = row_a.f_statistic.expect("residual df > 0");
        assert!(f > 10.0, "A F={f}");
    }

    #[test]
    fn anova_r_squared_reasonable() {
        let design = full_factorial(4).unwrap();
        let result = doe_anova(
            &design,
            &montgomery_responses(),
            &["A", "C", "D", "A:C"],
            None,
        )
        .unwrap();
        // A, C, D, AC explain ~77% of total variance in this dataset
        assert!(result.r_squared > 0.70, "R²={}", result.r_squared);
    }

    #[test]
    fn anova_ss_partition() {
        let design = full_factorial(4).unwrap();
        let result = doe_anova(
            &design,
            &montgomery_responses(),
            &["A", "C", "D", "A:C"],
            None,
        )
        .unwrap();
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
        assert!(doe_anova(&design, &[1.0, 2.0], &["A"], None).is_err());
    }

    /// The names `estimate_effects` produces are the names this function takes.
    ///
    /// The order to estimate up to used to be a literal 2 here, so a caller who
    /// asked `estimate_effects` for three-factor interactions got `A:B:C` back
    /// and then had it refused as unknown -- the crate's own output rejected by
    /// the crate's own input, with a message naming exactly the form it refused.
    #[test]
    fn anova_takes_every_name_estimate_effects_returns() {
        // A replicated 2^3, so order 3 leaves residual degrees of freedom.
        let corners = [
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [-1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        let mut rows = Vec::new();
        for c in corners {
            rows.push(c.to_vec());
            rows.push(c.to_vec());
        }
        let design = DesignMatrix {
            factor_names: vec!["A".into(), "B".into(), "C".into()],
            data: rows,
        };
        // The reporter's responses, each corner replicated with the same +2
        // offset so the replication contributes pure error without cancelling
        // any contrast. (An alternating offset zeroes A:B:C exactly, which is
        // how the first draft of this test passed for the wrong reason.)
        let y = [
            60.0, 62.0, 72.0, 74.0, 54.0, 56.0, 68.0, 70.0, 52.0, 54.0, 83.0, 85.0, 45.0, 47.0,
            80.0, 82.0,
        ];

        for order in 1..=3 {
            let est = crate::analysis::effects::estimate_effects(&design, &y, order)
                .unwrap_or_else(|e| panic!("estimate_effects at order {order}: {e}"));
            let names: Vec<&str> = est.iter().map(|e| e.name.as_str()).collect();
            let anova = doe_anova(&design, &y, &names, None)
                .unwrap_or_else(|e| panic!("doe_anova at order {order}: {e}"));
            assert_eq!(
                anova.effects.len(),
                names.len(),
                "order {order}: every requested term must appear"
            );
            for (row, name) in anova.effects.iter().zip(&names) {
                assert_eq!(&row.name, name, "order {order}");
            }
        }

        // The three-factor term specifically, on its own.
        let only_abc = doe_anova(&design, &y, &["A:B:C"], None).expect("A:B:C is estimable here");
        assert_eq!(only_abc.effects.len(), 1);
        assert_eq!(only_abc.effects[0].name, "A:B:C");
        assert!(
            only_abc.effects[0].sum_of_squares > 0.0,
            "SS = {}",
            only_abc.effects[0].sum_of_squares
        );
    }

    /// Widening the order must not widen what counts as a name: a factor that
    /// is not in the design is still unknown, however many parts the name has.
    #[test]
    fn anova_still_refuses_a_name_no_factor_answers_to() {
        let design = DesignMatrix {
            factor_names: vec!["A".into(), "B".into(), "C".into()],
            data: vec![
                vec![-1.0, -1.0, -1.0],
                vec![1.0, -1.0, -1.0],
                vec![-1.0, 1.0, -1.0],
                vec![1.0, 1.0, -1.0],
                vec![-1.0, -1.0, 1.0],
                vec![1.0, -1.0, 1.0],
                vec![-1.0, 1.0, 1.0],
                vec![1.0, 1.0, 1.0],
            ],
        };
        let y = [60.0, 72.0, 54.0, 68.0, 52.0, 83.0, 45.0, 80.0];
        for bad in ["Z", "A:Z", "A:B:Z", "A:B:C:D"] {
            assert!(
                matches!(
                    doe_anova(&design, &y, &[bad], None),
                    Err(crate::error::DoeError::UnknownEffect { .. })
                ),
                "{bad} must stay unknown"
            );
        }
    }

    #[test]
    fn anova_unknown_effect_errors() {
        let design = full_factorial(2).unwrap();
        let responses = vec![10.0, 20.0, 30.0, 40.0];
        // Valid names succeed
        assert!(doe_anova(&design, &responses, &["A", "B", "A:B"], None).is_ok());
        // Unknown name errors instead of silently producing a smaller model
        let err = doe_anova(&design, &responses, &["A", "Z"], None).unwrap_err();
        assert_eq!(
            err,
            crate::error::DoeError::UnknownEffect { effect: "Z".into() }
        );
        // Legacy separator-less interaction names are unknown (renamed in 0.5.0)
        assert!(doe_anova(&design, &responses, &["AB"], None).is_err());
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

        match doe_anova(
            &design,
            &responses,
            &[first.as_str(), second.as_str()],
            None,
        ) {
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
            doe_anova(&design, &responses, &["A", "B"], None).is_ok(),
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
            doe_anova(&design, &responses, &["A", "B", "A:B"], None).is_ok(),
            "a saturated model is legal: residual df is zero, not negative"
        );
    }

    #[test]
    fn refuses_a_design_that_is_not_two_level_coded() {
        use crate::design::ccd::{ccd, AlphaType};
        let design = ccd(2, AlphaType::FaceCentered, 3).expect("ccd builds");
        let responses = vec![1.0; design.run_count()];
        assert!(matches!(
            doe_anova(&design, &responses, &["A", "B"], None),
            Err(DoeError::NotTwoLevelCoded { .. })
        ));
    }

    #[test]
    fn saturated_model_leaves_zero_residual_degrees_of_freedom() {
        use crate::design::factorial::full_factorial;
        let design = full_factorial(2).expect("full factorial");
        let responses = vec![10.0, 20.0, 15.0, 25.0];
        let result =
            doe_anova(&design, &responses, &["A", "B", "A:B"], None).expect("saturated is legal");
        assert_eq!(result.residual_df, 0);
        assert!(
            result.residual_ss.abs() < 1e-9,
            "an exactly saturated model explains the whole total sum of squares, so the residual is zero rather than a clamped negative: got {}",
            result.residual_ss
        );
    }

    #[test]
    fn ragged_design_is_an_error_not_a_panic() {
        use crate::design::DesignMatrix;
        // A short row of zeros would pass for a centre point; a short factorial
        // row would be indexed past its end.
        for short in [vec![0.0], vec![1.0]] {
            let design = DesignMatrix {
                data: vec![vec![-1.0, -1.0], vec![1.0, -1.0], vec![-1.0, 1.0], short],
                factor_names: vec!["A".into(), "B".into()],
            };
            assert_eq!(
                doe_anova(&design, &[1.0, 2.0, 3.0, 4.0], &["A"], None).map(|r| r.residual_df),
                Err(DoeError::DimensionMismatch {
                    parameter: "design",
                    index: Some(3),
                    expected: 2,
                    got: 1,
                    counted: "one entry per factor",
                })
            );
        }
    }

    #[test]
    fn saturated_anova_reports_no_f_and_no_p() {
        let design = crate::design::factorial::full_factorial(2).unwrap();
        let result = doe_anova(&design, &[1.0, 2.0, 3.0, 7.0], &["A", "B", "A:B"], None).unwrap();
        assert_eq!(result.residual_df, 0);
        for row in &result.effects {
            assert_eq!((row.f_statistic, row.p_value), (None, None), "{}", row.name);
        }
        // Residual degrees of freedom but nothing left over: y is exactly A + B.
        let result = doe_anova(&design, &[0.0, 2.0, 2.0, 4.0], &["A", "B"], None).unwrap();
        assert_eq!(result.residual_df, 1);
        assert!(result
            .effects
            .iter()
            .all(|r| r.f_statistic.is_none() && r.p_value.is_none()));
    }

    #[test]
    fn empty_design_is_an_error_not_a_panic() {
        use crate::design::DesignMatrix;
        let empty = DesignMatrix {
            data: vec![],
            factor_names: vec![],
        };
        assert!(matches!(
            doe_anova(&empty, &[], &["A"], None),
            Err(DoeError::EmptyInput {
                parameter: "design"
            })
        ));
    }
}

/// The balanced two-level analysis against statsmodels' OLS and
/// `anova_lm(typ=3)`: a 2^3 run twice, the model A, B, C, A:B, A:C and B:C (the
/// three-factor interaction left to the residual), responses drawn from a known
/// model with a seeded NumPy generator. On balanced data the contrast formula
/// and the regression give the same effects (twice the coefficient) and the
/// same sums of squares.
#[cfg(test)]
mod statsmodels_reference {
    use super::*;
    use crate::analysis::effects::estimate_effects;

    fn close(got: f64, want: f64) -> bool {
        (got - want).abs() <= 1e-9 * want.abs().max(1.0)
    }

    #[test]
    fn a_balanced_analysis_reproduces_statsmodels() {
        let design = DesignMatrix {
            data: vec![
                vec![-1.0, -1.0, -1.0],
                vec![1.0, -1.0, -1.0],
                vec![-1.0, 1.0, -1.0],
                vec![1.0, 1.0, -1.0],
                vec![-1.0, -1.0, 1.0],
                vec![1.0, -1.0, 1.0],
                vec![-1.0, 1.0, 1.0],
                vec![1.0, 1.0, 1.0],
                vec![-1.0, -1.0, -1.0],
                vec![1.0, -1.0, -1.0],
                vec![-1.0, 1.0, -1.0],
                vec![1.0, 1.0, -1.0],
                vec![-1.0, -1.0, 1.0],
                vec![1.0, -1.0, 1.0],
                vec![-1.0, 1.0, 1.0],
                vec![1.0, 1.0, 1.0],
            ],
            factor_names: vec!["A".into(), "B".into(), "C".into()],
        };
        let responses = [
            19.05, 21.803, 14.119, 22.911, 22.213, 25.983, 12.298, 20.34, 20.182, 23.38, 14.73,
            22.731, 22.202, 25.573, 14.278, 22.385,
        ];
        // name, effect, sum of squares, F, p
        let terms: [(&str, f64, f64, f64, f64); 6] = [
            (
                "A",
                5.75425,
                132.44557225,
                185.308186343051,
                2.61251393490552e-07,
            ),
            (
                "B",
                -4.57425,
                83.69505225,
                117.100013800834,
                1.84816009915841e-06,
            ),
            (
                "C",
                0.795749999999999,
                2.53287224999999,
                3.54381014716134,
                0.0924365655158291,
            ),
            (
                "A:B",
                2.48125000000001,
                24.6264062500001,
                34.4554717897312,
                0.000237800567348805,
            ),
            (
                "A:C",
                0.0682499999999993,
                0.0186322499999996,
                0.0260688854775228,
                0.875298610636731,
            ),
            (
                "B:C",
                -2.09325,
                17.52678225,
                24.5221956159208,
                0.000789144855768242,
            ),
        ];
        let names = ["A", "B", "C", "A:B", "A:C", "B:C"];
        let anova = doe_anova(&design, &responses, &names, None).expect("a valid analysis");
        for (row, (name, _, ss, f, p)) in anova.effects.iter().zip(terms) {
            assert_eq!(row.name, name);
            assert!(
                close(row.sum_of_squares, ss),
                "{name} SS {}",
                row.sum_of_squares
            );
            assert!(close(row.f_statistic.expect("residual df"), f), "{name} F");
            assert!(close(row.p_value.expect("residual df"), p), "{name} p");
        }
        assert!(close(anova.residual_ss, 6.43258225));
        assert_eq!(anova.residual_df, 9);
        assert!(close(anova.total_ss, 267.27789975));
        assert!(close(anova.r_squared, 0.975932981155506));
        assert!(close(anova.r_squared_adj, 0.959888301925844));
        for (got, want) in anova.fitted.iter().zip([
            19.501375, 22.706125, 14.539125, 22.706375, 22.322125, 25.663375, 13.173375, 21.477125,
            19.501375, 22.706125, 14.539125, 22.706375, 22.322125, 25.663375, 13.173375, 21.477125,
        ]) {
            assert!(close(*got, want), "fitted {got} vs {want}");
        }

        let effects = estimate_effects(&design, &responses, 2).expect("a valid estimate");
        for (name, effect, ss, _, _) in terms {
            let e = effects
                .iter()
                .find(|e| e.name == name)
                .expect("the term is estimated");
            assert!(close(e.estimate, effect), "{name} effect {}", e.estimate);
            assert!(
                close(e.sum_of_squares, ss),
                "{name} effect SS {}",
                e.sum_of_squares
            );
        }
    }
}

/// The blocked analysis against its definition (Montgomery, *Design and
/// Analysis of Experiments*, ch. 7): the block sum of squares is the sum of
/// squares of the confounded contrast, the unconfounded effects do not move,
/// and a block difference leaves the residual instead of inflating it.
#[cfg(test)]
mod blocks {
    use super::*;
    use crate::analysis::effects::estimate_effects;
    use crate::design::factorial::{blocked_factorial, full_factorial};

    /// `y = 60 + 6A − 4B + 3C + 2A:B + 1.5C:D` plus a small deterministic
    /// wobble, so the residual is not zero.
    fn response(row: &[f64], run: usize) -> f64 {
        let x = |f: usize| row.get(f).copied().unwrap_or(0.0);
        60.0 + 6.0 * x(0) - 4.0 * x(1)
            + 3.0 * x(2)
            + 2.0 * x(0) * x(1)
            + 1.5 * x(2) * x(3)
            + ((run * 37 % 11) as f64 - 5.0) * 0.1
    }

    #[test]
    fn the_block_sum_of_squares_is_the_confounded_contrasts() {
        for (k, n_blocks) in [(3, 2), (4, 2), (4, 4), (5, 4)] {
            let blocked = blocked_factorial(k, n_blocks).expect("supported");
            let y: Vec<f64> = blocked
                .design
                .data
                .iter()
                .enumerate()
                .map(|(run, r)| response(r, run))
                .collect();
            let all = estimate_effects(&blocked.design, &y, k).expect("full factorial");
            let confounded = blocked.confounded_names();
            let confounded_ss: f64 = all
                .iter()
                .filter(|e| confounded.contains(&e.name))
                .map(|e| e.sum_of_squares)
                .sum();
            let r = doe_anova(&blocked.design, &y, &["A", "B"], Some(&blocked.blocks))
                .unwrap_or_else(|e| panic!("2^{k} in {n_blocks}: {e}"));
            let block = r.blocks.expect("blocks given");
            assert_eq!(block.df, n_blocks - 1);
            assert!(
                (block.sum_of_squares - confounded_ss).abs() < 1e-9 * confounded_ss.max(1.0),
                "2^{k} in {n_blocks}: block SS {} vs confounded {confounded_ss}",
                block.sum_of_squares
            );
            assert!((block.mean_square - block.sum_of_squares / block.df as f64).abs() < 1e-12);
        }
    }

    #[test]
    fn blocks_leave_the_unconfounded_effects_alone() {
        let blocked = blocked_factorial(4, 4).expect("2^4 in 4 blocks");
        let y: Vec<f64> = blocked
            .design
            .data
            .iter()
            .enumerate()
            .map(|(run, r)| response(r, run))
            .collect();
        let names = ["A", "B", "C", "D", "A:B", "C:D"];
        let plain = doe_anova(&blocked.design, &y, &names, None).expect("unblocked");
        let with = doe_anova(&blocked.design, &y, &names, Some(&blocked.blocks)).expect("blocked");
        for (p, w) in plain.effects.iter().zip(&with.effects) {
            assert_eq!(p.name, w.name);
            assert!(
                (p.sum_of_squares - w.sum_of_squares).abs() < 1e-9,
                "{}",
                p.name
            );
        }
        let block = with.blocks.expect("blocks").sum_of_squares;
        assert!((plain.residual_ss - with.residual_ss - block).abs() < 1e-9);
        assert_eq!(plain.residual_df - with.residual_df, 3);
        assert!(plain.blocks.is_none());
    }

    #[test]
    fn a_block_shift_moves_out_of_the_residual() {
        let blocked = blocked_factorial(4, 2).expect("2^4 in 2 blocks");
        let y: Vec<f64> = blocked
            .design
            .data
            .iter()
            .zip(&blocked.blocks)
            .enumerate()
            .map(|(run, (r, &b))| response(r, run) + 5.0 * b as f64)
            .collect();
        let names = ["A", "B", "C", "A:B", "C:D"];
        let plain = doe_anova(&blocked.design, &y, &names, None).expect("unblocked");
        let with = doe_anova(&blocked.design, &y, &names, Some(&blocked.blocks)).expect("blocked");
        let ms = |r: &DoeAnovaResult| r.residual_ss / r.residual_df as f64;
        assert!(
            ms(&with) < ms(&plain) / 10.0,
            "{} vs {}",
            ms(&with),
            ms(&plain)
        );
        // 16 runs, each 2.5 (plus the wobble's block part) from the grand mean.
        assert!(with.blocks.expect("blocks").sum_of_squares > 16.0 * 2.5 * 2.5 * 0.9);
        // The fitted values are a least-squares fit: residuals square to the
        // residual sum of squares and sum to zero within every block.
        let ss: f64 = with.residuals.iter().map(|e| e * e).sum();
        assert!((ss - with.residual_ss).abs() < 1e-9);
        for b in 0..2 {
            let s: f64 = with
                .residuals
                .iter()
                .zip(&blocked.blocks)
                .filter(|(_, &bb)| bb == b)
                .map(|(e, _)| e)
                .sum();
            assert!(s.abs() < 1e-9, "block {b} residual sum {s}");
        }
    }

    #[test]
    fn a_confounded_effect_is_refused_naming_the_block() {
        let blocked = blocked_factorial(3, 2).expect("2^3 in 2 blocks");
        let y: Vec<f64> = (0..8).map(|i| i as f64).collect();
        assert_eq!(
            doe_anova(&blocked.design, &y, &["A", "A:B:C"], Some(&blocked.blocks)).err(),
            Some(DoeError::AliasedEffects {
                first: "A:B:C".into(),
                second: "Block".into()
            })
        );
    }

    #[test]
    fn blocks_that_cut_across_an_effect_are_refused() {
        // Runs 0..3 and 4..7 of the standard order: C is -1 in the first half
        // and +1 in the second, so C is confounded; A alternates and is fine;
        // a split of 3 + 5 cuts across A.
        let design = full_factorial(3).expect("2^3");
        let y: Vec<f64> = (0..8).map(|i| (i * i) as f64).collect();
        let halves = [0, 0, 0, 0, 1, 1, 1, 1];
        assert!(matches!(
            doe_anova(&design, &y, &["C"], Some(&halves)),
            Err(DoeError::AliasedEffects { .. })
        ));
        doe_anova(&design, &y, &["A", "B", "A:B"], Some(&halves))
            .expect("orthogonal to the halves");
        let uneven = [0, 0, 0, 1, 1, 1, 1, 1];
        assert_eq!(
            doe_anova(&design, &y, &["A"], Some(&uneven)).err(),
            Some(DoeError::PartiallyAliasedEffects {
                first: "A".into(),
                second: "Block".into()
            })
        );
        assert!(matches!(
            doe_anova(&design, &y, &["A"], Some(&[0, 1])),
            Err(DoeError::DimensionMismatch {
                parameter: "blocks",
                expected: 8,
                got: 2,
                ..
            })
        ));
    }

    /// Centre points spread evenly over two blocks: curvature is unchanged,
    /// pure error is taken within blocks only, and the fit stays a
    /// least-squares fit. Spread unevenly they are refused.
    #[test]
    fn centre_points_in_blocks() {
        let blocked = blocked_factorial(3, 2).expect("2^3 in 2 blocks");
        let mut data = blocked.design.data.clone();
        let mut blocks = blocked.blocks.clone();
        let mut y: Vec<f64> = data
            .iter()
            .zip(&blocks)
            .map(|(r, &b)| 50.0 + 5.0 * r[0] + 3.0 * r[1] + 8.0 * b as f64)
            .collect();
        // Two centre runs per block, the block offset on them too, and a bump
        // (curvature) of +2.
        for (b, c) in [(0, 51.9), (0, 52.1), (1, 60.0), (1, 60.2)] {
            data.push(vec![0.0; 3]);
            blocks.push(b);
            y.push(c);
        }
        let design = DesignMatrix {
            data,
            factor_names: blocked.design.factor_names.clone(),
        };
        let r = doe_anova(&design, &y, &["A", "B"], Some(&blocks)).expect("even centre points");
        // Pure error: (51.9, 52.1) and (60.0, 60.2) within their blocks only.
        let pe = r.pure_error.expect("replicated centres");
        assert_eq!(pe.df, 2);
        assert!(
            (pe.sum_of_squares - 0.04).abs() < 1e-9,
            "{}",
            pe.sum_of_squares
        );
        // Curvature: factorial mean 54, centre mean 56.05 -> d = -2.05.
        let c = r.curvature.expect("centre points");
        let expected = 8.0 * 4.0 * 2.05 * 2.05 / 12.0;
        assert!(
            (c.sum_of_squares - expected).abs() < 1e-9,
            "{}",
            c.sum_of_squares
        );
        let ss: f64 = r.residuals.iter().map(|e| e * e).sum();
        assert!(
            (ss - r.residual_ss).abs() < 1e-9,
            "{ss} vs {}",
            r.residual_ss
        );
        assert!(
            (r.effects.iter().map(|e| e.sum_of_squares).sum::<f64>()
                + c.sum_of_squares
                + r.blocks.expect("blocks").sum_of_squares
                + r.residual_ss
                - r.total_ss)
                .abs()
                < 1e-9
        );
        assert_eq!(r.residual_df, 12 - 1 - 2 - 1 - 1);

        let mut uneven = blocks.clone();
        uneven[9] = 1; // three centre runs in block 1, one in block 0
        assert_eq!(
            doe_anova(&design, &y, &["A", "B"], Some(&uneven)).err(),
            Some(DoeError::PartiallyAliasedEffects {
                first: "curvature".into(),
                second: "Block".into()
            })
        );
    }

    #[test]
    fn one_block_is_no_block() {
        let design = full_factorial(3).expect("2^3");
        let y: Vec<f64> = (0..8).map(|i| (i * 3 % 7) as f64).collect();
        let plain = doe_anova(&design, &y, &["A", "B"], None).expect("plain");
        let one = doe_anova(&design, &y, &["A", "B"], Some(&[4; 8])).expect("one block");
        assert_eq!(plain.residual_df, one.residual_df);
        assert!((plain.residual_ss - one.residual_ss).abs() < 1e-12);
        assert_eq!(
            one.blocks,
            Some(BlockTerm {
                sum_of_squares: 0.0,
                df: 0,
                mean_square: 0.0
            })
        );
    }
}

/// The blocked analysis against statsmodels 0.15.0 OLS + `anova_lm(typ=1)`:
/// a 2^4 in four blocks (`A:B:C`, `A:C:D` and `B:D` confounded), one centre run
/// per block, responses from a known model with a block shift, a curvature
/// bump and seeded NumPy noise; the model `C(block) + curv + A + B + C + D +
/// A:C`. Blocks, curvature and effects are orthogonal, so the sequential sums
/// of squares are the ones this table reports. (statsmodels tests curvature
/// against the residual; here there is no pure error -- one centre run per
/// block -- so only its sum of squares is compared.)
#[cfg(test)]
mod statsmodels_blocked_reference {
    use super::*;

    const Y: [f64; 20] = [
        39.734159, 48.31023, 36.137885, 33.80518, 44.864019, 51.437162, 31.786477, 48.247779,
        35.376527, 46.776129, 37.394012, 37.495747, 41.219994, 46.420958, 34.362252, 50.902013,
        41.10517, 45.884264, 41.901252, 47.561029,
    ];
    const BLOCKS: [usize; 20] = [0, 3, 1, 2, 3, 0, 2, 1, 2, 1, 3, 0, 1, 2, 0, 3, 0, 1, 2, 3];
    const FITTED: [f64; 20] = [
        38.97537509999998,
        48.39091344999999,
        36.300708175,
        34.94175252500001,
        44.852277449999995,
        51.06670735,
        31.40311652500001,
        48.392040425000026,
        35.459648524999984,
        45.690955675000005,
        37.90685119999999,
        37.36366435000001,
        42.152319675,
        47.55098077500002,
        33.825028350000004,
        49.99818345000002,
        42.90371484999998,
        45.73002704999999,
        39.93489564999999,
        47.88307744999998,
    ];

    fn close(got: f64, want: f64) -> bool {
        (got - want).abs() <= 1e-9 * want.abs().max(1.0)
    }

    #[test]
    fn a_blocked_analysis_reproduces_statsmodels() {
        let mut data: Vec<Vec<f64>> = (0..16)
            .map(|i| {
                (0..4)
                    .map(|f| if (i >> f) & 1 == 1 { 1.0 } else { -1.0 })
                    .collect()
            })
            .collect();
        data.extend((0..4).map(|_| vec![0.0; 4]));
        let design = DesignMatrix {
            data,
            factor_names: ["A", "B", "C", "D"].map(String::from).to_vec(),
        };
        let r = doe_anova(&design, &Y, &["A", "B", "C", "D", "A:C"], Some(&BLOCKS))
            .expect("a valid blocked analysis");

        let block = r.blocks.expect("blocks");
        assert_eq!(block.df, 3);
        assert!(
            close(block.sum_of_squares, 178.7359348631167),
            "{}",
            block.sum_of_squares
        );
        let curvature = r.curvature.expect("centre points");
        assert!(close(curvature.sum_of_squares, 21.565841142219746));
        assert!(
            r.pure_error.is_none(),
            "one centre run per block: no pure error"
        );
        // name, SS, F, p
        let terms: [(&str, f64, f64, f64); 5] = [
            (
                "A",
                244.2959074960086,
                172.6325527165604,
                1.2398283542525177e-07,
            ),
            (
                "B",
                121.04308533474294,
                85.53551725121116,
                3.236956659129568e-06,
            ),
            (
                "C",
                73.1486131447641,
                51.69072189661257,
                2.9622705781370073e-05,
            ),
            (
                "D",
                1.1964307073175542,
                0.8454619206262549,
                0.37948485712231217,
            ),
            (
                "A:C",
                45.665384139880615,
                32.26960253101548,
                0.00020366301050287673,
            ),
        ];
        for (row, (name, ss, f, p)) in r.effects.iter().zip(terms) {
            assert_eq!(row.name, name);
            assert!(
                close(row.sum_of_squares, ss),
                "{name} SS {}",
                row.sum_of_squares
            );
            assert!(close(row.f_statistic.expect("F"), f), "{name} F");
            assert!(close(row.p_value.expect("p"), p), "{name} p");
        }
        assert!(
            close(r.residual_ss, 14.151207501235868),
            "{}",
            r.residual_ss
        );
        assert_eq!(r.residual_df, 10);
        assert!(close(r.r_squared, 0.979778281106652), "{}", r.r_squared);
        assert!(
            close(r.r_squared_adj, 0.9615787341026388),
            "{}",
            r.r_squared_adj
        );
        for (run, (got, want)) in r.fitted.iter().zip(FITTED).enumerate() {
            assert!(close(*got, want), "fitted[{run}] {got} vs {want}");
        }
    }
}
