//! Derringer-Suich desirability functions for multi-response optimization.
//!
//! Transforms each response into a dimensionless desirability value d ∈ [0, 1],
//! then combines them into an overall desirability D via geometric mean.
//!
//! # Individual Desirability
//!
//! Three goal types are supported: maximize, minimize, and hit-a-target.
//! The weight parameter s (or s1/s2) controls curve shape:
//! - s = 1: linear ramp
//! - s < 1: concave (tolerance for deviation)
//! - s > 1: convex (strict requirement)
//!
//! Each ramp divides by its width -- `target - lower` rising, `upper - target`
//! falling -- so a specification is only meaningful when the limits a goal
//! uses are finite and strictly ordered and the exponents are positive. A
//! [`ResponseSpec`] can only be built in that state: the constructors refuse
//! anything else, which is why [`ResponseSpec::desirability`] cannot fail.
//!
//! # Overall Desirability
//!
//! D = (∏ dᵢ^rᵢ)^(1/Σrᵢ)
//!
//! where rᵢ is each response's relative importance weight (default 1.0, giving the
//! equal-weighted geometric mean D = (∏ dᵢ)^(1/m)). If any dᵢ = 0 (with rᵢ > 0),
//! then D = 0.
//!
//! Reference: Derringer, G. & Suich, R. (1980). "Simultaneous Optimization
//! of Several Response Variables". *Journal of Quality Technology* 12(4),
//! pp. 214–219.

use std::ops::Bound::{self, Excluded, Unbounded};

use crate::error::DoeError;

/// What a curve-shape exponent takes: greater than 0.
const POSITIVE: (Bound<f64>, Bound<f64>) = (Excluded(0.0), Unbounded);

/// Optimization goal for a single response.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "wasm", derive(serde::Serialize))]
pub enum GoalType {
    /// Maximize the response (larger is better). Uses `lower` and `target`.
    Maximize,
    /// Minimize the response (smaller is better). Uses `target` and `upper`.
    Minimize,
    /// Hit a target value. Uses `lower`, `target` and `upper`.
    Target,
}

/// Specification for a single response's desirability function.
///
/// Built by [`ResponseSpec::new`] or the goal-specific constructors, which
/// refuse limits and exponents the ramp has no meaning for.
#[derive(Debug, Clone, PartialEq)]
pub struct ResponseSpec {
    goal: GoalType,
    lower: f64,
    target: f64,
    upper: f64,
    s1: f64,
    s2: f64,
    importance: f64,
}

impl ResponseSpec {
    /// Create a specification for `goal`.
    ///
    /// - `lower` (L): below it d = 0 for `Maximize` and `Target`.
    /// - `target` (T): d = 1 there (and beyond it, for the one-sided goals).
    /// - `upper` (U): above it d = 0 for `Minimize` and `Target`.
    /// - `s1`: curve-shape exponent of the rising ramp (`Maximize`, and
    ///   `Target` left of T) or of the falling ramp (`Minimize`).
    /// - `s2`: curve-shape exponent of `Target`'s falling ramp.
    ///
    /// A limit or exponent the goal does not use is not examined. Importance is
    /// 1.0; see [`ResponseSpec::with_importance`].
    ///
    /// # Errors
    /// [`DoeError::InvalidDesirabilityLimits`] unless the limits the goal uses
    /// are finite and `lower < target` (`Maximize`), `target < upper`
    /// (`Minimize`) or `lower < target < upper` (`Target`) -- a ramp of no width
    /// or running backwards makes d a step or a constant.
    /// [`DoeError::ParameterOutOfRange`] (`exclusive_min` 0) unless each
    /// exponent the goal uses is positive -- `s = 0` makes every point of the
    /// ramp d = 1, and a negative `s` sends d above 1 -- and
    /// [`DoeError::ValueNotFinite`] if one is NaN or infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use u_doe::optimization::desirability::{GoalType, ResponseSpec};
    /// let spec = ResponseSpec::new(GoalType::Target, 60.0, 80.0, 100.0, 1.0, 1.0).unwrap();
    /// assert!((spec.desirability(90.0) - 0.5).abs() < 1e-9);
    /// assert!(ResponseSpec::new(GoalType::Target, 10.0, 5.0, 0.0, 1.0, 1.0).is_err());
    /// ```
    pub fn new(
        goal: GoalType,
        lower: f64,
        target: f64,
        upper: f64,
        s1: f64,
        s2: f64,
    ) -> Result<Self, DoeError> {
        let rising = lower.is_finite() && target.is_finite() && lower < target;
        let falling = target.is_finite() && upper.is_finite() && target < upper;
        let ramps_ok = match goal {
            GoalType::Maximize => rising,
            GoalType::Minimize => falling,
            GoalType::Target => rising && falling,
        };
        if !ramps_ok {
            return Err(DoeError::InvalidDesirabilityLimits {
                index: None,
                goal,
                lower,
                target,
                upper,
            });
        }
        DoeError::check_real("s1", s1, POSITIVE)?;
        if goal == GoalType::Target {
            DoeError::check_real("s2", s2, POSITIVE)?;
        }
        Ok(Self {
            goal,
            lower,
            target,
            upper,
            s1,
            s2,
            importance: 1.0,
        })
    }

    /// Create a "larger-is-better" specification.
    ///
    /// d = 0 for y < lower; d = 1 for y ≥ target; a ramp with exponent `s` in
    /// between.
    ///
    /// # Errors
    /// As [`ResponseSpec::new`].
    ///
    /// # Examples
    ///
    /// ```
    /// use u_doe::optimization::desirability::ResponseSpec;
    /// let spec = ResponseSpec::maximize(50.0, 100.0, 1.0).unwrap();
    /// assert!((spec.desirability(75.0) - 0.5).abs() < 1e-9);
    /// ```
    pub fn maximize(lower: f64, target: f64, s: f64) -> Result<Self, DoeError> {
        Self::new(GoalType::Maximize, lower, target, target, s, s)
    }

    /// Create a "smaller-is-better" specification.
    ///
    /// d = 1 for y ≤ target; d = 0 for y > upper; a ramp with exponent `s` in
    /// between.
    ///
    /// # Errors
    /// As [`ResponseSpec::new`].
    ///
    /// # Examples
    ///
    /// ```
    /// use u_doe::optimization::desirability::ResponseSpec;
    /// let spec = ResponseSpec::minimize(0.0, 50.0, 1.0).unwrap();
    /// assert!((spec.desirability(25.0) - 0.5).abs() < 1e-9);
    /// ```
    pub fn minimize(target: f64, upper: f64, s: f64) -> Result<Self, DoeError> {
        Self::new(GoalType::Minimize, target, target, upper, s, s)
    }

    /// Create a "hit-a-target" specification.
    ///
    /// d = 1 at y = target; d = 0 outside [lower, upper].
    /// s1 controls the left ramp shape; s2 controls the right ramp.
    ///
    /// # Errors
    /// As [`ResponseSpec::new`].
    ///
    /// # Examples
    ///
    /// ```
    /// use u_doe::optimization::desirability::ResponseSpec;
    /// let spec = ResponseSpec::target(60.0, 80.0, 100.0, 1.0, 1.0).unwrap();
    /// assert!((spec.desirability(80.0) - 1.0).abs() < 1e-10);
    /// ```
    pub fn target(lower: f64, target: f64, upper: f64, s1: f64, s2: f64) -> Result<Self, DoeError> {
        Self::new(GoalType::Target, lower, target, upper, s1, s2)
    }

    /// Set this response's relative importance weight rᵢ (Derringer & Suich 1980),
    /// used by [`overall_desirability`]. Defaults to 1.0; 0 leaves the response
    /// out of the overall value.
    ///
    /// # Errors
    /// [`DoeError::ValueNotFinite`] if `importance` is NaN or infinite,
    /// [`DoeError::ParameterOutOfRange`] (`min` 0) if it is negative.
    ///
    /// # Examples
    ///
    /// ```
    /// use u_doe::optimization::desirability::ResponseSpec;
    /// let spec = ResponseSpec::maximize(0.0, 100.0, 1.0)
    ///     .and_then(|s| s.with_importance(3.0))
    ///     .unwrap();
    /// assert!((spec.importance() - 3.0).abs() < 1e-12);
    /// ```
    pub fn with_importance(mut self, importance: f64) -> Result<Self, DoeError> {
        DoeError::check_real("importance", importance, 0.0..)?;
        self.importance = importance;
        Ok(self)
    }

    /// Optimization goal.
    pub fn goal(&self) -> GoalType {
        self.goal
    }

    /// Lower acceptability limit (L); `None` for `Minimize`, which has none.
    pub fn lower(&self) -> Option<f64> {
        (self.goal != GoalType::Minimize).then_some(self.lower)
    }

    /// Target value (T).
    pub fn target_value(&self) -> f64 {
        self.target
    }

    /// Upper acceptability limit (U); `None` for `Maximize`, which has none.
    pub fn upper(&self) -> Option<f64> {
        (self.goal != GoalType::Maximize).then_some(self.upper)
    }

    /// Exponent of the rising ramp (the falling one for `Minimize`).
    pub fn s1(&self) -> f64 {
        self.s1
    }

    /// Exponent of the falling ramp of `Target`; `None` for the other goals,
    /// which have one ramp.
    pub fn s2(&self) -> Option<f64> {
        (self.goal == GoalType::Target).then_some(self.s2)
    }

    /// Relative importance weight rᵢ.
    pub fn importance(&self) -> f64 {
        self.importance
    }

    /// Compute the individual desirability for observed value `y`.
    pub fn desirability(&self, y: f64) -> f64 {
        match self.goal {
            GoalType::Maximize => {
                if y < self.lower {
                    0.0
                } else if y >= self.target {
                    1.0
                } else {
                    ((y - self.lower) / (self.target - self.lower)).powf(self.s1)
                }
            }
            GoalType::Minimize => {
                if y <= self.target {
                    1.0
                } else if y > self.upper {
                    0.0
                } else {
                    ((self.upper - y) / (self.upper - self.target)).powf(self.s1)
                }
            }
            GoalType::Target => {
                if y < self.lower || y > self.upper {
                    0.0
                } else if y <= self.target {
                    ((y - self.lower) / (self.target - self.lower)).powf(self.s1)
                } else {
                    ((self.upper - y) / (self.upper - self.target)).powf(self.s2)
                }
            }
        }
    }
}

/// Compute the overall desirability D from multiple response specifications.
///
/// D = (∏ dᵢ^rᵢ)^(1/Σrᵢ), the Derringer-Suich importance-weighted geometric mean,
/// where rᵢ = `spec.importance()` (default 1.0). With all weights equal this
/// reduces to the plain geometric mean (∏ dᵢ)^(1/m). D is 0 if any response with
/// a positive weight has d = 0; a response with importance 0 is left out.
///
/// # Errors
/// [`DoeError::EmptyInput`] (`specs`) if there are no specifications;
/// [`DoeError::DimensionMismatch`] unless there is one response per
/// specification; [`DoeError::NoWeightedResponse`] if every importance is 0.
///
/// # Examples
///
/// ```
/// use u_doe::optimization::desirability::{ResponseSpec, overall_desirability};
/// let specs = vec![
///     ResponseSpec::maximize(0.0, 100.0, 1.0).unwrap(),
///     ResponseSpec::maximize(0.0, 100.0, 1.0).unwrap(),
/// ];
/// let d = overall_desirability(&specs, &[100.0, 100.0]).unwrap();
/// assert!((d - 1.0).abs() < 1e-10);
/// ```
pub fn overall_desirability(specs: &[ResponseSpec], responses: &[f64]) -> Result<f64, DoeError> {
    if specs.is_empty() {
        return Err(DoeError::EmptyInput { parameter: "specs" });
    }
    if specs.len() != responses.len() {
        return Err(DoeError::DimensionMismatch {
            parameter: "responses",
            index: None,
            expected: specs.len(),
            got: responses.len(),
            counted: "one response per specification",
        });
    }
    let sum_r: f64 = specs.iter().map(|spec| spec.importance).sum();
    if sum_r <= 0.0 {
        return Err(DoeError::NoWeightedResponse);
    }
    // Weighted geometric mean via logs: D = exp( (Σ rᵢ·ln dᵢ) / Σrᵢ ). A dᵢ = 0
    // with positive weight vetoes the whole design (D = 0), matching the
    // unweighted convention; a zero weight drops the response entirely.
    let mut weighted_log_sum = 0.0;
    for (spec, &y) in specs.iter().zip(responses.iter()) {
        let r = spec.importance;
        if r == 0.0 {
            continue;
        }
        let d = spec.desirability(y);
        if d <= 0.0 {
            return Ok(0.0);
        }
        weighted_log_sum += r * d.ln();
    }
    Ok((weighted_log_sum / sum_r).exp())
}

/// The outcome of a search for the most desirable candidate.
///
/// Two cases, not one with a sentinel, because an overall desirability of 0 is
/// **not an optimum**: D is a weighted geometric mean, so one response that
/// misses its limits everywhere makes every candidate score exactly 0. A search
/// that kept the maximum would then return whichever candidate it scored first
/// -- a setting that satisfies no specification, presented as the recommended
/// one. The type makes that unrepresentable.
#[derive(Debug, Clone, PartialEq)]
pub enum DesirabilityOptimum {
    /// A candidate whose overall desirability is above 0.
    Best {
        /// Position of the candidate in the input.
        index: usize,
        /// Its overall desirability D.
        overall: f64,
        /// Each response's individual desirability dᵢ, in specification order.
        individual: Vec<f64>,
    },
    /// No candidate scores above 0.
    Infeasible {
        /// Specifications whose desirability is 0 at *every* candidate, in
        /// specification order. Empty when each response is met somewhere but
        /// never all at once -- the distinction that tells a user whether to
        /// move one response's limits or to relax the trade-off.
        ///
        /// A specification with importance 0 is left out of D, so it is never
        /// listed here however it scores.
        unreachable: Vec<usize>,
    },
}

/// Find the candidate with the highest overall desirability.
///
/// `candidates` is one response vector per candidate, each in specification
/// order. Ties go to the first candidate that is strictly best.
///
/// Scoring one vector is [`overall_desirability`]; this is the search over
/// settings that Derringer & Suich's method ends in, and the reason it lives
/// here rather than in each caller is [`DesirabilityOptimum`]'s doc comment.
///
/// # Errors
/// [`DoeError::EmptyInput`] (`parameter` naming which) if `specs` or
/// `candidates` is empty;
/// [`DoeError::DimensionMismatch`] if a candidate is not one response per
/// specification; [`DoeError::ValueNotFinite`] (`parameter` `"candidates"`,
/// `index` the response's position in the candidates laid end to end) if a
/// response is not finite --
/// a NaN is not a score, and letting it fall through as d = 0 would report a
/// reachable response as unreachable; [`DoeError::NoWeightedResponse`] if every
/// importance is 0.
///
/// # Examples
/// ```
/// use u_doe::optimization::desirability::{
///     optimize_desirability, DesirabilityOptimum, ResponseSpec,
/// };
///
/// // Yield to maximise over [0, 100]; the third candidate is the best of three.
/// let specs = vec![ResponseSpec::maximize(0.0, 100.0, 1.0).unwrap()];
/// let found = optimize_desirability(&specs, &[vec![10.0], vec![50.0], vec![90.0]]).unwrap();
/// assert!(matches!(found, DesirabilityOptimum::Best { index: 2, .. }));
///
/// // The same candidates against a target no one reaches: not an optimum.
/// let specs = vec![ResponseSpec::maximize(200.0, 300.0, 1.0).unwrap()];
/// let found = optimize_desirability(&specs, &[vec![10.0], vec![50.0], vec![90.0]]).unwrap();
/// assert_eq!(found, DesirabilityOptimum::Infeasible { unreachable: vec![0] });
/// ```
pub fn optimize_desirability(
    specs: &[ResponseSpec],
    candidates: &[Vec<f64>],
) -> Result<DesirabilityOptimum, DoeError> {
    if specs.is_empty() {
        return Err(DoeError::EmptyInput { parameter: "specs" });
    }
    if candidates.is_empty() {
        return Err(DoeError::EmptyInput {
            parameter: "candidates",
        });
    }
    for (row, candidate) in candidates.iter().enumerate() {
        if candidate.len() != specs.len() {
            return Err(DoeError::DimensionMismatch {
                parameter: "candidates",
                index: Some(row),
                expected: specs.len(),
                got: candidate.len(),
                counted: "one response per specification",
            });
        }
        // Reported as a position in the candidates laid end to end, the
        // flat list the WebAssembly binding takes.
        if let Some(column) = candidate.iter().position(|y| !y.is_finite()) {
            return Err(DoeError::ValueNotFinite {
                parameter: "candidates",
                index: Some(row * specs.len() + column),
            });
        }
    }

    // `overall_desirability` refuses an all-zero weighting; ask it once so the
    // search does not report an infeasible region for a weighting problem.
    let mut best: Option<(usize, f64)> = None;
    for (index, candidate) in candidates.iter().enumerate() {
        let overall = overall_desirability(specs, candidate)?;
        if best.is_none_or(|(_, b)| overall > b) {
            best = Some((index, overall));
        }
    }
    let (index, overall) = best.expect("candidates is not empty");

    if overall > 0.0 {
        return Ok(DesirabilityOptimum::Best {
            index,
            overall,
            individual: specs
                .iter()
                .zip(candidates[index].iter())
                .map(|(spec, &y)| spec.desirability(y))
                .collect(),
        });
    }

    // Nothing scores. Name the responses that are the reason -- those that are
    // 0 everywhere -- and say nothing when there are none, which means each is
    // satisfiable alone and the candidates never satisfy them together.
    let unreachable = specs
        .iter()
        .enumerate()
        .filter(|(_, spec)| spec.importance() > 0.0)
        .filter(|(i, spec)| {
            candidates
                .iter()
                .all(|candidate| spec.desirability(candidate[*i]) <= 0.0)
        })
        .map(|(i, _)| i)
        .collect();
    Ok(DesirabilityOptimum::Infeasible { unreachable })
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── the search, and what it refuses to call an optimum ──────────────────

    /// The reporter's first row: Yield can never reach its lower limit, so
    /// every candidate scores 0. A loop keeping the maximum returns candidate
    /// 0 and calls it optimal; this says which response is the reason.
    #[test]
    fn a_response_that_is_never_met_is_named_rather_than_averaged_away() {
        let specs = vec![ResponseSpec::maximize(200.0, 300.0, 1.0).expect("ordered")];
        let candidates = vec![vec![10.0], vec![50.0], vec![90.0]];
        for c in &candidates {
            assert_eq!(overall_desirability(&specs, c).expect("valid"), 0.0);
        }
        assert_eq!(
            optimize_desirability(&specs, &candidates).expect("valid"),
            DesirabilityOptimum::Infeasible {
                unreachable: vec![0]
            }
        );
    }

    /// The reporter's second row: each response is met somewhere, never both at
    /// once. `unreachable` is empty, which is the distinction that tells a user
    /// to relax the trade-off rather than move a limit.
    #[test]
    fn responses_met_apart_but_never_together_name_no_single_reason() {
        let specs = vec![
            ResponseSpec::maximize(60.0, 95.0, 1.0).expect("ordered"),
            ResponseSpec::minimize(4.0, 8.0, 1.0).expect("ordered"),
        ];
        // Yield meets its limits only at run 3, Cost only at run 1.
        let candidates = vec![vec![10.0, 5.0], vec![50.0, 10.0], vec![90.0, 15.0]];
        assert_eq!(
            optimize_desirability(&specs, &candidates).expect("valid"),
            DesirabilityOptimum::Infeasible {
                unreachable: vec![]
            }
        );
    }

    #[test]
    fn the_best_candidate_is_the_first_strictly_best_one() {
        let specs = vec![ResponseSpec::maximize(0.0, 100.0, 1.0).expect("ordered")];
        let candidates = vec![vec![10.0], vec![90.0], vec![90.0], vec![50.0]];
        let found = optimize_desirability(&specs, &candidates).expect("valid");
        match found {
            DesirabilityOptimum::Best {
                index,
                overall,
                individual,
            } => {
                assert_eq!(index, 1, "a tie goes to the first");
                assert!((overall - 0.9).abs() < 1e-12, "D = {overall}");
                assert_eq!(individual.len(), 1);
                assert!((individual[0] - 0.9).abs() < 1e-12);
            }
            other => panic!("{other:?}"),
        }
    }

    /// A response left out of D cannot be the reason D is 0.
    #[test]
    fn a_response_with_no_importance_is_never_the_reason() {
        let specs = vec![
            ResponseSpec::maximize(0.0, 100.0, 1.0)
                .expect("ordered")
                .with_importance(1.0)
                .expect("positive"),
            ResponseSpec::maximize(200.0, 300.0, 1.0)
                .expect("ordered")
                .with_importance(0.0)
                .expect("zero is allowed"),
        ];
        // The second response is unreachable, but weightless: D comes from the
        // first alone and is positive.
        let found = optimize_desirability(&specs, &[vec![50.0, 10.0]]).expect("valid");
        assert!(
            matches!(found, DesirabilityOptimum::Best { index: 0, .. }),
            "{found:?}"
        );

        // And when nothing scores, it is still not named.
        let specs = vec![
            ResponseSpec::maximize(200.0, 300.0, 1.0).expect("ordered"),
            ResponseSpec::maximize(200.0, 300.0, 1.0)
                .expect("ordered")
                .with_importance(0.0)
                .expect("zero is allowed"),
        ];
        assert_eq!(
            optimize_desirability(&specs, &[vec![10.0, 10.0]]).expect("valid"),
            DesirabilityOptimum::Infeasible {
                unreachable: vec![0]
            }
        );
    }

    #[test]
    fn the_search_refuses_what_it_cannot_score() {
        let specs = vec![ResponseSpec::maximize(0.0, 100.0, 1.0).expect("ordered")];
        assert!(matches!(
            optimize_desirability(&specs, &[]),
            Err(DoeError::EmptyInput {
                parameter: "candidates"
            })
        ));
        assert!(matches!(
            optimize_desirability(&[], &[vec![1.0]]),
            Err(DoeError::EmptyInput { parameter: "specs" })
        ));
        assert!(matches!(
            optimize_desirability(&specs, &[vec![1.0, 2.0]]),
            Err(DoeError::DimensionMismatch {
                parameter: "candidates",
                index: Some(0),
                expected: 1,
                got: 2,
                ..
            })
        ));
        for bad in [f64::NAN, f64::INFINITY] {
            assert_eq!(
                optimize_desirability(&specs, &[vec![50.0], vec![bad]]),
                Err(DoeError::ValueNotFinite {
                    parameter: "candidates",
                    index: Some(1),
                }),
                "{bad} must not be scored"
            );
        }
    }

    fn maximize(lower: f64, target: f64, s: f64) -> ResponseSpec {
        ResponseSpec::maximize(lower, target, s).expect("valid maximize spec")
    }

    fn minimize(target: f64, upper: f64, s: f64) -> ResponseSpec {
        ResponseSpec::minimize(target, upper, s).expect("valid minimize spec")
    }

    fn target(lower: f64, t: f64, upper: f64, s1: f64, s2: f64) -> ResponseSpec {
        ResponseSpec::target(lower, t, upper, s1, s2).expect("valid target spec")
    }

    #[test]
    fn maximize_at_target_is_one() {
        let spec = maximize(50.0, 100.0, 1.0);
        assert!((spec.desirability(100.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn maximize_below_lower_is_zero() {
        let spec = maximize(50.0, 100.0, 1.0);
        assert!(spec.desirability(49.9).abs() < 1e-10);
    }

    #[test]
    fn maximize_midpoint_linear() {
        // s=1, linear: d(75) = (75-50)/(100-50) = 0.5
        let spec = maximize(50.0, 100.0, 1.0);
        assert!((spec.desirability(75.0) - 0.5).abs() < 1e-9);
    }

    #[test]
    fn minimize_at_target_is_one() {
        let spec = minimize(0.0, 50.0, 1.0);
        assert!((spec.desirability(0.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn minimize_above_upper_is_zero() {
        let spec = minimize(0.0, 50.0, 1.0);
        assert!(spec.desirability(50.1).abs() < 1e-10);
    }

    #[test]
    fn target_at_target_is_one() {
        let spec = target(60.0, 80.0, 100.0, 1.0, 1.0);
        assert!((spec.desirability(80.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn target_below_lower_is_zero() {
        let spec = target(60.0, 80.0, 100.0, 1.0, 1.0);
        assert!(spec.desirability(59.9).abs() < 1e-10);
    }

    #[test]
    fn target_above_upper_is_zero() {
        let spec = target(60.0, 80.0, 100.0, 1.0, 1.0);
        assert!(spec.desirability(100.1).abs() < 1e-10);
    }

    #[test]
    fn target_left_midpoint() {
        // s1=1, linear: d(70) = (70-60)/(80-60) = 0.5
        let spec = target(60.0, 80.0, 100.0, 1.0, 1.0);
        assert!((spec.desirability(70.0) - 0.5).abs() < 1e-9);
    }

    #[test]
    fn target_right_midpoint() {
        // s2=1: d(90) = (100-90)/(100-80) = 0.5
        let spec = target(60.0, 80.0, 100.0, 1.0, 1.0);
        assert!((spec.desirability(90.0) - 0.5).abs() < 1e-9);
    }

    #[test]
    fn overall_geometric_mean() {
        let specs = vec![maximize(0.0, 100.0, 1.0), maximize(0.0, 100.0, 1.0)];
        // d1=0.5 (at 50), d2=1.0 (at 100) → D = (0.5 * 1.0)^(1/2) = sqrt(0.5) ≈ 0.707
        let d = overall_desirability(&specs, &[50.0, 100.0]).unwrap();
        assert!((d - 0.5_f64.sqrt()).abs() < 1e-9);
    }

    #[test]
    fn overall_zero_if_any_zero() {
        let specs = vec![maximize(50.0, 100.0, 1.0), maximize(50.0, 100.0, 1.0)];
        // Second response below L → d=0 → D=0
        let d = overall_desirability(&specs, &[80.0, 10.0]).unwrap();
        assert!(d.abs() < 1e-10);
    }

    #[test]
    fn desirability_weight_s2_convex() {
        // s=2 (convex): d(75) = ((75-50)/50)^2 = 0.25
        let spec = maximize(50.0, 100.0, 2.0);
        assert!((spec.desirability(75.0) - 0.25).abs() < 1e-9);
    }

    // upstream-015: importance-weighted overall desirability.

    #[test]
    fn default_importance_is_one() {
        let specs = vec![maximize(0.0, 100.0, 1.0), maximize(0.0, 100.0, 1.0)];
        assert!((specs[0].importance() - 1.0).abs() < 1e-12);
        let d = overall_desirability(&specs, &[50.0, 100.0]).unwrap();
        assert!((d - 0.5_f64.sqrt()).abs() < 1e-9, "D = {d}");
    }

    #[test]
    fn equal_desirabilities_give_that_value_regardless_of_weights() {
        // Both dᵢ = 0.5; any positive weights → weighted geometric mean = 0.5.
        let specs = vec![
            maximize(0.0, 10.0, 1.0).with_importance(10.0).unwrap(),
            maximize(0.0, 10.0, 1.0).with_importance(1.0).unwrap(),
        ];
        let d = overall_desirability(&specs, &[5.0, 5.0]).unwrap();
        assert!((d - 0.5).abs() < 1e-9, "equal d's must give 0.5, got {d}");
    }

    #[test]
    fn higher_importance_pulls_overall_toward_that_response() {
        // d1 = 0.25 (weighted heavily), d2 = 1.0. As r1 grows, D → 0.25.
        let base = maximize(0.0, 100.0, 1.0);
        let make = |r1: f64| {
            vec![
                base.clone().with_importance(r1).unwrap(),
                base.clone().with_importance(1.0).unwrap(),
            ]
        };
        let d_equal = overall_desirability(&make(1.0), &[25.0, 100.0]).unwrap();
        let d_heavy = overall_desirability(&make(9.0), &[25.0, 100.0]).unwrap();
        // Equal weights: sqrt(0.25*1) = 0.5. Heavy on response 1: (0.25^9)^(1/10) ≈ 0.28.
        assert!((d_equal - 0.5).abs() < 1e-9, "equal → 0.5, got {d_equal}");
        assert!(
            d_heavy < d_equal,
            "heavier weight on the worse response must lower D"
        );
        assert!((d_heavy - 0.25_f64.powf(0.9)).abs() < 1e-9, "D = {d_heavy}");
    }

    #[test]
    fn zero_desirability_still_vetoes_with_weights() {
        let specs = vec![
            maximize(50.0, 100.0, 1.0).with_importance(5.0).unwrap(),
            maximize(50.0, 100.0, 1.0).with_importance(1.0).unwrap(),
        ];
        let d = overall_desirability(&specs, &[80.0, 10.0]).unwrap(); // second → d=0
        assert!(
            d.abs() < 1e-12,
            "any zero d with positive weight → D=0, got {d}"
        );
    }

    #[test]
    fn zero_importance_leaves_a_response_out() {
        let specs = vec![
            maximize(50.0, 100.0, 1.0),
            maximize(50.0, 100.0, 1.0).with_importance(0.0).unwrap(),
        ];
        // The second response is below L (d = 0) but carries no weight.
        let d = overall_desirability(&specs, &[75.0, 10.0]).unwrap();
        assert!((d - 0.5).abs() < 1e-9, "D = {d}");
    }

    /// The six rows of the report: each spec whose ramp has no width or runs
    /// backwards used to score as a step or a constant.
    #[test]
    fn refuses_ramps_with_no_width_or_running_backwards() {
        let cases = [
            (GoalType::Maximize, 10.0, 0.0, 0.0),
            (GoalType::Maximize, 0.0, 0.0, 10.0),
            (GoalType::Minimize, 10.0, 10.0, 0.0),
            (GoalType::Target, 0.0, 15.0, 10.0),
            (GoalType::Target, 10.0, 5.0, 0.0),
            (GoalType::Target, 0.0, 0.0, 10.0),
            (GoalType::Target, 0.0, 10.0, 10.0),
            (GoalType::Minimize, 0.0, f64::NEG_INFINITY, 10.0),
            (GoalType::Maximize, f64::NAN, 10.0, 10.0),
        ];
        for (goal, lower, t, upper) in cases {
            match ResponseSpec::new(goal, lower, t, upper, 1.0, 1.0) {
                Err(DoeError::InvalidDesirabilityLimits {
                    index: None,
                    goal: g,
                    ..
                }) => assert_eq!(g, goal),
                other => panic!("{goal:?} {lower}/{t}/{upper}: {other:?}"),
            }
        }
    }

    #[test]
    fn a_limit_the_goal_does_not_use_is_not_examined() {
        // Maximize reads L and T; Minimize reads T and U.
        assert!(ResponseSpec::maximize(0.0, 10.0, 1.0).is_ok());
        assert!(ResponseSpec::minimize(0.0, 10.0, 1.0).is_ok());
        // Nor is the exponent it does not use.
        assert!(ResponseSpec::new(GoalType::Maximize, 0.0, 10.0, 10.0, 1.0, 0.0).is_ok());
    }

    #[test]
    fn refuses_exponents_and_weights_outside_their_range() {
        let positive = |parameter, got| DoeError::ParameterOutOfRange {
            parameter,
            index: None,
            min: None,
            max: None,
            exclusive_min: Some(0.0),
            exclusive_max: None,
            got,
        };
        for s in [0.0, -1.0] {
            assert_eq!(
                ResponseSpec::maximize(0.0, 10.0, s).map(|_| ()),
                Err(positive("s1", s))
            );
            assert_eq!(
                ResponseSpec::target(0.0, 5.0, 10.0, 1.0, s).map(|_| ()),
                Err(positive("s2", s))
            );
        }
        let not_finite = |parameter| DoeError::ValueNotFinite {
            parameter,
            index: None,
        };
        for s in [f64::INFINITY, f64::NAN] {
            assert_eq!(
                ResponseSpec::maximize(0.0, 10.0, s).map(|_| ()),
                Err(not_finite("s1"))
            );
            assert_eq!(
                ResponseSpec::target(0.0, 5.0, 10.0, 1.0, s).map(|_| ()),
                Err(not_finite("s2"))
            );
        }
        assert_eq!(
            maximize(0.0, 10.0, 1.0).with_importance(-1.0).map(|_| ()),
            Err(DoeError::ParameterOutOfRange {
                parameter: "importance",
                index: None,
                min: Some(0.0),
                max: None,
                exclusive_min: None,
                exclusive_max: None,
                got: -1.0,
            })
        );
        for r in [f64::INFINITY, f64::NAN] {
            assert_eq!(
                maximize(0.0, 10.0, 1.0).with_importance(r).map(|_| ()),
                Err(not_finite("importance"))
            );
        }
    }

    #[test]
    fn overall_refuses_what_it_cannot_combine() {
        let spec = maximize(0.0, 10.0, 1.0);
        assert_eq!(
            overall_desirability(&[], &[]),
            Err(DoeError::EmptyInput { parameter: "specs" })
        );
        assert_eq!(
            overall_desirability(std::slice::from_ref(&spec), &[1.0, 2.0]),
            Err(DoeError::DimensionMismatch {
                parameter: "responses",
                index: None,
                expected: 1,
                got: 2,
                counted: "one response per specification",
            })
        );
        let unweighted = spec.with_importance(0.0).unwrap();
        assert_eq!(
            overall_desirability(&[unweighted], &[5.0]),
            Err(DoeError::NoWeightedResponse)
        );
    }
}

/// The individual and overall desirabilities against the Derringer-Suich
/// formulas in NumPy: a maximize (importance 2), a minimize (1) and a target
/// goal with unequal exponents (3), six points -- inside the ramps, on the
/// limits, and outside, where one vetoes the whole design. And a target
/// specification in units of 1e-13, where a fixed tolerance around the target
/// would call every value on target.
#[cfg(test)]
mod numpy_reference {
    use super::*;

    fn close(got: f64, want: f64) -> bool {
        (got - want).abs() <= 1e-12 * want.abs().max(1.0)
    }

    #[test]
    fn individual_and_overall_desirability_reproduce_the_formulas() {
        let specs = [
            ResponseSpec::maximize(20.0, 80.0, 1.5)
                .and_then(|s| s.with_importance(2.0))
                .expect("valid"),
            ResponseSpec::minimize(5.0, 25.0, 0.7).expect("valid"),
            ResponseSpec::target(40.0, 50.0, 60.0, 2.0, 0.5)
                .and_then(|s| s.with_importance(3.0))
                .expect("valid"),
        ];
        let cases: [([f64; 3], [f64; 3], f64); 6] = [
            (
                [55.0, 12.0, 47.0],
                [0.445528192565151, 0.739672290708735, 0.49],
                0.508428236549523,
            ),
            (
                [79.9, 5.0, 50.0],
                [0.9975010419562, 1.0, 1.0],
                0.999166319154791,
            ),
            (
                [30.0, 10.0, 45.0],
                [0.0680413817439772, 0.817603768177013, 0.25],
                0.197386839456212,
            ),
            (
                [85.0, 24.0, 58.0],
                [1.0, 0.122822802611579, 0.447213595499958],
                0.471488024155257,
            ),
            ([60.0, 26.0, 50.0], [0.544331053951817, 0.0, 1.0], 0.0),
            (
                [50.0, 15.0, 41.5],
                [0.353553390593274, 0.615572206672458, 0.0225],
                0.0978264052324173,
            ),
        ];
        for (y, d, overall) in cases {
            for ((spec, y), want) in specs.iter().zip(y).zip(d) {
                assert!(
                    close(spec.desirability(y), want),
                    "d({y}): got {}, want {want}",
                    spec.desirability(y)
                );
            }
            let got = overall_desirability(&specs, &y).expect("one response per specification");
            assert!(close(got, overall), "{y:?}: got {got}, want {overall}");
        }
    }

    #[test]
    fn a_target_specification_does_not_depend_on_the_unit() {
        let spec = ResponseSpec::target(0.0, 4e-13, 1e-12, 1.0, 1.0).expect("valid");
        assert!(close(spec.desirability(1e-13), 0.25));
        assert!(close(spec.desirability(7e-13), 0.5));
        assert_eq!(spec.desirability(4e-13), 1.0);
    }

    #[test]
    fn the_search_returns_the_candidate_with_the_largest_overall_desirability() {
        let specs = [
            ResponseSpec::maximize(20.0, 80.0, 1.5)
                .and_then(|s| s.with_importance(2.0))
                .expect("valid"),
            ResponseSpec::minimize(5.0, 25.0, 0.7).expect("valid"),
            ResponseSpec::target(40.0, 50.0, 60.0, 2.0, 0.5)
                .and_then(|s| s.with_importance(3.0))
                .expect("valid"),
        ];
        // The six points of the first test; the second has the largest D.
        let candidates = vec![
            vec![55.0, 12.0, 47.0],
            vec![79.9, 5.0, 50.0],
            vec![30.0, 10.0, 45.0],
            vec![85.0, 24.0, 58.0],
            vec![60.0, 26.0, 50.0],
            vec![50.0, 15.0, 41.5],
        ];
        match optimize_desirability(&specs, &candidates).expect("valid") {
            DesirabilityOptimum::Best {
                index,
                overall,
                individual,
            } => {
                assert_eq!(index, 1);
                assert!(close(overall, 0.999166319154791));
                for (got, want) in individual.iter().zip([0.9975010419562, 1.0, 1.0]) {
                    assert!(close(*got, want));
                }
            }
            other => panic!("expected a best candidate, got {other:?}"),
        }
    }
}
