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

use crate::error::DoeError;

const POSITIVE: &str = "finite and greater than 0";

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
    /// [`DoeError::ValueOutOfDomain`] unless each exponent the goal
    /// uses is finite and positive: `s = 0` makes every point of the ramp
    /// d = 1, and a negative `s` sends d above 1.
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
        let exponent_ok = |s: f64| s.is_finite() && s > 0.0;
        if !exponent_ok(s1) {
            return Err(DoeError::ValueOutOfDomain {
                index: None,
                parameter: "s1",
                value: s1,
                domain: POSITIVE,
            });
        }
        if goal == GoalType::Target && !exponent_ok(s2) {
            return Err(DoeError::ValueOutOfDomain {
                index: None,
                parameter: "s2",
                value: s2,
                domain: POSITIVE,
            });
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
    /// between. `upper` is not used.
    ///
    /// # Errors
    /// As [`ResponseSpec::new`].
    ///
    /// # Examples
    ///
    /// ```
    /// use u_doe::optimization::desirability::ResponseSpec;
    /// let spec = ResponseSpec::maximize(50.0, 100.0, 100.0, 1.0).unwrap();
    /// assert!((spec.desirability(75.0) - 0.5).abs() < 1e-9);
    /// ```
    pub fn maximize(lower: f64, target: f64, upper: f64, s: f64) -> Result<Self, DoeError> {
        Self::new(GoalType::Maximize, lower, target, upper, s, s)
    }

    /// Create a "smaller-is-better" specification.
    ///
    /// d = 1 for y ≤ target; d = 0 for y > upper; a ramp with exponent `s` in
    /// between. `lower` is not used.
    ///
    /// # Errors
    /// As [`ResponseSpec::new`].
    ///
    /// # Examples
    ///
    /// ```
    /// use u_doe::optimization::desirability::ResponseSpec;
    /// let spec = ResponseSpec::minimize(0.0, 0.0, 50.0, 1.0).unwrap();
    /// assert!((spec.desirability(25.0) - 0.5).abs() < 1e-9);
    /// ```
    pub fn minimize(lower: f64, target: f64, upper: f64, s: f64) -> Result<Self, DoeError> {
        Self::new(GoalType::Minimize, lower, target, upper, s, s)
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
    /// [`DoeError::ValueOutOfDomain`] unless `importance` is finite
    /// and not negative.
    ///
    /// # Examples
    ///
    /// ```
    /// use u_doe::optimization::desirability::ResponseSpec;
    /// let spec = ResponseSpec::maximize(0.0, 100.0, 100.0, 1.0)
    ///     .and_then(|s| s.with_importance(3.0))
    ///     .unwrap();
    /// assert!((spec.importance() - 3.0).abs() < 1e-12);
    /// ```
    pub fn with_importance(mut self, importance: f64) -> Result<Self, DoeError> {
        if !(importance.is_finite() && importance >= 0.0) {
            return Err(DoeError::ValueOutOfDomain {
                index: None,
                parameter: "importance",
                value: importance,
                domain: "finite and not negative",
            });
        }
        self.importance = importance;
        Ok(self)
    }

    /// Optimization goal.
    pub fn goal(&self) -> GoalType {
        self.goal
    }

    /// Lower acceptability limit (L).
    pub fn lower(&self) -> f64 {
        self.lower
    }

    /// Target value (T).
    pub fn target_value(&self) -> f64 {
        self.target
    }

    /// Upper acceptability limit (U).
    pub fn upper(&self) -> f64 {
        self.upper
    }

    /// Exponent of the rising ramp (the falling one for `Minimize`).
    pub fn s1(&self) -> f64 {
        self.s1
    }

    /// Exponent of `Target`'s falling ramp.
    pub fn s2(&self) -> f64 {
        self.s2
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
                } else if (y - self.target).abs() < 1e-12 {
                    1.0
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
/// [`DoeError::EmptyResponses`] if there are no specifications;
/// [`DoeError::ResponseCountMismatch`] unless there is one response per
/// specification; [`DoeError::NoWeightedResponse`] if every importance is 0.
///
/// # Examples
///
/// ```
/// use u_doe::optimization::desirability::{ResponseSpec, overall_desirability};
/// let specs = vec![
///     ResponseSpec::maximize(0.0, 100.0, 100.0, 1.0).unwrap(),
///     ResponseSpec::maximize(0.0, 100.0, 100.0, 1.0).unwrap(),
/// ];
/// let d = overall_desirability(&specs, &[100.0, 100.0]).unwrap();
/// assert!((d - 1.0).abs() < 1e-10);
/// ```
pub fn overall_desirability(specs: &[ResponseSpec], responses: &[f64]) -> Result<f64, DoeError> {
    if specs.is_empty() {
        return Err(DoeError::EmptyResponses);
    }
    if specs.len() != responses.len() {
        return Err(DoeError::ResponseCountMismatch {
            expected: specs.len(),
            got: responses.len(),
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

#[cfg(test)]
mod tests {
    use super::*;

    fn maximize(lower: f64, target: f64, upper: f64, s: f64) -> ResponseSpec {
        ResponseSpec::maximize(lower, target, upper, s).expect("valid maximize spec")
    }

    fn minimize(lower: f64, target: f64, upper: f64, s: f64) -> ResponseSpec {
        ResponseSpec::minimize(lower, target, upper, s).expect("valid minimize spec")
    }

    fn target(lower: f64, t: f64, upper: f64, s1: f64, s2: f64) -> ResponseSpec {
        ResponseSpec::target(lower, t, upper, s1, s2).expect("valid target spec")
    }

    #[test]
    fn maximize_at_target_is_one() {
        let spec = maximize(50.0, 100.0, 100.0, 1.0);
        assert!((spec.desirability(100.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn maximize_below_lower_is_zero() {
        let spec = maximize(50.0, 100.0, 100.0, 1.0);
        assert!(spec.desirability(49.9).abs() < 1e-10);
    }

    #[test]
    fn maximize_midpoint_linear() {
        // s=1, linear: d(75) = (75-50)/(100-50) = 0.5
        let spec = maximize(50.0, 100.0, 100.0, 1.0);
        assert!((spec.desirability(75.0) - 0.5).abs() < 1e-9);
    }

    #[test]
    fn minimize_at_target_is_one() {
        let spec = minimize(0.0, 0.0, 50.0, 1.0);
        assert!((spec.desirability(0.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn minimize_above_upper_is_zero() {
        let spec = minimize(0.0, 0.0, 50.0, 1.0);
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
        let specs = vec![
            maximize(0.0, 100.0, 100.0, 1.0),
            maximize(0.0, 100.0, 100.0, 1.0),
        ];
        // d1=0.5 (at 50), d2=1.0 (at 100) → D = (0.5 * 1.0)^(1/2) = sqrt(0.5) ≈ 0.707
        let d = overall_desirability(&specs, &[50.0, 100.0]).unwrap();
        assert!((d - 0.5_f64.sqrt()).abs() < 1e-9);
    }

    #[test]
    fn overall_zero_if_any_zero() {
        let specs = vec![
            maximize(50.0, 100.0, 100.0, 1.0),
            maximize(50.0, 100.0, 100.0, 1.0),
        ];
        // Second response below L → d=0 → D=0
        let d = overall_desirability(&specs, &[80.0, 10.0]).unwrap();
        assert!(d.abs() < 1e-10);
    }

    #[test]
    fn desirability_weight_s2_convex() {
        // s=2 (convex): d(75) = ((75-50)/50)^2 = 0.25
        let spec = maximize(50.0, 100.0, 100.0, 2.0);
        assert!((spec.desirability(75.0) - 0.25).abs() < 1e-9);
    }

    // upstream-015: importance-weighted overall desirability.

    #[test]
    fn default_importance_is_one() {
        let specs = vec![
            maximize(0.0, 100.0, 100.0, 1.0),
            maximize(0.0, 100.0, 100.0, 1.0),
        ];
        assert!((specs[0].importance() - 1.0).abs() < 1e-12);
        let d = overall_desirability(&specs, &[50.0, 100.0]).unwrap();
        assert!((d - 0.5_f64.sqrt()).abs() < 1e-9, "D = {d}");
    }

    #[test]
    fn equal_desirabilities_give_that_value_regardless_of_weights() {
        // Both dᵢ = 0.5; any positive weights → weighted geometric mean = 0.5.
        let specs = vec![
            maximize(0.0, 10.0, 10.0, 1.0)
                .with_importance(10.0)
                .unwrap(),
            maximize(0.0, 10.0, 10.0, 1.0).with_importance(1.0).unwrap(),
        ];
        let d = overall_desirability(&specs, &[5.0, 5.0]).unwrap();
        assert!((d - 0.5).abs() < 1e-9, "equal d's must give 0.5, got {d}");
    }

    #[test]
    fn higher_importance_pulls_overall_toward_that_response() {
        // d1 = 0.25 (weighted heavily), d2 = 1.0. As r1 grows, D → 0.25.
        let base = maximize(0.0, 100.0, 100.0, 1.0);
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
            maximize(50.0, 100.0, 100.0, 1.0)
                .with_importance(5.0)
                .unwrap(),
            maximize(50.0, 100.0, 100.0, 1.0)
                .with_importance(1.0)
                .unwrap(),
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
            maximize(50.0, 100.0, 100.0, 1.0),
            maximize(50.0, 100.0, 100.0, 1.0)
                .with_importance(0.0)
                .unwrap(),
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
        assert!(ResponseSpec::maximize(0.0, 10.0, -5.0, 1.0).is_ok());
        assert!(ResponseSpec::minimize(99.0, 0.0, 10.0, 1.0).is_ok());
        // Nor is the exponent it does not use.
        assert!(ResponseSpec::new(GoalType::Maximize, 0.0, 10.0, 10.0, 1.0, 0.0).is_ok());
    }

    #[test]
    fn refuses_exponents_and_weights_outside_their_domain() {
        for s in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(
                matches!(
                    ResponseSpec::maximize(0.0, 10.0, 10.0, s),
                    Err(DoeError::ValueOutOfDomain {
                        parameter: "s1",
                        ..
                    })
                ),
                "s1 = {s}"
            );
            assert!(
                matches!(
                    ResponseSpec::target(0.0, 5.0, 10.0, 1.0, s),
                    Err(DoeError::ValueOutOfDomain {
                        parameter: "s2",
                        ..
                    })
                ),
                "s2 = {s}"
            );
        }
        for r in [-1.0, f64::INFINITY, f64::NAN] {
            assert!(
                matches!(
                    maximize(0.0, 10.0, 10.0, 1.0).with_importance(r),
                    Err(DoeError::ValueOutOfDomain {
                        parameter: "importance",
                        ..
                    })
                ),
                "importance = {r}"
            );
        }
    }

    #[test]
    fn overall_refuses_what_it_cannot_combine() {
        let spec = maximize(0.0, 10.0, 10.0, 1.0);
        assert_eq!(
            overall_desirability(&[], &[]),
            Err(DoeError::EmptyResponses)
        );
        assert_eq!(
            overall_desirability(std::slice::from_ref(&spec), &[1.0, 2.0]),
            Err(DoeError::ResponseCountMismatch {
                expected: 1,
                got: 2
            })
        );
        let unweighted = spec.with_importance(0.0).unwrap();
        assert_eq!(
            overall_desirability(&[unweighted], &[5.0]),
            Err(DoeError::NoWeightedResponse)
        );
    }
}
