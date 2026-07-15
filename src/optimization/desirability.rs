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

/// Optimization goal for a single response.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GoalType {
    /// Maximize the response (larger is better).
    Maximize,
    /// Minimize the response (smaller is better).
    Minimize,
    /// Hit a target value.
    Target,
}

/// Specification for a single response's desirability function.
#[derive(Debug, Clone)]
pub struct ResponseSpec {
    /// Optimization goal.
    pub goal: GoalType,
    /// Lower acceptability limit (L). Below this: d = 0 for Maximize/Target.
    pub lower: f64,
    /// Target value (T). At this value: d = 1.
    pub target: f64,
    /// Upper acceptability limit (U). Above this: d = 0 for Minimize/Target.
    pub upper: f64,
    /// Curve-shape exponent for the ascending side (left of target).
    pub s1: f64,
    /// Curve-shape exponent for the descending side (right of target). Equals s1 for Maximize/Minimize.
    pub s2: f64,
    /// Relative importance weight rᵢ for this response in the overall desirability
    /// (Derringer & Suich 1980). Defaults to 1.0 (equal weighting). A larger value
    /// makes this response count more toward `overall_desirability`; it is a
    /// *weight*, distinct from the curve-shape exponents `s1`/`s2`.
    pub importance: f64,
}

impl ResponseSpec {
    /// Create a "larger-is-better" specification.
    ///
    /// d = 0 for y < lower; d = 1 for y ≥ target (= upper); linear ramp in between.
    ///
    /// # Examples
    ///
    /// ```
    /// use u_doe::optimization::desirability::ResponseSpec;
    /// let spec = ResponseSpec::maximize(50.0, 100.0, 100.0, 1.0);
    /// assert!((spec.desirability(75.0) - 0.5).abs() < 1e-9);
    /// ```
    pub fn maximize(lower: f64, target: f64, upper: f64, s: f64) -> Self {
        Self {
            goal: GoalType::Maximize,
            lower,
            target,
            upper,
            s1: s,
            s2: s,
            importance: 1.0,
        }
    }

    /// Create a "smaller-is-better" specification.
    ///
    /// d = 1 for y ≤ target (= lower); d = 0 for y > upper; linear ramp in between.
    ///
    /// # Examples
    ///
    /// ```
    /// use u_doe::optimization::desirability::ResponseSpec;
    /// let spec = ResponseSpec::minimize(0.0, 0.0, 50.0, 1.0);
    /// assert!((spec.desirability(25.0) - 0.5).abs() < 1e-9);
    /// ```
    pub fn minimize(lower: f64, target: f64, upper: f64, s: f64) -> Self {
        Self {
            goal: GoalType::Minimize,
            lower,
            target,
            upper,
            s1: s,
            s2: s,
            importance: 1.0,
        }
    }

    /// Create a "hit-a-target" specification.
    ///
    /// d = 1 at y = target; d = 0 outside [lower, upper].
    /// s1 controls the left ramp shape; s2 controls the right ramp.
    ///
    /// # Examples
    ///
    /// ```
    /// use u_doe::optimization::desirability::ResponseSpec;
    /// let spec = ResponseSpec::target(60.0, 80.0, 100.0, 1.0, 1.0);
    /// assert!((spec.desirability(80.0) - 1.0).abs() < 1e-10);
    /// ```
    pub fn target(lower: f64, target: f64, upper: f64, s1: f64, s2: f64) -> Self {
        Self {
            goal: GoalType::Target,
            lower,
            target,
            upper,
            s1,
            s2,
            importance: 1.0,
        }
    }

    /// Set this response's relative importance weight rᵢ (Derringer & Suich 1980),
    /// used by [`overall_desirability`]. Builder-style; defaults to 1.0.
    ///
    /// # Examples
    ///
    /// ```
    /// use u_doe::optimization::desirability::ResponseSpec;
    /// let spec = ResponseSpec::maximize(0.0, 100.0, 100.0, 1.0).with_importance(3.0);
    /// assert!((spec.importance - 3.0).abs() < 1e-12);
    /// ```
    pub fn with_importance(mut self, importance: f64) -> Self {
        self.importance = importance;
        self
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
/// where rᵢ = `spec.importance` (default 1.0). With all weights equal this reduces
/// to the plain geometric mean (∏ dᵢ)^(1/m).
///
/// Returns 0.0 if any individual desirability is 0 (with positive weight), if the
/// inputs are empty / length-mismatched, or if the weights sum to ≤ 0. Responses
/// with a non-positive `importance` are excluded from the product.
///
/// # Examples
///
/// ```
/// use u_doe::optimization::desirability::{ResponseSpec, overall_desirability};
/// let specs = vec![
///     ResponseSpec::maximize(0.0, 100.0, 100.0, 1.0),
///     ResponseSpec::maximize(0.0, 100.0, 100.0, 1.0),
/// ];
/// let d = overall_desirability(&specs, &[100.0, 100.0]);
/// assert!((d - 1.0).abs() < 1e-10);
/// ```
pub fn overall_desirability(specs: &[ResponseSpec], responses: &[f64]) -> f64 {
    if specs.is_empty() || specs.len() != responses.len() {
        return 0.0;
    }
    let sum_r: f64 = specs.iter().map(|spec| spec.importance.max(0.0)).sum();
    if sum_r <= 0.0 {
        return 0.0;
    }
    // Weighted geometric mean via logs: D = exp( (Σ rᵢ·ln dᵢ) / Σrᵢ ). A dᵢ = 0
    // with positive weight vetoes the whole design (D = 0), matching the
    // unweighted convention; a non-positive weight drops the response entirely.
    let mut weighted_log_sum = 0.0;
    for (spec, &y) in specs.iter().zip(responses.iter()) {
        let r = spec.importance;
        if r <= 0.0 {
            continue;
        }
        let d = spec.desirability(y);
        if d <= 0.0 {
            return 0.0;
        }
        weighted_log_sum += r * d.ln();
    }
    (weighted_log_sum / sum_r).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn maximize_at_target_is_one() {
        let spec = ResponseSpec::maximize(50.0, 100.0, 100.0, 1.0);
        assert!((spec.desirability(100.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn maximize_below_lower_is_zero() {
        let spec = ResponseSpec::maximize(50.0, 100.0, 100.0, 1.0);
        assert!(spec.desirability(49.9).abs() < 1e-10);
    }

    #[test]
    fn maximize_midpoint_linear() {
        // s=1, linear: d(75) = (75-50)/(100-50) = 0.5
        let spec = ResponseSpec::maximize(50.0, 100.0, 100.0, 1.0);
        assert!((spec.desirability(75.0) - 0.5).abs() < 1e-9);
    }

    #[test]
    fn minimize_at_target_is_one() {
        let spec = ResponseSpec::minimize(0.0, 0.0, 50.0, 1.0);
        assert!((spec.desirability(0.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn minimize_above_upper_is_zero() {
        let spec = ResponseSpec::minimize(0.0, 0.0, 50.0, 1.0);
        assert!(spec.desirability(50.1).abs() < 1e-10);
    }

    #[test]
    fn target_at_target_is_one() {
        let spec = ResponseSpec::target(60.0, 80.0, 100.0, 1.0, 1.0);
        assert!((spec.desirability(80.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn target_below_lower_is_zero() {
        let spec = ResponseSpec::target(60.0, 80.0, 100.0, 1.0, 1.0);
        assert!(spec.desirability(59.9).abs() < 1e-10);
    }

    #[test]
    fn target_above_upper_is_zero() {
        let spec = ResponseSpec::target(60.0, 80.0, 100.0, 1.0, 1.0);
        assert!(spec.desirability(100.1).abs() < 1e-10);
    }

    #[test]
    fn target_left_midpoint() {
        // s1=1, linear: d(70) = (70-60)/(80-60) = 0.5
        let spec = ResponseSpec::target(60.0, 80.0, 100.0, 1.0, 1.0);
        assert!((spec.desirability(70.0) - 0.5).abs() < 1e-9);
    }

    #[test]
    fn target_right_midpoint() {
        // s2=1: d(90) = (100-90)/(100-80) = 0.5
        let spec = ResponseSpec::target(60.0, 80.0, 100.0, 1.0, 1.0);
        assert!((spec.desirability(90.0) - 0.5).abs() < 1e-9);
    }

    #[test]
    fn overall_geometric_mean() {
        let specs = vec![
            ResponseSpec::maximize(0.0, 100.0, 100.0, 1.0),
            ResponseSpec::maximize(0.0, 100.0, 100.0, 1.0),
        ];
        // d1=0.5 (at 50), d2=1.0 (at 100) → D = (0.5 * 1.0)^(1/2) = sqrt(0.5) ≈ 0.707
        let d = overall_desirability(&specs, &[50.0, 100.0]);
        assert!((d - 0.5_f64.sqrt()).abs() < 1e-9);
    }

    #[test]
    fn overall_zero_if_any_zero() {
        let specs = vec![
            ResponseSpec::maximize(50.0, 100.0, 100.0, 1.0),
            ResponseSpec::maximize(50.0, 100.0, 100.0, 1.0),
        ];
        // Second response below L → d=0 → D=0
        let d = overall_desirability(&specs, &[80.0, 10.0]);
        assert!(d.abs() < 1e-10);
    }

    #[test]
    fn desirability_weight_s2_convex() {
        // s=2 (convex): d(75) = ((75-50)/50)^2 = 0.25
        let spec = ResponseSpec::maximize(50.0, 100.0, 100.0, 2.0);
        assert!((spec.desirability(75.0) - 0.25).abs() < 1e-9);
    }

    // upstream-015: importance-weighted overall desirability.

    #[test]
    fn default_importance_is_one_and_backward_compatible() {
        // Default importance = 1.0 → unchanged plain geometric mean.
        let specs = vec![
            ResponseSpec::maximize(0.0, 100.0, 100.0, 1.0),
            ResponseSpec::maximize(0.0, 100.0, 100.0, 1.0),
        ];
        assert!((specs[0].importance - 1.0).abs() < 1e-12);
        let d = overall_desirability(&specs, &[50.0, 100.0]);
        assert!((d - 0.5_f64.sqrt()).abs() < 1e-9, "D = {d}");
    }

    #[test]
    fn equal_desirabilities_give_that_value_regardless_of_weights() {
        // Both dᵢ = 0.5; any positive weights → weighted geometric mean = 0.5.
        // (The bug symptom in the report: r=[10,1] must NOT collapse this to ~0.)
        let specs = vec![
            ResponseSpec::maximize(0.0, 10.0, 10.0, 1.0).with_importance(10.0),
            ResponseSpec::maximize(0.0, 10.0, 10.0, 1.0).with_importance(1.0),
        ];
        let d = overall_desirability(&specs, &[5.0, 5.0]);
        assert!((d - 0.5).abs() < 1e-9, "equal d's must give 0.5, got {d}");
    }

    #[test]
    fn higher_importance_pulls_overall_toward_that_response() {
        // d1 = 0.25 (weighted heavily), d2 = 1.0. As r1 grows, D → 0.25.
        let base = ResponseSpec::maximize(0.0, 100.0, 100.0, 1.0);
        let make = |r1: f64| {
            vec![
                base.clone().with_importance(r1),
                base.clone().with_importance(1.0),
            ]
        };
        let d_equal = overall_desirability(&make(1.0), &[25.0, 100.0]);
        let d_heavy = overall_desirability(&make(9.0), &[25.0, 100.0]);
        // Equal weights: sqrt(0.25*1) = 0.5. Heavy on response 1: (0.25^9)^(1/10) ≈ 0.28.
        assert!((d_equal - 0.5).abs() < 1e-9, "equal → 0.5, got {d_equal}");
        assert!(d_heavy < d_equal, "heavier weight on the worse response must lower D");
        assert!((d_heavy - 0.25_f64.powf(0.9)).abs() < 1e-9, "D = {d_heavy}");
    }

    #[test]
    fn zero_desirability_still_vetoes_with_weights() {
        let specs = vec![
            ResponseSpec::maximize(50.0, 100.0, 100.0, 1.0).with_importance(5.0),
            ResponseSpec::maximize(50.0, 100.0, 100.0, 1.0).with_importance(1.0),
        ];
        let d = overall_desirability(&specs, &[80.0, 10.0]); // second → d=0
        assert!(d.abs() < 1e-12, "any zero d with positive weight → D=0, got {d}");
    }
}
