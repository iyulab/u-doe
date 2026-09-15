/// Errors produced by u-doe operations.
///
/// Every variant carries the values that explain it as fields rather than as
/// text, so a caller can tell refusals apart and say what to change without
/// reading the message. With the `wasm` feature the variants serialize as a
/// flat object whose `code` is the variant name in `snake_case`
/// (`OverSpecifiedModel` → `"over_specified_model"`) next to those fields; the
/// WebAssembly bindings attach that object to the `Error` they throw.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(
    feature = "wasm",
    derive(serde::Serialize),
    serde(tag = "code", rename_all = "snake_case")
)]
pub enum DoeError {
    /// Invalid number of factors for this design type.
    InvalidFactorCount { min: usize, max: usize, got: usize },
    /// A numeric parameter other than the factor count is outside the range the
    /// function accepts. `max` is `None` when there is no upper bound.
    ParameterOutOfRange {
        parameter: &'static str,
        min: usize,
        max: Option<usize>,
        got: usize,
    },
    /// The 2^(k-p) fraction is not in the standard generator table. `supported`
    /// lists every `(k, p)` that is.
    UnsupportedFraction {
        k: usize,
        p: usize,
        supported: Vec<(usize, usize)>,
    },
    /// No Taguchi orthogonal array has this name. `supported` lists the names.
    UnknownArray {
        array: String,
        supported: &'static [&'static str],
    },
    /// A design with no runs or no factors: there is nothing to analyse.
    EmptyDesign { runs: usize, factors: usize },
    /// A run has a different number of entries from the number of factor names.
    DesignShapeMismatch {
        run: usize,
        expected: usize,
        got: usize,
    },
    /// The number of responses differs from the number of runs.
    ResponseCountMismatch { expected: usize, got: usize },
    /// Effect name not found among the effects estimable for the design.
    UnknownEffect { effect: String },
    /// An effect refers to a factor column the design does not have.
    EffectOutsideDesign { effect: String, factors: usize },
    /// A design matrix that must be two-level coded (every entry -1 or +1)
    /// carries some other value.
    NotTwoLevelCoded {
        run: usize,
        factor: usize,
        value: f64,
    },
    /// Two requested effects share the same contrast column, so they are
    /// aliased in this design and cannot be estimated separately. `second` is
    /// `I` when the contrast is constant, i.e. aliased with the overall mean.
    AliasedEffects { first: String, second: String },
    /// Two requested effects -- or an effect and the overall mean, written `I` --
    /// have contrasts that are correlated without being identical. The per-term
    /// sums of squares then overlap, so they do not add up to the model sum of
    /// squares and no ANOVA can be built from them. Plackett-Burman designs
    /// produce this for two-factor interactions.
    PartiallyAliasedEffects { first: String, second: String },
    /// The requested model has more terms than the design has degrees of
    /// freedom to spend on them. `terms` counts the model terms (and the
    /// curvature term, when there is one) but not the mean; the design offers
    /// `runs - 1` degrees of freedom.
    OverSpecifiedModel { terms: usize, runs: usize },
    /// The model matrix is singular: its columns are linearly dependent over
    /// the runs given, so the coefficients are not identifiable.
    SingularModel,
    /// A response-surface model whose coefficient list does not have the
    /// length a full quadratic model in `factors` factors has.
    CoefficientCountMismatch {
        factors: usize,
        expected: usize,
        got: usize,
    },
    /// Lenth's method needs at least `needed` distinct contrasts.
    TooFewContrasts { needed: usize, got: usize },
    /// A coding range whose low end is not below its high end.
    InvalidCodingRange { low: f64, high: f64 },
    /// No responses were given.
    EmptyResponses,
    /// A run has fewer replicate measurements than the calculation needs.
    TooFewReplicates {
        run: usize,
        needed: usize,
        got: usize,
    },
    /// A response that must be strictly positive is not.
    NonPositiveResponse { run: usize, value: f64 },
    /// A desirability specification whose ramp has no width or runs
    /// backwards: the limits `goal` uses are not finite, or not ordered
    /// `lower < target` (rising) and `target < upper` (falling). `index` is the
    /// position of the specification in its list, when there is one.
    InvalidDesirabilityLimits {
        index: Option<usize>,
        goal: crate::optimization::desirability::GoalType,
        lower: f64,
        target: f64,
        upper: f64,
    },
    /// A desirability exponent (`s1`, `s2`) that is not finite and positive, or
    /// an `importance` weight that is not finite and non-negative.
    InvalidDesirabilityParameter {
        index: Option<usize>,
        parameter: &'static str,
        value: f64,
    },
    /// Every desirability specification has importance 0, so there is nothing
    /// to combine.
    NoWeightedResponse,
}

impl DoeError {
    /// The same error, placed at position `index` of the input list it came
    /// from. Only errors about one element of a list carry an index; others
    /// are returned unchanged.
    pub fn at_index(self, index: usize) -> Self {
        match self {
            DoeError::InvalidDesirabilityLimits {
                goal,
                lower,
                target,
                upper,
                ..
            } => DoeError::InvalidDesirabilityLimits {
                index: Some(index),
                goal,
                lower,
                target,
                upper,
            },
            DoeError::InvalidDesirabilityParameter {
                parameter, value, ..
            } => DoeError::InvalidDesirabilityParameter {
                index: Some(index),
                parameter,
                value,
            },
            other => other,
        }
    }

    /// Stable, machine-readable reason: the variant name in `snake_case`.
    ///
    /// ```
    /// use u_doe::error::DoeError;
    /// assert_eq!(DoeError::SingularModel.code(), "singular_model");
    /// ```
    pub fn code(&self) -> &'static str {
        match self {
            DoeError::InvalidFactorCount { .. } => "invalid_factor_count",
            DoeError::ParameterOutOfRange { .. } => "parameter_out_of_range",
            DoeError::UnsupportedFraction { .. } => "unsupported_fraction",
            DoeError::UnknownArray { .. } => "unknown_array",
            DoeError::EmptyDesign { .. } => "empty_design",
            DoeError::DesignShapeMismatch { .. } => "design_shape_mismatch",
            DoeError::ResponseCountMismatch { .. } => "response_count_mismatch",
            DoeError::UnknownEffect { .. } => "unknown_effect",
            DoeError::EffectOutsideDesign { .. } => "effect_outside_design",
            DoeError::NotTwoLevelCoded { .. } => "not_two_level_coded",
            DoeError::AliasedEffects { .. } => "aliased_effects",
            DoeError::PartiallyAliasedEffects { .. } => "partially_aliased_effects",
            DoeError::OverSpecifiedModel { .. } => "over_specified_model",
            DoeError::SingularModel => "singular_model",
            DoeError::CoefficientCountMismatch { .. } => "coefficient_count_mismatch",
            DoeError::TooFewContrasts { .. } => "too_few_contrasts",
            DoeError::InvalidCodingRange { .. } => "invalid_coding_range",
            DoeError::EmptyResponses => "empty_responses",
            DoeError::TooFewReplicates { .. } => "too_few_replicates",
            DoeError::NonPositiveResponse { .. } => "non_positive_response",
            DoeError::InvalidDesirabilityLimits { .. } => "invalid_desirability_limits",
            DoeError::InvalidDesirabilityParameter { .. } => "invalid_desirability_parameter",
            DoeError::NoWeightedResponse => "no_weighted_response",
        }
    }
}

impl std::fmt::Display for DoeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DoeError::InvalidFactorCount { min, max, got } => {
                write!(f, "invalid factor count: expected {min}..={max}, got {got}")
            }
            DoeError::ParameterOutOfRange {
                parameter,
                min,
                max: Some(max),
                got,
            } => write!(f, "{parameter} must be in {min}..={max}, got {got}"),
            DoeError::ParameterOutOfRange {
                parameter,
                min,
                max: None,
                got,
            } => write!(f, "{parameter} must be at least {min}, got {got}"),
            DoeError::UnsupportedFraction { k, p, supported } => {
                write!(
                    f,
                    "unsupported design: 2^({k}-{p}) not in standard table; supported (k,p): "
                )?;
                for (i, (k, p)) in supported.iter().enumerate() {
                    if i > 0 {
                        f.write_str(",")?;
                    }
                    write!(f, "({k},{p})")?;
                }
                Ok(())
            }
            DoeError::UnknownArray { array, supported } => write!(
                f,
                "unknown Taguchi array '{array}'; supported: {}",
                supported.join(", ")
            ),
            DoeError::EmptyDesign { runs, factors } => write!(
                f,
                "design has {runs} runs and {factors} factors; both must be non-zero"
            ),
            DoeError::DesignShapeMismatch { run, expected, got } => write!(
                f,
                "run {run} has {got} entries but the design names {expected} factors"
            ),
            DoeError::ResponseCountMismatch { expected, got } => {
                write!(f, "got {got} responses where {expected} were expected")
            }
            DoeError::UnknownEffect { effect } => write!(
                f,
                "unknown effect '{effect}': effect names are factor names joined with ':' \
                 (e.g. \"A\", \"A:B\")"
            ),
            DoeError::EffectOutsideDesign { effect, factors } => write!(
                f,
                "effect {effect} names a column the design does not have ({factors} factors)"
            ),
            DoeError::NotTwoLevelCoded { run, factor, value } => write!(
                f,
                "design is not two-level coded: run {run}, factor {factor} is {value}, \
                 expected -1 or +1. Effects are estimated from contrasts formed by \
                 multiplying factor columns, which is only defined for two-level runs. \
                 Centre points (every factor 0) are accepted by \
                 `analysis::anova::doe_anova`, which tests them for curvature; a design \
                 with axial or three-level points (central composite, Box-Behnken, \
                 definitive screening) must be fitted with `analysis::rsm::fit_rsm`"
            ),
            DoeError::AliasedEffects { first, second } if second == "I" => write!(
                f,
                "effect '{first}' is aliased with the overall mean: its contrast is the \
                 same in every run, so it would estimate the mean rather than an effect; \
                 lower the interaction order or drop the term"
            ),
            DoeError::AliasedEffects { first, second } => write!(
                f,
                "effects '{first}' and '{second}' share one contrast column and are \
                 therefore aliased in this design: they cannot both enter the model, \
                 because their sum of squares would be counted twice"
            ),
            DoeError::PartiallyAliasedEffects { first, second } if second == "I" => write!(
                f,
                "effect '{first}' is correlated with the overall mean: its contrast does not \
                 sum to zero, as happens when the runs are unbalanced (a run left out, or \
                 unequal replication), so its estimate absorbs part of the mean; balance the \
                 runs, or fit the terms together by regression"
            ),
            DoeError::PartiallyAliasedEffects { first, second } => write!(
                f,
                "effects '{first}' and '{second}' are partially aliased: their contrasts \
                 are correlated, so their sums of squares overlap; drop one, or fit the \
                 terms together by regression"
            ),
            DoeError::OverSpecifiedModel { terms, runs } => write!(
                f,
                "model has {terms} terms but the design has only {} degrees of freedom \
                 ({runs} runs); drop terms or add runs",
                runs.saturating_sub(1)
            ),
            DoeError::SingularModel => write!(
                f,
                "model matrix is singular: its columns are linearly dependent over these \
                 runs, so the coefficients cannot be identified; add runs or drop terms"
            ),
            DoeError::CoefficientCountMismatch {
                factors,
                expected,
                got,
            } => write!(
                f,
                "a quadratic model in {factors} factors has {expected} coefficients, got {got}"
            ),
            DoeError::TooFewContrasts { needed, got } => write!(
                f,
                "Lenth's method needs at least {needed} distinct contrasts, got {got}"
            ),
            DoeError::InvalidCodingRange { low, high } => {
                write!(f, "low ({low}) must be strictly less than high ({high})")
            }
            DoeError::EmptyResponses => write!(f, "responses must not be empty"),
            DoeError::TooFewReplicates { run, needed, got } => write!(
                f,
                "run {run}: needs at least {needed} replicates, got {got}"
            ),
            DoeError::NonPositiveResponse { run, value } => write!(
                f,
                "run {run}: responses must be greater than 0, got {value}"
            ),
            DoeError::InvalidDesirabilityLimits {
                index,
                goal,
                lower,
                target,
                upper,
            } => {
                if let Some(i) = index {
                    write!(f, "specs[{i}]: ")?;
                }
                let order = match goal {
                    crate::optimization::desirability::GoalType::Maximize => "lower < target",
                    crate::optimization::desirability::GoalType::Minimize => "target < upper",
                    crate::optimization::desirability::GoalType::Target => "lower < target < upper",
                };
                write!(
                    f,
                    "{goal:?} needs finite limits with {order}, got lower {lower}, target \
                     {target}, upper {upper}: the ramp has no width or runs backwards"
                )
            }
            DoeError::InvalidDesirabilityParameter {
                index,
                parameter,
                value,
            } => {
                if let Some(i) = index {
                    write!(f, "specs[{i}]: ")?;
                }
                let domain = if *parameter == "importance" {
                    "finite and not negative"
                } else {
                    "finite and greater than 0"
                };
                write!(f, "{parameter} must be {domain}, got {value}")
            }
            DoeError::NoWeightedResponse => write!(
                f,
                "every response has importance 0, so there is nothing to combine"
            ),
        }
    }
}

impl std::error::Error for DoeError {}

#[cfg(test)]
mod tests {
    use super::*;

    fn every_variant() -> Vec<DoeError> {
        vec![
            DoeError::InvalidFactorCount {
                min: 2,
                max: 7,
                got: 1,
            },
            DoeError::ParameterOutOfRange {
                parameter: "n_center",
                min: 1,
                max: None,
                got: 0,
            },
            DoeError::ParameterOutOfRange {
                parameter: "m",
                min: 1,
                max: Some(4),
                got: 0,
            },
            DoeError::UnsupportedFraction {
                k: 4,
                p: 2,
                supported: vec![(4, 1), (5, 1)],
            },
            DoeError::UnknownArray {
                array: "L5".into(),
                supported: &["L4", "L8"],
            },
            DoeError::EmptyDesign {
                runs: 0,
                factors: 0,
            },
            DoeError::DesignShapeMismatch {
                run: 1,
                expected: 3,
                got: 2,
            },
            DoeError::ResponseCountMismatch {
                expected: 8,
                got: 7,
            },
            DoeError::UnknownEffect { effect: "Z".into() },
            DoeError::EffectOutsideDesign {
                effect: "D".into(),
                factors: 3,
            },
            DoeError::NotTwoLevelCoded {
                run: 3,
                factor: 1,
                value: 0.5,
            },
            DoeError::AliasedEffects {
                first: "A".into(),
                second: "B:C".into(),
            },
            DoeError::AliasedEffects {
                first: "A:B:C".into(),
                second: "I".into(),
            },
            DoeError::PartiallyAliasedEffects {
                first: "A:B".into(),
                second: "C".into(),
            },
            DoeError::PartiallyAliasedEffects {
                first: "A".into(),
                second: "I".into(),
            },
            DoeError::OverSpecifiedModel { terms: 9, runs: 8 },
            DoeError::SingularModel,
            DoeError::CoefficientCountMismatch {
                factors: 2,
                expected: 6,
                got: 2,
            },
            DoeError::TooFewContrasts { needed: 3, got: 2 },
            DoeError::InvalidCodingRange {
                low: 2.0,
                high: 1.0,
            },
            DoeError::EmptyResponses,
            DoeError::TooFewReplicates {
                run: 0,
                needed: 2,
                got: 1,
            },
            DoeError::NonPositiveResponse {
                run: 2,
                value: -1.0,
            },
            DoeError::InvalidDesirabilityLimits {
                index: Some(1),
                goal: crate::optimization::desirability::GoalType::Target,
                lower: 10.0,
                target: 5.0,
                upper: 0.0,
            },
            DoeError::InvalidDesirabilityParameter {
                index: None,
                parameter: "importance",
                value: -1.0,
            },
            DoeError::NoWeightedResponse,
        ]
    }

    /// Line continuations without a trailing backslash put runs of spaces into
    /// these messages. A caller reads them verbatim.
    #[test]
    fn messages_carry_no_runs_of_spaces() {
        for e in every_variant() {
            let msg = e.to_string();
            assert!(!msg.contains("  "), "{msg}");
        }
    }

    #[test]
    fn unsupported_fraction_message_lists_the_table() {
        let e = DoeError::UnsupportedFraction {
            k: 4,
            p: 2,
            supported: vec![(4, 1), (5, 2)],
        };
        assert_eq!(
            e.to_string(),
            "unsupported design: 2^(4-2) not in standard table; supported (k,p): (4,1),(5,2)"
        );
    }

    /// `code()` and the serialized `code` are two spellings of one rule; a
    /// variant added to one and not the other would drift silently.
    #[cfg(feature = "wasm")]
    #[test]
    fn serialized_code_matches_code_method() {
        for e in every_variant() {
            let v = serde_json::to_value(&e).expect("DoeError serializes");
            assert_eq!(v["code"], e.code(), "{e:?}");
        }
    }

    /// The fields are copied onto a JavaScript `Error`, whose own `message`,
    /// `name`, `stack` and `cause` a field of the same name would overwrite.
    #[cfg(feature = "wasm")]
    #[test]
    fn serialized_fields_do_not_shadow_error_properties() {
        for e in every_variant() {
            let v = serde_json::to_value(&e).expect("DoeError serializes");
            let obj = v.as_object().expect("an object");
            for reserved in ["message", "name", "stack", "cause"] {
                assert!(!obj.contains_key(reserved), "{e:?} carries `{reserved}`");
            }
        }
    }
}
