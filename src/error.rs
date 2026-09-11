/// Errors produced by u-doe operations.
#[derive(Debug, Clone, PartialEq)]
pub enum DoeError {
    /// Invalid number of factors for this design type.
    InvalidFactorCount { min: usize, max: usize, got: usize },
    /// Insufficient response data.
    InsufficientResponses { expected: usize, got: usize },
    /// Specification limits are invalid.
    InvalidSpecification(String),
    /// Effect name not found among the effects estimable for the design.
    UnknownEffect { name: String },
    /// Matrix operation failed.
    MatrixError(String),
    /// Requested design not supported in v0.1.
    UnsupportedDesign(String),
    /// A design matrix that must be two-level coded (every entry -1 or +1)
    /// carries some other value.
    NotTwoLevelCoded {
        run: usize,
        factor: usize,
        value: f64,
    },
    /// Two requested effects share the same contrast column, so they are
    /// aliased in this design and cannot be estimated separately.
    AliasedEffects { first: String, second: String },
    /// Two requested effects -- or an effect and the overall mean, written `I` --
    /// have contrasts that are correlated without being identical. The per-term
    /// sums of squares then overlap, so they do not add up to the model sum of
    /// squares and no ANOVA can be built from them. Plackett-Burman designs
    /// produce this for two-factor interactions.
    PartiallyAliasedEffects { first: String, second: String },
    /// The requested model has more terms than the design has degrees of
    /// freedom to spend on them.
    OverSpecifiedModel { terms: usize, runs: usize },
}

impl std::fmt::Display for DoeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DoeError::InvalidFactorCount { min, max, got } => {
                write!(f, "invalid factor count: expected {min}..={max}, got {got}")
            }
            DoeError::InsufficientResponses { expected, got } => {
                write!(f, "insufficient responses: expected {expected}, got {got}")
            }
            DoeError::InvalidSpecification(msg) => write!(f, "invalid specification: {msg}"),
            DoeError::UnknownEffect { name } => write!(
                f,
                "unknown effect '{name}': effect names are factor names joined with ':' \
                 (e.g. \"A\", \"A:B\")"
            ),
            DoeError::MatrixError(msg) => write!(f, "matrix error: {msg}"),
            DoeError::UnsupportedDesign(msg) => write!(f, "unsupported design: {msg}"),
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
            DoeError::AliasedEffects { first, second } => write!(
                f,
                "effects '{first}' and '{second}' share one contrast column and are \
                 therefore aliased in this design: they cannot both enter the model, \
                 because their sum of squares would be counted twice"
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
        }
    }
}

impl std::error::Error for DoeError {}

#[cfg(test)]
mod tests {
    use super::*;

    /// Line continuations without a trailing backslash put runs of spaces into
    /// these messages. A caller reads them verbatim.
    #[test]
    fn messages_carry_no_runs_of_spaces() {
        let errors = [
            DoeError::NotTwoLevelCoded {
                run: 3,
                factor: 1,
                value: 0.5,
            },
            DoeError::AliasedEffects {
                first: "A".into(),
                second: "B:C".into(),
            },
            DoeError::PartiallyAliasedEffects {
                first: "A:B".into(),
                second: "C".into(),
            },
            DoeError::OverSpecifiedModel { terms: 9, runs: 8 },
            DoeError::UnknownEffect { name: "Z".into() },
        ];
        for e in errors {
            let msg = e.to_string();
            assert!(!msg.contains("  "), "{msg}");
        }
    }
}
