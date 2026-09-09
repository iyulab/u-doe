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
            DoeError::UnknownEffect { name } => {
                write!(
                    f,
                    "unknown effect '{name}': effect names are factor names joined with ':' \
                     (e.g. \"A\", \"A:B\")"
                )
            }
            DoeError::MatrixError(msg) => write!(f, "matrix error: {msg}"),
            DoeError::UnsupportedDesign(msg) => write!(f, "unsupported design: {msg}"),
            DoeError::NotTwoLevelCoded { run, factor, value } => write!(
                f,
                "design is not two-level coded: run {run}, factor {factor} is {value},                  expected -1 or +1. Effect estimation multiplies factor columns to                  form contrasts, which is only defined for a two-level design; a                  design carrying centre points or axial points (central composite,                  Box-Behnken, definitive screening) must be fitted with                  `analysis::rsm::fit_rsm` instead"
            ),
            DoeError::AliasedEffects { first, second } => write!(
                f,
                "effects '{first}' and '{second}' share one contrast column and are                  therefore aliased in this design: they cannot both enter the model,                  because their sum of squares would be counted twice"
            ),
            DoeError::OverSpecifiedModel { terms, runs } => write!(
                f,
                "model has {terms} terms but the design has only {} degrees of freedom                  ({runs} runs); drop terms or add runs",
                runs.saturating_sub(1)
            ),
        }
    }
}

impl std::error::Error for DoeError {}
