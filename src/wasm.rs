//! WASM bindings for u-doe.
//!
//! Exposes DOE design generation and analysis functions to JavaScript/TypeScript
//! via `wasm-bindgen`. Only enabled when the `wasm` feature is active.
//!
//! # Feature
//!
//! Only compiled when the `wasm` feature is enabled:
//! ```toml
//! [dependencies]
//! u-doe = { version = "...", features = ["wasm"] }
//! ```

#![cfg(feature = "wasm")]

use serde::Serialize;
use wasm_bindgen::prelude::*;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn js_err(e: impl std::fmt::Display) -> JsValue {
    JsValue::from_str(&e.to_string())
}

fn to_js<T: Serialize>(value: &T) -> Result<JsValue, JsValue> {
    serde_wasm_bindgen::to_value(value).map_err(js_err)
}

// ---------------------------------------------------------------------------
// DTO types
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct DesignMatrixDto {
    data: Vec<Vec<f64>>,
    factor_names: Vec<String>,
    run_count: usize,
    factor_count: usize,
}

impl From<crate::design::DesignMatrix> for DesignMatrixDto {
    fn from(d: crate::design::DesignMatrix) -> Self {
        let run_count = d.run_count();
        let factor_count = d.factor_count();
        DesignMatrixDto {
            data: d.data,
            factor_names: d.factor_names,
            run_count,
            factor_count,
        }
    }
}

#[derive(Serialize)]
struct EffectRowDto {
    name: String,
    sum_of_squares: f64,
    df: usize,
    mean_square: f64,
    f_statistic: f64,
    p_value: Option<f64>,
}

#[derive(Serialize)]
struct DoeAnovaResultDto {
    effects: Vec<EffectRowDto>,
    residual_ss: f64,
    residual_df: usize,
    total_ss: f64,
    r_squared: f64,
    r_squared_adj: f64,
}

// ---------------------------------------------------------------------------
// Design generation
// ---------------------------------------------------------------------------

/// Generate a 2^k full factorial design.
///
/// Returns `{ data: [[f64]], factor_names: [str], run_count: usize, factor_count: usize }`.
///
/// # Errors
/// Returns an error string if `k` is out of range (1..=7).
#[wasm_bindgen]
pub fn full_factorial(k: usize) -> Result<JsValue, JsValue> {
    let design = crate::design::factorial::full_factorial(k).map_err(js_err)?;
    to_js(&DesignMatrixDto::from(design))
}

/// Generate a Central Composite Design (CCD).
///
/// `design_type`: `"FaceCentered"` | `"Rotatable"` | `"Inscribed"`
///
/// `n_center`: number of center point replicates (≥ 1, default 3 if 0 is passed is rejected).
///
/// Returns `{ data: [[f64]], factor_names: [str], run_count: usize, factor_count: usize }`.
///
/// # Errors
/// Returns an error string if `k` is out of range (2..=6), `n_center == 0`,
/// or `design_type` is unrecognised.
#[wasm_bindgen]
pub fn ccd(k: usize, design_type: &str, n_center: usize) -> Result<JsValue, JsValue> {
    let alpha_type = match design_type {
        "FaceCentered" => crate::design::ccd::AlphaType::FaceCentered,
        "Rotatable" => crate::design::ccd::AlphaType::Rotatable,
        "Inscribed" => crate::design::ccd::AlphaType::Inscribed,
        other => {
            return Err(js_err(format!(
                "unknown design_type '{}'; expected FaceCentered, Rotatable, or Inscribed",
                other
            )))
        }
    };
    let design = crate::design::ccd::ccd(k, alpha_type, n_center).map_err(js_err)?;
    to_js(&DesignMatrixDto::from(design))
}

/// Generate a Box-Behnken Design (BBD).
///
/// Supported: k = 3, 4, or 5.
///
/// Returns `{ data: [[f64]], factor_names: [str], run_count: usize, factor_count: usize }`.
///
/// # Errors
/// Returns an error string if `k` is not 3, 4, or 5, or `n_center == 0`.
#[wasm_bindgen]
pub fn box_behnken(k: usize, n_center: usize) -> Result<JsValue, JsValue> {
    let design = crate::design::box_behnken::box_behnken(k, n_center).map_err(js_err)?;
    to_js(&DesignMatrixDto::from(design))
}

/// Get a Taguchi orthogonal array.
///
/// `name`: `"L4"` | `"L8"` | `"L9"` | `"L12"` | `"L16"` | `"L18"` | `"L27"`
///
/// `k`: number of factors to use (must be ≤ max columns for the array).
///
/// Returns `{ data: [[f64]], factor_names: [str], run_count: usize, factor_count: usize }`.
///
/// # Errors
/// Returns an error string if the array name is unknown or `k` exceeds capacity.
#[wasm_bindgen]
pub fn taguchi_array(name: &str, k: usize) -> Result<JsValue, JsValue> {
    let design = crate::design::taguchi::taguchi_array(name, k).map_err(js_err)?;
    to_js(&DesignMatrixDto::from(design))
}

/// Generate a Definitive Screening Design (DSD).
///
/// Returns `{ data: [[f64]], factor_names: [str], run_count: usize, factor_count: usize }`.
///
/// # Errors
/// Returns an error string if `k` is out of the supported range.
#[wasm_bindgen]
pub fn definitive_screening(k: usize) -> Result<JsValue, JsValue> {
    let design = crate::design::definitive_screening(k).map_err(js_err)?;
    to_js(&DesignMatrixDto::from(design))
}

// ---------------------------------------------------------------------------
// Analysis
// ---------------------------------------------------------------------------

/// Perform DOE ANOVA.
///
/// `design_json`: JSON array-of-arrays `[[f64]]` — the coded design matrix (rows = runs).
/// `responses`: flat array of response values, one per run.
/// `factor_names_json`: JSON array of factor name strings used as column labels.
/// `effect_names_json`: JSON array of effect names to include (e.g. `["A","B","AB"]`).
///
/// Returns an ANOVA result object with `effects`, `residual_ss`, `residual_df`,
/// `total_ss`, `r_squared`, `r_squared_adj`.
///
/// # Errors
/// Returns an error string if dimensions do not match or JSON is malformed.
#[wasm_bindgen]
pub fn doe_anova(
    design_json: JsValue,
    responses: &[f64],
    factor_names_json: JsValue,
    effect_names_json: JsValue,
) -> Result<JsValue, JsValue> {
    let data: Vec<Vec<f64>> = serde_wasm_bindgen::from_value(design_json).map_err(js_err)?;
    let factor_names: Vec<String> =
        serde_wasm_bindgen::from_value(factor_names_json).map_err(js_err)?;
    let effect_names: Vec<String> =
        serde_wasm_bindgen::from_value(effect_names_json).map_err(js_err)?;

    let design = crate::design::DesignMatrix { data, factor_names };
    let effect_refs: Vec<&str> = effect_names.iter().map(|s| s.as_str()).collect();

    let result =
        crate::analysis::anova::doe_anova(&design, responses, &effect_refs).map_err(js_err)?;

    let dto = DoeAnovaResultDto {
        effects: result
            .effects
            .into_iter()
            .map(|e| EffectRowDto {
                name: e.name,
                sum_of_squares: e.sum_of_squares,
                df: e.df,
                mean_square: e.mean_square,
                f_statistic: e.f_statistic,
                p_value: e.p_value,
            })
            .collect(),
        residual_ss: result.residual_ss,
        residual_df: result.residual_df,
        total_ss: result.total_ss,
        r_squared: result.r_squared,
        r_squared_adj: result.r_squared_adj,
    };
    to_js(&dto)
}

/// Compute Taguchi Signal-to-Noise ratios.
///
/// `responses_json`: JSON array-of-arrays `[[f64]]` — one inner array per run,
///   containing the replicate measurements for that run.
/// `goal`: `"LargerIsBetter"` | `"SmallerIsBetter"` | `"NominalIsBest"`
///
/// Returns a flat `[f64]` of SN values in dB, one per run.
///
/// # Errors
/// Returns an error string if the goal string is unrecognised, or if the
/// response data violates the requirements for the chosen goal.
#[wasm_bindgen]
pub fn signal_to_noise(responses_json: JsValue, goal: &str) -> Result<JsValue, JsValue> {
    let responses: Vec<Vec<f64>> =
        serde_wasm_bindgen::from_value(responses_json).map_err(js_err)?;

    let sn_goal = match goal {
        "LargerIsBetter" => crate::analysis::taguchi_sn::SnGoal::LargerIsBetter,
        "SmallerIsBetter" => crate::analysis::taguchi_sn::SnGoal::SmallerIsBetter,
        "NominalIsBest" => crate::analysis::taguchi_sn::SnGoal::NominalIsBest,
        other => {
            return Err(js_err(format!(
                "unknown goal '{}'; expected LargerIsBetter, SmallerIsBetter, or NominalIsBest",
                other
            )))
        }
    };

    let sn_values =
        crate::analysis::taguchi_sn::signal_to_noise(&responses, sn_goal).map_err(js_err)?;
    to_js(&sn_values)
}

// ---------------------------------------------------------------------------
// Power analysis
// ---------------------------------------------------------------------------

/// Compute the statistical power of a 2^(k-p) factorial design.
///
/// Uses the normal approximation. Returns power in [0, 1], or 0.0 for invalid
/// inputs (k=0, p≥k, n_replicates=0, effect_size≤0, sigma≤0).
///
/// # Arguments
/// * `k` — total number of factors
/// * `p` — number of generators (p=0 = full factorial)
/// * `n_replicates` — number of replicates
/// * `effect_size` — detectable effect δ (in response units)
/// * `sigma` — process standard deviation σ
/// * `alpha` — type-I error rate (e.g. 0.05)
#[wasm_bindgen]
pub fn two_level_factorial_power(
    k: usize,
    p: usize,
    n_replicates: usize,
    effect_size: f64,
    sigma: f64,
    alpha: f64,
) -> f64 {
    crate::power::two_level_factorial_power(k, p, n_replicates, effect_size, sigma, alpha)
}
