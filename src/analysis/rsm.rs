//! Response Surface Methodology (RSM) — second-order model fitting.
//!
//! Fits a full quadratic model to experimental data using ordinary least
//! squares. Provides prediction and steepest ascent path computation.
//!
//! # Model
//!
//! `ŷ = β₀ + Σ βᵢxᵢ + Σ βᵢᵢxᵢ² + Σᵢ<ⱼ βᵢⱼxᵢxⱼ`
//!
//! Coefficients are estimated as β = (X'X)⁻¹ X'y.
//!
//! Reference: Montgomery, D.C. (2019). *Introduction to Statistical Quality
//! Control*, 8th ed., Sections 11.2–11.3. Wiley.

use u_numflow::matrix::Matrix;

use crate::{design::DesignMatrix, error::DoeError};

/// Fitted second-order response surface model.
#[derive(Debug, Clone)]
pub struct RsmModel {
    /// Coefficients in order: [β₀, β₁..βₖ, β₁₁..βₖₖ, β₁₂..β(k−1)k].
    pub coefficients: Vec<f64>,
    /// Coefficient of determination R².
    pub r_squared: f64,
    /// Number of factors.
    pub factor_count: usize,
}

impl RsmModel {
    /// Predict the response at coded factor levels `x`.
    ///
    /// `x.len()` must equal `factor_count`.
    pub fn predict(&self, x: &[f64]) -> f64 {
        let row = model_row(x);
        row.iter()
            .zip(self.coefficients.iter())
            .map(|(xi, bi)| xi * bi)
            .sum()
    }
}

/// One step along the steepest ascent path.
#[derive(Debug, Clone)]
pub struct AscentStep {
    /// Coded factor values at this step.
    pub coded: Vec<f64>,
    /// Step number (1-indexed).
    pub step_number: usize,
}

/// Fit a second-order RSM model using OLS.
///
/// # Arguments
/// * `design`    — Experimental design (must have enough runs to estimate all coefficients)
/// * `responses` — Observed responses
///
/// # Errors
///
/// Returns `Err` if `responses.len() != run_count`, or if the model matrix
/// is singular (not enough unique support points).
///
/// # Examples
///
/// ```
/// use u_doe::design::ccd::{ccd, AlphaType};
/// use u_doe::analysis::rsm::fit_rsm;
///
/// let design = ccd(2, AlphaType::FaceCentered, 3).unwrap();
/// let n = design.run_count();
/// let responses: Vec<f64> = (0..n).map(|i| i as f64).collect();
/// let model = fit_rsm(&design, &responses).unwrap();
/// assert_eq!(model.factor_count, 2);
/// ```
pub fn fit_rsm(design: &DesignMatrix, responses: &[f64]) -> Result<RsmModel, DoeError> {
    let n = design.run_count();
    let k = design.factor_count();

    if responses.len() != n {
        return Err(DoeError::InsufficientResponses {
            expected: n,
            got: responses.len(),
        });
    }

    // Build model matrix X (n × p)
    let x_rows: Vec<Vec<f64>> = (0..n)
        .map(|run| {
            let factors: Vec<f64> = (0..k).map(|j| design.get(run, j)).collect();
            model_row(&factors)
        })
        .collect();

    let x_refs: Vec<&[f64]> = x_rows.iter().map(|r| r.as_slice()).collect();
    let x_mat = Matrix::from_rows(&x_refs);
    let xt = x_mat.transpose();

    // X'X
    let xtx = xt
        .mul_mat(&x_mat)
        .map_err(|e| DoeError::MatrixError(e.to_string()))?;

    // X'y  (xt is p×n, responses is n-vector)
    let xty = xt
        .mul_vec(responses)
        .map_err(|e| DoeError::MatrixError(e.to_string()))?;

    // β = (X'X)⁻¹ X'y
    let xtx_inv = xtx
        .inverse()
        .map_err(|e| DoeError::MatrixError(e.to_string()))?;

    let coefficients = xtx_inv
        .mul_vec(&xty)
        .map_err(|e| DoeError::MatrixError(e.to_string()))?;

    // R² = 1 - SS_res / SS_tot
    let y_mean = responses.iter().sum::<f64>() / n as f64;
    let ss_tot: f64 = responses.iter().map(|&y| (y - y_mean).powi(2)).sum();
    let ss_res: f64 = (0..n)
        .map(|i| {
            let factors: Vec<f64> = (0..k).map(|j| design.get(i, j)).collect();
            let row = model_row(&factors);
            let pred: f64 = row
                .iter()
                .zip(coefficients.iter())
                .map(|(xi, bi)| xi * bi)
                .sum();
            (responses[i] - pred).powi(2)
        })
        .sum();

    let r_squared = if ss_tot > 1e-12 {
        1.0 - ss_res / ss_tot
    } else {
        1.0
    };

    Ok(RsmModel {
        coefficients,
        r_squared,
        factor_count: k,
    })
}

/// Compute steepest ascent path from a fitted model.
///
/// Uses the first-order (linear) terms of the model to determine the
/// ascent direction. Each step moves `step_size` coded units in the
/// direction of the gradient normalised by the largest absolute linear
/// coefficient.
///
/// Returns an empty `Vec` when `n_steps == 0` or all linear coefficients
/// are essentially zero.
pub fn steepest_ascent(model: &RsmModel, n_steps: usize, step_size: f64) -> Vec<AscentStep> {
    let k = model.factor_count;
    // Linear coefficients: indices 1..=k in coefficients vector
    let linear: Vec<f64> = model.coefficients[1..=k].to_vec();

    // Normalise by largest absolute coefficient
    let max_abs = linear.iter().map(|b| b.abs()).fold(0.0_f64, f64::max);
    if max_abs < 1e-12 || n_steps == 0 {
        return Vec::new();
    }
    let direction: Vec<f64> = linear.iter().map(|b| b / max_abs).collect();

    (1..=n_steps)
        .map(|step| {
            let t = step as f64 * step_size;
            AscentStep {
                coded: direction.iter().map(|d| d * t).collect(),
                step_number: step,
            }
        })
        .collect()
}

/// Build a model matrix row for factor values `x`.
///
/// Order: [1, x₁..xₖ, x₁²..xₖ², x₁x₂..x(k-1)xₖ]
fn model_row(x: &[f64]) -> Vec<f64> {
    let k = x.len();
    let mut row = Vec::with_capacity(1 + k + k + k * (k - 1) / 2);
    row.push(1.0); // intercept
    row.extend_from_slice(x); // linear
    for &xi in x {
        row.push(xi * xi); // pure quadratic
    }
    for i in 0..k {
        for j in (i + 1)..k {
            row.push(x[i] * x[j]); // interaction
        }
    }
    row
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::design::ccd::{ccd, AlphaType};

    #[test]
    fn rsm_fit_k2_coefficient_count() {
        // k=2: 1(intercept) + 2(linear) + 2(quadratic) + 1(interaction) = 6
        let design = ccd(2, AlphaType::FaceCentered, 5).unwrap();
        let n = design.run_count();
        let responses: Vec<f64> = (0..n).map(|i| (i as f64) * 2.0 + 1.0).collect();
        let model = fit_rsm(&design, &responses).unwrap();
        assert_eq!(model.coefficients.len(), 6);
        assert_eq!(model.factor_count, 2);
    }

    #[test]
    fn rsm_fit_k3_coefficient_count() {
        // k=3: 1 + 3 + 3 + 3 = 10
        let design = ccd(3, AlphaType::FaceCentered, 3).unwrap();
        let n = design.run_count();
        let responses: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let model = fit_rsm(&design, &responses).unwrap();
        assert_eq!(model.coefficients.len(), 10);
    }

    #[test]
    fn rsm_predict_at_center() {
        // ccd(2, FaceCentered, 5): 4 factorial + 4 axial + 5 center = 13 runs
        let design = ccd(2, AlphaType::FaceCentered, 5).unwrap();
        assert_eq!(design.run_count(), 13);
        let responses = vec![
            39.3, 40.0, 40.9, 41.5, // factorial (4)
            40.3, 40.5, 41.3, 40.0, // axial (4)
            40.8, 40.9, 41.0, 40.5, 40.7, // center points (5)
        ];
        let model = fit_rsm(&design, &responses).unwrap();
        let pred_center = model.predict(&[0.0, 0.0]);
        // Should be near mean of center points (~40.78)
        assert!(
            (pred_center - 40.78).abs() < 1.0,
            "pred={pred_center}"
        );
    }

    #[test]
    fn rsm_r_squared_in_range() {
        let design = ccd(2, AlphaType::FaceCentered, 3).unwrap();
        let n = design.run_count();
        let responses: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let model = fit_rsm(&design, &responses).unwrap();
        assert!(
            model.r_squared >= 0.0 && model.r_squared <= 1.0 + 1e-9
        );
    }

    #[test]
    fn rsm_predict_consistency() {
        let design = ccd(2, AlphaType::FaceCentered, 1).unwrap();
        let n = design.run_count();
        let responses: Vec<f64> = (0..n).map(|i| (i as f64).powi(2)).collect();
        let model = fit_rsm(&design, &responses).unwrap();
        // Prediction should be finite for any coded input
        assert!(model.predict(&[0.5, -0.5]).is_finite());
        assert!(model.predict(&[-1.0, 1.0]).is_finite());
    }

    #[test]
    fn steepest_ascent_returns_steps() {
        let design = ccd(2, AlphaType::FaceCentered, 3).unwrap();
        let n = design.run_count();
        let responses: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let model = fit_rsm(&design, &responses).unwrap();
        let steps = steepest_ascent(&model, 5, 0.5);
        assert_eq!(steps.len(), 5);
        assert_eq!(steps[0].step_number, 1);
        assert_eq!(steps[0].coded.len(), 2);
    }

    #[test]
    fn rsm_insufficient_data() {
        let design = ccd(2, AlphaType::FaceCentered, 1).unwrap();
        // Fewer responses than runs
        let responses = vec![1.0, 2.0];
        assert!(fit_rsm(&design, &responses).is_err());
    }
}
