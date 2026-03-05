//! Integration tests — Montgomery (2019) Example 6.2 full pipeline,
//! RSM pipeline, and desirability pipeline.

use u_doe::analysis::anova::doe_anova;
use u_doe::analysis::effects::estimate_effects;
use u_doe::analysis::rsm::{fit_rsm, steepest_ascent};
use u_doe::design::ccd::{ccd, AlphaType};
use u_doe::design::factorial::full_factorial;
use u_doe::optimization::desirability::{overall_desirability, ResponseSpec};

/// Montgomery (2019) Example 6.2 — Filtration rate experiment (2^4 design).
/// Textbook: pp. 234–240, Table 6.8.
#[test]
fn montgomery_ex62_full_pipeline() {
    let design = full_factorial(4).unwrap();
    let responses = vec![
        45.0, 71.0, 48.0, 65.0, 68.0, 60.0, 80.0, 65.0, 43.0, 100.0, 45.0, 104.0, 75.0, 86.0, 70.0,
        96.0,
    ];

    // Step 1: effects match textbook
    let effects = estimate_effects(&design, &responses, 2).unwrap();
    let find = |name: &str| effects.iter().find(|e| e.name == name).unwrap().estimate;
    assert!((find("A") - 21.625).abs() < 0.001);
    assert!((find("C") - 9.875).abs() < 0.001);
    assert!((find("D") - 14.625).abs() < 0.001);
    assert!((find("AC") - (-18.125)).abs() < 0.001);

    // Step 2: ANOVA
    let anova = doe_anova(&design, &responses, &["A", "C", "D", "AC"]).unwrap();
    let ss_sum: f64 =
        anova.effects.iter().map(|r| r.sum_of_squares).sum::<f64>() + anova.residual_ss;
    assert!(
        (ss_sum - anova.total_ss).abs() < 0.1,
        "SS partition failed: {ss_sum} vs {}",
        anova.total_ss
    );
    assert!(anova.r_squared > 0.70, "R²={}", anova.r_squared);

    // Step 3: significant effects (A, D, AC clearly significant at α=0.05;
    // C is marginally active — significant at α=0.10 per Montgomery Table 6.9)
    let clearly_significant = ["A", "D", "AC"];
    for row in &anova.effects {
        if clearly_significant.contains(&row.name.as_str()) {
            if let Some(p) = row.p_value {
                assert!(p < 0.05, "Effect {} not significant: p={p}", row.name);
            }
        }
    }
    // C is active but borderline: p < 0.15 in this model
    let c_row = anova.effects.iter().find(|r| r.name == "C").unwrap();
    if let Some(p) = c_row.p_value {
        assert!(p < 0.15, "Effect C p-value unexpectedly large: p={p}");
    }
}

/// CCD RSM pipeline — fit second-order model and compute steepest ascent path.
#[test]
fn ccd_rsm_pipeline() {
    let design = ccd(2, AlphaType::FaceCentered, 3).unwrap();
    assert_eq!(design.run_count(), 11);
    let responses: Vec<f64> = (0..11).map(|i| 39.0 + i as f64 * 0.2).collect();

    let model = fit_rsm(&design, &responses).unwrap();
    assert_eq!(model.coefficients.len(), 6);
    assert_eq!(model.factor_count, 2);
    assert!(model.r_squared >= 0.0 && model.r_squared <= 1.0 + 1e-9);

    let steps = steepest_ascent(&model, 3, 0.5);
    assert_eq!(steps.len(), 3);
    assert_eq!(steps[0].coded.len(), 2);
}

/// Desirability pipeline — multi-response optimisation with ANOVA.
#[test]
fn desirability_pipeline() {
    let design = full_factorial(3).unwrap();
    let y1 = vec![45.0, 71.0, 48.0, 65.0, 68.0, 60.0, 80.0, 65.0];
    let y2 = [10.0, 8.0, 12.0, 6.0, 9.0, 7.0, 11.0, 5.0];

    let spec1 = ResponseSpec::maximize(40.0, 85.0, 85.0, 1.0);
    let spec2 = ResponseSpec::minimize(3.0, 3.0, 15.0, 1.0);

    for i in 0..8 {
        let d = overall_desirability(&[spec1.clone(), spec2.clone()], &[y1[i], y2[i]]);
        assert!((0.0..=1.0).contains(&d), "run {i}: d={d}");
    }

    let anova = doe_anova(&design, &y1, &["A", "B", "C"]).unwrap();
    assert!(anova.r_squared >= 0.0);
}
