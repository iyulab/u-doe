pub mod anova;
pub mod effects;
pub mod least_squares;
pub mod lenth;
pub mod rsm;
pub mod taguchi_sn;

/// Whether `responses` have no variation to explain: their total sum of
/// squares is below what rounding the mean leaves, judged against the
/// responses' own size. A fixed number would call a response measured in small
/// units constant and report an exact fit, or no effect, whatever the data say.
pub(crate) fn has_no_variation(total_ss: f64, responses: &[f64]) -> bool {
    let size: f64 = responses.iter().map(|y| y * y).sum();
    total_ss <= 1e-20 * size
}

#[cfg(test)]
mod unit_invariance {
    use crate::analysis::anova::doe_anova;
    use crate::analysis::effects::estimate_effects;
    use crate::analysis::least_squares::fit_least_squares;
    use crate::analysis::rsm::fit_rsm;
    use crate::design::DesignMatrix;

    /// A response in other units has the same fit: R^2 and the share of each
    /// effect are ratios, so the unit cannot change them.
    #[test]
    fn r_squared_and_effect_shares_do_not_depend_on_the_unit() {
        let b = 2.0_f64.sqrt();
        let two_level = DesignMatrix {
            data: vec![
                vec![-1.0, -1.0],
                vec![1.0, -1.0],
                vec![-1.0, 1.0],
                vec![1.0, 1.0],
                vec![-1.0, -1.0],
                vec![1.0, -1.0],
                vec![-1.0, 1.0],
                vec![1.0, 1.0],
            ],
            factor_names: vec!["A".into(), "B".into()],
        };
        let central = DesignMatrix {
            data: vec![
                vec![-1.0, -1.0],
                vec![1.0, -1.0],
                vec![-1.0, 1.0],
                vec![1.0, 1.0],
                vec![b, 0.0],
                vec![-b, 0.0],
                vec![0.0, b],
                vec![0.0, -b],
                vec![0.0, 0.0],
                vec![0.0, 0.0],
            ],
            factor_names: vec!["A".into(), "B".into()],
        };
        let y2 = [10.1, 12.3, 9.8, 15.2, 10.4, 12.0, 10.0, 15.6];
        let yc = [5.1, 7.9, 4.2, 9.8, 7.0, 3.1, 6.4, 5.5, 6.1, 5.9];
        let scaled = |y: &[f64]| y.iter().map(|v| v * 1e-9).collect::<Vec<f64>>();
        let names = ["A", "B", "A:B"];

        let r2 = |y: &[f64]| doe_anova(&two_level, y, &names).expect("valid").r_squared;
        let ls = |y: &[f64]| fit_least_squares(&two_level, y, &names).expect("valid");
        let rsm = |y: &[f64]| fit_rsm(&central, y).expect("valid").r_squared;
        let share = |y: &[f64]| {
            estimate_effects(&two_level, y, 2)
                .expect("valid")
                .iter()
                .map(|e| e.percent_contribution)
                .collect::<Vec<f64>>()
        };

        assert!(
            r2(&y2) < 1.0 && r2(&y2) > 0.5,
            "the fixture has an imperfect fit"
        );
        assert!((r2(&y2) - r2(&scaled(&y2))).abs() < 1e-9);
        assert!((ls(&y2).r_squared - ls(&scaled(&y2)).r_squared).abs() < 1e-9);
        assert!((ls(&y2).r_squared_adj - ls(&scaled(&y2)).r_squared_adj).abs() < 1e-9);
        assert!((rsm(&yc) - rsm(&scaled(&yc))).abs() < 1e-9);
        for (a, b) in share(&y2).iter().zip(share(&scaled(&y2))) {
            assert!((a - b).abs() < 1e-9, "effect share {a} vs {b}");
        }
    }
}
