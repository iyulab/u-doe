//! Design generators against the properties that define them.
//!
//! No library generates these designs, so the reference is the definition,
//! restated here without reading the generator: a full factorial is every
//! combination of the levels; a central composite design is a cube, a star
//! and a centre, and it is rotatable when Box and Hunter's fourth-moment
//! condition holds; a Box-Behnken design is every pair of factors at the four
//! sign combinations with the rest at the middle; a definitive screening design
//! is a fold-over whose main effects are orthogonal to every quadratic and
//! two-factor-interaction column (Jones & Nachtsheim 2011).

use std::collections::BTreeSet;
use u_doe::design::box_behnken::box_behnken_design;
use u_doe::design::ccd::{ccd_design, AlphaType};
use u_doe::design::definitive_screening;
use u_doe::design::factorial::full_factorial;
use u_doe::design::PointType;

const TOL: f64 = 1e-9;

fn column(rows: &[Vec<f64>], j: usize) -> Vec<f64> {
    rows.iter().map(|r| r[j]).collect()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn n_choose_2(k: usize) -> usize {
    k * (k - 1) / 2
}

/// Every combination of two levels appears exactly once, and every effect --
/// each non-empty product of columns -- sums to zero.
#[test]
fn full_factorial_is_every_combination_and_every_effect_is_balanced() {
    for k in 2..=7usize {
        let d = full_factorial(k).unwrap();
        assert_eq!(d.data.len(), 1 << k, "k={k}");
        let codes: BTreeSet<usize> = d
            .data
            .iter()
            .map(|r| {
                r.iter().enumerate().fold(0usize, |acc, (j, &v)| {
                    assert!(v == 1.0 || v == -1.0, "k={k}: level {v}");
                    acc | (usize::from(v > 0.0) << j)
                })
            })
            .collect();
        assert_eq!(codes.len(), 1 << k, "k={k}: a combination repeats");
        if k <= 5 {
            for mask in 1usize..(1 << k) {
                let sum: f64 = d
                    .data
                    .iter()
                    .map(|r| {
                        (0..k)
                            .filter(|j| mask >> j & 1 == 1)
                            .map(|j| r[j])
                            .product::<f64>()
                    })
                    .sum();
                assert!(sum.abs() < TOL, "k={k} effect {mask:b} sums to {sum}");
            }
        }
    }
}

/// Counts, block coordinates, and (for rotatable/inscribed) Box-Hunter
/// rotatability: odd moments vanish and `[iiii] = 3 [iijj]`.
#[test]
fn ccd_blocks_have_the_defined_coordinates_and_rotatable_alpha_satisfies_box_hunter() {
    for k in 2..=6usize {
        let f = (1usize << k) as f64;
        for n_center in [1usize, 3] {
            for alpha in [
                AlphaType::FaceCentered,
                AlphaType::Rotatable,
                AlphaType::Inscribed,
            ] {
                let d = ccd_design(k, alpha, n_center).unwrap();
                let rows = &d.matrix.data;
                assert_eq!(rows.len(), (1 << k) + 2 * k + n_center, "k={k}");
                let a_rot = f.powf(0.25);
                let (cube, star) = match alpha {
                    AlphaType::FaceCentered => (1.0, 1.0),
                    AlphaType::Rotatable => (1.0, a_rot),
                    AlphaType::Inscribed => (1.0 / a_rot, 1.0),
                };
                for (row, t) in rows.iter().zip(&d.point_types) {
                    match t {
                        PointType::Factorial => {
                            assert!(row.iter().all(|v| (v.abs() - cube).abs() < TOL))
                        }
                        PointType::Axial => {
                            let nz: Vec<f64> =
                                row.iter().copied().filter(|v| v.abs() > TOL).collect();
                            assert_eq!(nz.len(), 1);
                            assert!((nz[0].abs() - star).abs() < TOL);
                        }
                        PointType::Center => assert!(row.iter().all(|&v| v == 0.0)),
                        PointType::Edge => panic!("a CCD has no edge points"),
                    }
                }
                // each axis has one star point at each side
                for j in 0..k {
                    let hi = rows
                        .iter()
                        .zip(&d.point_types)
                        .any(|(r, t)| *t == PointType::Axial && (r[j] - star).abs() < TOL);
                    let lo = rows
                        .iter()
                        .zip(&d.point_types)
                        .any(|(r, t)| *t == PointType::Axial && (r[j] + star).abs() < TOL);
                    assert!(hi && lo, "k={k} axis {j}");
                }
                // moments
                let m = |f: &dyn Fn(&[f64]) -> f64| rows.iter().map(|r| f(r)).sum::<f64>();
                for i in 0..k {
                    assert!(m(&|r| r[i]).abs() < TOL);
                    assert!(m(&|r| r[i].powi(3)).abs() < TOL);
                    for j in 0..k {
                        if i != j {
                            assert!(m(&|r| r[i] * r[j]).abs() < TOL);
                        }
                    }
                }
                let rotatable = !matches!(alpha, AlphaType::FaceCentered);
                let iiii = m(&|r| r[0].powi(4));
                let iijj = m(&|r| r[0].powi(2) * r[1].powi(2));
                if rotatable {
                    assert!(
                        (iiii - 3.0 * iijj).abs() < 1e-8 * iiii.max(1.0),
                        "k={k} {alpha:?}: [iiii]={iiii} 3[iijj]={}",
                        3.0 * iijj
                    );
                } else {
                    assert!((iiii - 3.0 * iijj).abs() > 1e-3, "k={k} face-centred");
                }
            }
        }
    }
}

/// Every pair of factors takes the four sign combinations once with the others
/// at the middle; nothing else but the centre rows.
#[test]
fn box_behnken_is_every_pair_at_four_sign_combinations() {
    for k in 3..=5usize {
        for n_center in [1usize, 3, 6] {
            let d = box_behnken_design(k, n_center).unwrap();
            let rows = &d.matrix.data;
            assert_eq!(rows.len(), 4 * n_choose_2(k) + n_center, "k={k}");
            let mut seen: BTreeSet<(usize, usize, i8, i8)> = BTreeSet::new();
            for (row, t) in rows.iter().zip(&d.point_types) {
                let nz: Vec<usize> = (0..k).filter(|&j| row[j] != 0.0).collect();
                match t {
                    PointType::Center => assert!(nz.is_empty()),
                    PointType::Edge => {
                        assert_eq!(nz.len(), 2);
                        assert!(nz.iter().all(|&j| row[j].abs() == 1.0));
                        let key = (nz[0], nz[1], row[nz[0]] as i8, row[nz[1]] as i8);
                        assert!(seen.insert(key), "k={k} repeats {key:?}");
                    }
                    other => panic!("a Box-Behnken design has no {other:?} point"),
                }
            }
            assert_eq!(seen.len(), 4 * n_choose_2(k), "k={k}: a pair is missing");
            for i in 0..k {
                for j in (i + 1)..k {
                    assert!(
                        dot(&column(rows, i), &column(rows, j)).abs() < TOL,
                        "k={k} columns {i},{j}"
                    );
                }
            }
        }
    }
}

/// Fold-over, a centre run, and main effects orthogonal to every quadratic
/// and every two-factor interaction and to each other. Even k takes 2k+1 runs;
/// odd k takes 2k+3 (the conference matrix of order k+1 less one column), so a
/// pair of runs has no zero.
#[test]
fn definitive_screening_is_a_fold_over_with_main_effects_clear_of_second_order_terms() {
    for k in 2..=12usize {
        let d = definitive_screening(k).unwrap();
        let rows = &d.data;
        let n = rows.len();
        let m = k + k % 2;
        assert_eq!(n, 2 * m + 1, "k={k}");
        assert!(rows
            .iter()
            .all(|r| r.iter().all(|&v| v == -1.0 || v == 0.0 || v == 1.0)));
        assert!(
            rows[n - 1].iter().all(|&v| v == 0.0),
            "k={k}: last row is the centre"
        );
        for i in 0..m {
            for (&v, &w) in rows[i].iter().zip(&rows[m + i]) {
                assert_eq!(v, -w, "k={k} fold-over row {i}");
            }
            let zeros = rows[i].iter().filter(|&&v| v == 0.0).count();
            // a conference-matrix row has one zero (its diagonal), unless that
            // is the dropped column
            assert_eq!(
                zeros,
                usize::from(i < k),
                "k={k} row {i}: factors at the middle"
            );
        }
        for a in 0..k {
            let main = column(rows, a);
            for b in 0..k {
                let sq: Vec<f64> = column(rows, b).iter().map(|v| v * v).collect();
                assert!(
                    dot(&main, &sq).abs() < TOL,
                    "k={k} main {a} vs quadratic {b}"
                );
                for c in (b + 1)..k {
                    let inter: Vec<f64> = rows.iter().map(|r| r[b] * r[c]).collect();
                    assert!(
                        dot(&main, &inter).abs() < TOL,
                        "k={k} main {a} vs interaction {b}{c}"
                    );
                }
            }
            for b in (a + 1)..k {
                assert!(
                    dot(&main, &column(rows, b)).abs() < TOL,
                    "k={k} mains {a},{b} not orthogonal"
                );
            }
        }
    }
}
