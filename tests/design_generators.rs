//! Fractional factorials, Plackett-Burman designs, Taguchi arrays and mixture
//! designs against the properties that define them (no library generates
//! these either; the reference is the definition, restated without reading
//! the generator).

use std::collections::BTreeSet;
use u_doe::design::factorial::{
    fractional_factorial, fractional_factorial_info, standard_fractions,
};
use u_doe::design::mixture::{simplex_centroid, simplex_lattice};
use u_doe::design::plackett_burman::plackett_burman;
use u_doe::design::taguchi::{max_factors, taguchi_array};

const TOL: f64 = 1e-9;

fn column(rows: &[Vec<f64>], j: usize) -> Vec<f64> {
    rows.iter().map(|r| r[j]).collect()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn word_mask(word: &str) -> usize {
    word.bytes().fold(0, |m, b| m | 1 << (b - b'A'))
}

/// The defining contrast subgroup read off the matrix: every set of columns
/// whose product is the all-ones column. The reported resolution, defining
/// relation and generators must all say what the matrix does.
#[test]
fn fractional_factorial_info_describes_the_matrix_it_comes_with() {
    let table = standard_fractions();
    assert!(table
        .windows(2)
        .all(|w| (w[0].k, w[0].p) < (w[1].k, w[1].p)));
    for f in table {
        let (k, p) = (f.k, f.p);
        let d = fractional_factorial(k, p).unwrap();
        let rows = &d.data;
        assert_eq!(rows.len(), 1 << (k - p), "{k}-{p}");
        assert!(rows
            .iter()
            .all(|r| r.len() == k && r.iter().all(|&v| v == 1.0 || v == -1.0)));
        let distinct: BTreeSet<Vec<i8>> = rows
            .iter()
            .map(|r| r.iter().map(|&v| v as i8).collect())
            .collect();
        assert_eq!(distinct.len(), rows.len(), "{k}-{p}: a run repeats");

        let words: BTreeSet<usize> = (1usize..(1 << k))
            .filter(|&mask| {
                rows.iter().all(|r| {
                    (0..k)
                        .filter(|j| mask >> j & 1 == 1)
                        .map(|j| r[j])
                        .product::<f64>()
                        == 1.0
                })
            })
            .collect();
        assert_eq!(
            words.len(),
            (1 << p) - 1,
            "{k}-{p}: defining contrast subgroup size"
        );
        let shortest = words.iter().map(|m| m.count_ones() as usize).min().unwrap();
        assert_eq!(f.resolution.order(), shortest, "{k}-{p}: resolution");

        let printed: BTreeSet<usize> = f
            .defining_relation
            .split('=')
            .skip(1)
            .map(word_mask)
            .collect();
        assert_eq!(printed, words, "{k}-{p}: {}", f.defining_relation);
        let info = fractional_factorial_info(k, p).unwrap();
        assert_eq!(info.defining_relation, f.defining_relation);

        assert_eq!(f.generators.len(), p);
        for (i, g) in f.generators.iter().enumerate() {
            let (lhs, rhs) = g.split_once('=').unwrap();
            let target = (lhs.as_bytes()[0] - b'A') as usize;
            assert_eq!(target, k - p + i, "{k}-{p}: generator order");
            let cols: Vec<usize> = rhs.bytes().map(|b| (b - b'A') as usize).collect();
            for r in rows {
                assert_eq!(
                    r[target],
                    cols.iter().map(|&c| r[c]).product::<f64>(),
                    "{k}-{p} {g}"
                );
            }
        }
        for a in 0..k {
            assert!(column(rows, a).iter().sum::<f64>().abs() < TOL);
            for b in (a + 1)..k {
                assert!(
                    dot(&column(rows, a), &column(rows, b)).abs() < TOL,
                    "{k}-{p}"
                );
            }
        }
    }
}

/// An N x (N-1) Plackett-Burman design with the all-ones column added is a
/// Hadamard matrix; any k <= N-1 columns of it are balanced and orthogonal.
#[test]
fn plackett_burman_columns_are_balanced_and_orthogonal_in_the_smallest_run_size() {
    for k in 1..=19usize {
        let d = plackett_burman(k).unwrap();
        let n = [8usize, 12, 16, 20].into_iter().find(|&n| n > k).unwrap();
        assert_eq!(d.data.len(), n, "k={k}");
        assert!(d
            .data
            .iter()
            .all(|r| r.len() == k && r.iter().all(|&v| v == 1.0 || v == -1.0)));
        for a in 0..k {
            let ca = column(&d.data, a);
            assert!(ca.iter().sum::<f64>().abs() < TOL, "k={k} col {a}");
            for b in (a + 1)..k {
                assert!(
                    dot(&ca, &column(&d.data, b)).abs() < TOL,
                    "k={k} cols {a},{b}"
                );
            }
        }
        if k == n - 1 {
            let h: Vec<Vec<f64>> = d
                .data
                .iter()
                .map(|r| std::iter::once(1.0).chain(r.iter().copied()).collect())
                .collect();
            for a in 0..n {
                for b in (a + 1)..n {
                    assert!(dot(&h[a], &h[b]).abs() < TOL, "k={k} rows {a},{b}");
                }
            }
        }
    }
}

/// An orthogonal array of strength two: for every pair of columns each pair of
/// levels appears the same number of times.
#[test]
fn taguchi_arrays_are_orthogonal_arrays_of_strength_two() {
    let arrays = [
        ("L4", 4usize),
        ("L8", 8),
        ("L9", 9),
        ("L12", 12),
        ("L16", 16),
        ("L18", 18),
        ("L27", 27),
    ];
    for (name, n) in arrays {
        let k = max_factors(name).unwrap();
        let d = taguchi_array(name, k).unwrap();
        assert_eq!(d.data.len(), n, "{name}");
        let levels: Vec<BTreeSet<i64>> = (0..k)
            .map(|j| column(&d.data, j).iter().map(|&v| v as i64).collect())
            .collect();
        for (j, l) in levels.iter().enumerate() {
            let ok = *l == BTreeSet::from([-1, 1]) || *l == BTreeSet::from([-1, 0, 1]);
            assert!(ok, "{name} col {j}: levels {l:?}");
        }
        for a in 0..k {
            for b in (a + 1)..k {
                let expect = n / (levels[a].len() * levels[b].len());
                for &x in &levels[a] {
                    for &y in &levels[b] {
                        let count = d
                            .data
                            .iter()
                            .filter(|r| r[a] as i64 == x && r[b] as i64 == y)
                            .count();
                        assert_eq!(count, expect, "{name} cols {a},{b} levels ({x},{y})");
                    }
                }
            }
        }
        for take in 1..k {
            let part = taguchi_array(name, take).unwrap();
            for (r, full) in part.data.iter().zip(&d.data) {
                assert_eq!(r[..], full[..take], "{name} k={take}");
            }
        }
    }
}

fn binomial(n: usize, r: usize) -> usize {
    (0..r).fold(1usize, |acc, i| acc * (n - i) / (i + 1))
}

/// {q, m} lattice: every composition of m into q parts, over m; simplex
/// centroid: the equal-share mixture of every non-empty subset.
#[test]
fn mixture_designs_are_the_compositions_and_the_subset_centroids() {
    for q in 2..=8usize {
        for m in 1..=6usize {
            let d = simplex_lattice(q, m).unwrap();
            assert_eq!(d.data.len(), binomial(q + m - 1, m), "q={q} m={m}");
            let keys: BTreeSet<Vec<i64>> = d
                .data
                .iter()
                .map(|r| {
                    assert!((r.iter().sum::<f64>() - 1.0).abs() < 1e-12);
                    r.iter()
                        .map(|&v| {
                            let s = v * m as f64;
                            assert!((s - s.round()).abs() < 1e-9, "q={q} m={m}: {v}");
                            s.round() as i64
                        })
                        .collect()
                })
                .collect();
            assert_eq!(keys.len(), d.data.len(), "q={q} m={m}: a point repeats");
            assert!(keys.iter().all(|k| k.iter().sum::<i64>() == m as i64));
        }
    }
    for q in 2..=12usize {
        let d = simplex_centroid(q).unwrap();
        assert_eq!(d.data.len(), (1 << q) - 1, "q={q}");
        let subsets: BTreeSet<Vec<bool>> = d
            .data
            .iter()
            .map(|r| {
                let on: Vec<f64> = r.iter().copied().filter(|&v| v > 0.0).collect();
                assert!(!on.is_empty());
                assert!(
                    on.iter()
                        .all(|&v| (v - 1.0 / on.len() as f64).abs() < 1e-12),
                    "q={q}"
                );
                assert!((r.iter().sum::<f64>() - 1.0).abs() < 1e-12);
                r.iter().map(|&v| v > 0.0).collect()
            })
            .collect();
        assert_eq!(subsets.len(), (1 << q) - 1, "q={q}: a subset is missing");
    }
}
