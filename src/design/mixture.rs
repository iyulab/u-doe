//! Mixture designs for formulation experiments.
//!
//! Mixture designs are used when the factors are proportions of a whole
//! (e.g. ingredient percentages) and must sum to 1.0 for every run.
//!
//! # Simplex Lattice Design {q, m}
//!
//! Generates all points on a (q-1)-dimensional simplex where each component
//! takes values in {0, 1/m, 2/m, ..., 1}.  Equivalent to enumerating all
//! compositions of the integer m into q non-negative parts and dividing by m.
//!
//! Run count: C(q + m - 1, m)
//!
//! Reference: Cornell (2002), *Experiments with Mixtures*, 3rd ed., Ch. 2
//!
//! # Simplex Centroid Design
//!
//! Generates the centroids of every non-empty subset of the q components
//! (2^q − 1 points total).  Each point gives equal weight to the components
//! in the subset and zero to the rest.
//!
//! Reference: Cornell (2002), Ch. 3

use super::DesignMatrix;

/// Generate a Simplex Lattice Design {q, m}.
///
/// Each component takes values from {0, 1/m, 2/m, ..., 1} subject to
/// the constraint that all components sum to 1.
///
/// # Arguments
/// * `q` — Number of mixture components (≥ 2).
/// * `m` — Degree of the lattice (≥ 1).
///
/// # Examples
///
/// ```
/// use u_doe::design::mixture::simplex_lattice;
/// let d = simplex_lattice(3, 2);
/// assert_eq!(d.run_count(), 6);
/// assert_eq!(d.factor_count(), 3);
/// ```
pub fn simplex_lattice(q: usize, m: usize) -> DesignMatrix {
    let mut rows: Vec<Vec<f64>> = Vec::new();
    let mut point = vec![0usize; q];
    generate_compositions(m, q, 0, &mut point, &mut rows, m);
    DesignMatrix {
        data: rows,
        factor_names: mixture_names(q),
    }
}

/// Generate a Simplex Centroid Design for `q` components.
///
/// Produces 2^q − 1 points: the centroids of every non-empty subset of the
/// q components.
///
/// # Arguments
/// * `q` — Number of mixture components (≥ 2).
///
/// # Examples
///
/// ```
/// use u_doe::design::mixture::simplex_centroid;
/// let d = simplex_centroid(3);
/// assert_eq!(d.run_count(), 7);
/// assert_eq!(d.factor_count(), 3);
/// ```
pub fn simplex_centroid(q: usize) -> DesignMatrix {
    let n_points = (1usize << q) - 1; // 2^q - 1
    let mut rows: Vec<Vec<f64>> = Vec::with_capacity(n_points);

    for mask in 1..=(n_points as u64) {
        let count = mask.count_ones() as f64;
        let row: Vec<f64> = (0..q)
            .map(|i| if (mask >> i) & 1 == 1 { 1.0 / count } else { 0.0 })
            .collect();
        rows.push(row);
    }

    DesignMatrix {
        data: rows,
        factor_names: mixture_names(q),
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Factor names X1, X2, ..., Xq for mixture components.
fn mixture_names(q: usize) -> Vec<String> {
    (1..=q).map(|i| format!("X{}", i)).collect()
}

/// Recursive enumeration of all compositions of `remaining` into
/// `q - index` non-negative integer parts.  Each leaf divides by `m`
/// to produce a proportion in [0, 1].
fn generate_compositions(
    remaining: usize,
    parts: usize,
    index: usize,
    point: &mut Vec<usize>,
    results: &mut Vec<Vec<f64>>,
    m: usize,
) {
    if index == parts - 1 {
        point[index] = remaining;
        let row: Vec<f64> = point.iter().map(|&v| v as f64 / m as f64).collect();
        results.push(row);
        return;
    }
    for v in 0..=remaining {
        point[index] = v;
        generate_compositions(remaining - v, parts, index + 1, point, results, m);
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn rows_sum_to_one(d: &DesignMatrix) {
        for (i, row) in d.data.iter().enumerate() {
            let s: f64 = row.iter().sum();
            assert!(
                (s - 1.0).abs() < 1e-10,
                "run {i} sums to {s}, expected 1.0"
            );
        }
    }

    #[test]
    fn simplex_lattice_3_2_run_count() {
        let d = simplex_lattice(3, 2);
        assert_eq!(d.run_count(), 6);
        assert_eq!(d.factor_count(), 3);
    }

    #[test]
    fn simplex_lattice_3_2_sums_to_one() {
        let d = simplex_lattice(3, 2);
        rows_sum_to_one(&d);
    }

    #[test]
    fn simplex_lattice_3_3_run_count() {
        let d = simplex_lattice(3, 3);
        assert_eq!(d.run_count(), 10);
    }

    #[test]
    fn simplex_centroid_3_run_count() {
        let d = simplex_centroid(3);
        assert_eq!(d.run_count(), 7);
        assert_eq!(d.factor_count(), 3);
    }

    #[test]
    fn simplex_centroid_3_sums_to_one() {
        let d = simplex_centroid(3);
        rows_sum_to_one(&d);
    }

    #[test]
    fn simplex_centroid_4_run_count() {
        let d = simplex_centroid(4);
        assert_eq!(d.run_count(), 15);
    }
}
