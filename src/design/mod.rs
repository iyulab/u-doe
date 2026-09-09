//! Design generation modules.

pub mod box_behnken;
pub mod ccd;
pub mod dsd;
pub mod factorial;
pub mod mixture;
pub mod plackett_burman;
pub mod taguchi;

pub use dsd::definitive_screening;

/// DOE design matrix in coded units.
///
/// Each row represents one experimental run.
/// Each column represents one factor.
/// Values are coded: ±1 for 2-level, -1/0/+1 for 3-level.
#[derive(Debug, Clone)]
pub struct DesignMatrix {
    /// Row-major: `data[run][factor]`
    pub data: Vec<Vec<f64>>,
    /// Human-readable factor labels (e.g., "A", "B", "Temperature")
    pub factor_names: Vec<String>,
}

impl DesignMatrix {
    /// Number of experimental runs.
    pub fn run_count(&self) -> usize {
        self.data.len()
    }

    /// Number of factors.
    pub fn factor_count(&self) -> usize {
        self.factor_names.len()
    }

    /// Get the value of factor `j` in run `i`.
    pub fn get(&self, run: usize, factor: usize) -> f64 {
        self.data[run][factor]
    }

    /// First entry that is neither `-1` nor `+1`, as `(run, factor, value)`.
    ///
    /// Returns `None` when every entry is two-level coded. Effect estimation
    /// builds contrast columns by multiplying factor columns together, which is
    /// only defined for a two-level design — a design carrying centre points or
    /// axial points (central composite, Box-Behnken, definitive screening) has
    /// to be fitted by least squares instead.
    ///
    /// # Examples
    ///
    /// ```
    /// use u_doe::design::factorial::full_factorial;
    /// use u_doe::design::ccd::{ccd, AlphaType};
    ///
    /// assert!(full_factorial(2).unwrap().two_level_violation().is_none());
    /// assert!(ccd(2, AlphaType::FaceCentered, 1).unwrap().two_level_violation().is_some());
    /// ```
    pub fn two_level_violation(&self) -> Option<(usize, usize, f64)> {
        const TOL: f64 = 1e-9;
        for (run, row) in self.data.iter().enumerate() {
            for (factor, &value) in row.iter().enumerate() {
                if (value.abs() - 1.0).abs() > TOL {
                    return Some((run, factor, value));
                }
            }
        }
        None
    }

    /// Default factor names: "A", "B", "C", ...
    pub fn default_names(k: usize) -> Vec<String> {
        (0..k)
            .map(|i| {
                if i < 26 {
                    ((b'A' + i as u8) as char).to_string()
                } else {
                    format!("X{}", i + 1)
                }
            })
            .collect()
    }
}
