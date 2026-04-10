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
