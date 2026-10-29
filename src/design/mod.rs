//! Design generation modules.

pub mod box_behnken;
pub mod ccd;
pub mod dsd;
pub mod factorial;
pub mod mixture;
pub mod plackett_burman;
pub mod taguchi;

pub use dsd::definitive_screening;

/// Most center points a design takes. A count is a number of rows to build,
/// and an unbounded one exhausts memory before it is refused; center-point
/// replicates beyond a few dozen add no information a design uses.
pub const MAX_CENTER_POINTS: usize = 1_000;

use crate::error::DoeError;

/// The block of a response-surface design a run belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointType {
    /// A corner of the factorial cube (all coordinates at the factorial level).
    Factorial,
    /// An axial (star) point: one coordinate off the center, the rest at 0.
    Axial,
    /// A run at the center of the region (all coordinates 0).
    Center,
    /// An edge midpoint (Box-Behnken): two coordinates at ±1, the rest at 0.
    Edge,
}

impl PointType {
    /// The wire name (`"factorial"`, `"axial"`, `"center"`, `"edge"`).
    pub fn name(self) -> &'static str {
        match self {
            Self::Factorial => "factorial",
            Self::Axial => "axial",
            Self::Center => "center",
            Self::Edge => "edge",
        }
    }
}

/// A design matrix with the block each run was built in -- what the
/// construction knows and the coordinates only hint at.
#[derive(Debug, Clone)]
pub struct ClassifiedDesign {
    /// The runs.
    pub matrix: DesignMatrix,
    /// The type of each run, in the order of `matrix.data`.
    pub point_types: Vec<PointType>,
}

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

/// A design with its runs in the order they are to be performed.
///
/// The standard-order runs of a generated design are not the order to run them
/// in: time-ordered disturbances (a drifting instrument, a warming furnace)
/// would be confounded with the factors. Running in random order spreads them
/// over the factor settings (Montgomery, *Design and Analysis of Experiments*,
/// ch. 1: randomization).
#[derive(Debug, Clone)]
pub struct RandomizedDesign {
    /// The runs, in the order to perform them.
    pub design: DesignMatrix,
    /// `run_order[i]` is the position in the standard-order design of the run
    /// performed `i`-th (0-based): `design.data[i]` is the original's
    /// `data[run_order[i]]`. A response recorded in run order is put back in
    /// standard order by placing the `i`-th one at `run_order[i]`.
    pub run_order: Vec<usize>,
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

    /// The design with its runs in a random order, the same for the same `seed`.
    ///
    /// # Examples
    ///
    /// ```
    /// use u_doe::design::factorial::full_factorial;
    ///
    /// let design = full_factorial(3).unwrap();
    /// let planned = design.randomized(42);
    /// assert_eq!(planned.design.run_count(), 8);
    /// // The same seed, the same order.
    /// assert_eq!(planned.run_order, design.randomized(42).run_order);
    /// // Every standard-order run appears once.
    /// let mut sorted = planned.run_order.clone();
    /// sorted.sort_unstable();
    /// assert_eq!(sorted, (0..8).collect::<Vec<_>>());
    /// ```
    pub fn randomized(&self, seed: u64) -> RandomizedDesign {
        let mut rng = u_numflow::random::create_rng(seed);
        let run_order = u_numflow::random::shuffled_indices(self.run_count(), &mut rng);
        let data = run_order.iter().map(|&i| self.data[i].clone()).collect();
        RandomizedDesign {
            design: DesignMatrix {
                data,
                factor_names: self.factor_names.clone(),
            },
            run_order,
        }
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

    /// Refuses a design with no runs or no factors: there is nothing to
    /// analyse. It is reachable -- `data` and `factor_names` are public fields
    /// and the WebAssembly bindings fill them straight from caller input, so
    /// `[]` arrives as a design.
    ///
    /// # Errors
    /// [`DoeError::EmptyInput`] with `parameter` `design` when there are no
    /// runs, else `factor_names` when there are no factors.
    pub fn check_not_empty(&self) -> Result<(), DoeError> {
        if self.data.is_empty() {
            Err(DoeError::EmptyInput {
                parameter: "design",
            })
        } else if self.factor_names.is_empty() {
            Err(DoeError::EmptyInput {
                parameter: "factor_names",
            })
        } else {
            Ok(())
        }
    }

    /// Refuses a run whose length differs from the number of factor names.
    ///
    /// `data` and `factor_names` are public fields and the WebAssembly bindings
    /// fill them straight from caller input, so nothing else keeps them in
    /// step. Analyses index `data[run][factor]` for every named factor; a short
    /// row would panic there, and a long one would be read only in part.
    ///
    /// # Errors
    /// [`DoeError::DimensionMismatch`] naming the first such run.
    pub fn check_shape(&self) -> Result<(), DoeError> {
        let expected = self.factor_count();
        match self
            .data
            .iter()
            .enumerate()
            .find(|(_, row)| row.len() != expected)
        {
            Some((run, row)) => Err(DoeError::DimensionMismatch {
                parameter: "design",
                index: Some(run),
                expected,
                got: row.len(),
                counted: "one entry per factor",
            }),
            None => Ok(()),
        }
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
