//! Actual-to-coded variable transformations.
//!
//! Converts factor levels between engineering units (actual) and
//! dimensionless coded values (−1 to +1).
//!
//! # Formula
//!
//! `coded = (actual − center) / half_range`
//! where `center = (high + low) / 2`, `half_range = (high − low) / 2`
//!
//! # Examples
//!
//! ```
//! use u_doe::coding::Encoder;
//!
//! let enc = Encoder::new(100.0, 200.0);
//! assert!((enc.encode(150.0) - 0.0).abs() < 1e-10);
//! assert!((enc.encode(100.0) - (-1.0)).abs() < 1e-10);
//! assert!((enc.decode(1.0) - 200.0).abs() < 1e-10);
//! ```
//!
//! Reference: Montgomery (2019), *Introduction to Statistical Quality Control*, 8th ed., Section 5.3.1.

use crate::error::DoeError;

/// Linear encoder between actual and coded values for one factor.
///
/// Maps `[low, high]` linearly to `[-1, +1]`.
#[derive(Debug, Clone, Copy)]
pub struct Encoder {
    center: f64,
    half_range: f64,
}

impl Encoder {
    /// Create an encoder.
    ///
    /// # Panics
    ///
    /// Panics if `low >= high`. Use [`try_new`] for fallible creation.
    ///
    /// [`try_new`]: Encoder::try_new
    pub fn new(low: f64, high: f64) -> Self {
        Self::try_new(low, high).expect("low must be strictly less than high")
    }

    /// Create an encoder, returning `Err` if `low >= high`.
    pub fn try_new(low: f64, high: f64) -> Result<Self, DoeError> {
        if low >= high {
            return Err(DoeError::InvalidSpecification(format!(
                "low ({low}) must be strictly less than high ({high})"
            )));
        }
        Ok(Self {
            center: (high + low) / 2.0,
            half_range: (high - low) / 2.0,
        })
    }

    /// Convert an actual value to its coded equivalent.
    ///
    /// Returns `−1.0` for `low`, `0.0` for center, `+1.0` for `high`.
    pub fn encode(&self, actual: f64) -> f64 {
        (actual - self.center) / self.half_range
    }

    /// Convert a coded value back to actual units.
    pub fn decode(&self, coded: f64) -> f64 {
        coded * self.half_range + self.center
    }

    /// The low level in actual units (encoded as −1).
    pub fn low(&self) -> f64 {
        self.center - self.half_range
    }

    /// The high level in actual units (encoded as +1).
    pub fn high(&self) -> f64 {
        self.center + self.half_range
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_center_is_zero() {
        let enc = Encoder::new(50.0, 70.0);
        let coded = enc.encode(60.0);
        assert!((coded - 0.0).abs() < 1e-10);
    }

    #[test]
    fn encode_low_is_minus_one() {
        let enc = Encoder::new(100.0, 200.0);
        assert!((enc.encode(100.0) - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn encode_high_is_plus_one() {
        let enc = Encoder::new(100.0, 200.0);
        assert!((enc.encode(200.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn decode_roundtrip() {
        let enc = Encoder::new(50.0, 70.0);
        let actual = 63.5;
        let roundtrip = enc.decode(enc.encode(actual));
        assert!((roundtrip - actual).abs() < 1e-10);
    }

    #[test]
    fn try_new_invalid_range() {
        assert!(Encoder::try_new(5.0, 5.0).is_err());
        assert!(Encoder::try_new(10.0, 5.0).is_err());
    }

    #[test]
    fn try_new_valid() {
        assert!(Encoder::try_new(0.0, 100.0).is_ok());
    }

    #[test]
    fn low_high_accessors() {
        let enc = Encoder::new(100.0, 200.0);
        assert!((enc.low() - 100.0).abs() < 1e-10);
        assert!((enc.high() - 200.0).abs() < 1e-10);
    }
}
