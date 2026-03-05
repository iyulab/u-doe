//! # u-doe
//!
//! Design of Experiments (DOE) framework for industrial and scientific experimentation.
//!
//! ## Modules
//!
//! - [`design`] — Design generation (factorial, PB, CCD, BBD, Taguchi)
//! - [`analysis`] — Effects estimation, ANOVA, RSM
//! - [`optimization`] — Desirability functions for multi-response optimization
//! - [`coding`] — Actual ↔ coded variable transformations
//! - [`error`] — Error types

pub mod analysis;
pub mod coding;
pub mod design;
pub mod error;
pub mod optimization;
pub mod power;
