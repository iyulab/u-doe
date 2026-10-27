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

/// This crate's version, as published (`CARGO_PKG_VERSION`) — for a consumer
/// that reports what it was built against without reading a lock file.
///
/// ```
/// assert!(!u_doe::VERSION.is_empty());
/// ```
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

pub mod analysis;
pub mod coding;
pub mod design;
pub mod error;
pub mod optimization;
pub mod power;

#[cfg(feature = "wasm")]
pub mod wasm;

#[cfg(all(test, feature = "wasm"))]
mod engine_card;

// The README's Rust examples are the first code most users copy, so they are
// compiled and run with the doc-tests. Without this they were checked by
// nothing.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
pub struct ReadmeDoctests;
