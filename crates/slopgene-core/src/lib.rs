//! SlopGene core: the sequence document model, file formats and molecular-biology
//! algorithms shared by the command line, the scripting engine and the desktop app.
//!
//! Coordinates are 0-based and half-open inside this crate ([`Span`]). Every user-facing
//! surface (CLI, JSON, scripts, GUI) shows 1-based inclusive coordinates, as SnapGene and
//! GenBank do; convert with [`Span::from_1based`] and [`Span::to_1based`].

pub mod align;
pub mod alphabet;
pub mod annotate;
pub mod clock;
pub mod cloning;
pub mod codon;
pub mod digest;
pub mod edit;
pub mod enzymes;
pub mod error;
pub mod excise;
pub mod gel;
pub mod io;
pub mod json;
pub mod model;
pub mod pcr;
pub mod primer_design;
pub mod primers;
pub mod protocol;
pub mod render;
pub mod search;
pub mod span;
pub mod stats;
pub mod translate;
pub mod xml;

pub use error::{Error, Result};
pub use model::*;
pub use span::Span;
