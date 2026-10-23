//! Colours for maps and views.
//!
//! The categorical colours come from the user's figure palette (indigo, gold, jade, cyan, azure,
//! lime, spring green, deep navy, teal), so maps match the rest of their figures.

use super::scene::Rgba;

pub const INDIGO: Rgba = Rgba::rgb(0x51, 0x5F, 0xAA);
pub const GOLD: Rgba = Rgba::rgb(0xFF, 0xC3, 0x39);
pub const JADE: Rgba = Rgba::rgb(0x4B, 0xBD, 0x94);
pub const CYAN: Rgba = Rgba::rgb(0x0B, 0xB0, 0xD5);
pub const AZURE: Rgba = Rgba::rgb(0x0B, 0x8E, 0xCE);
pub const LIME: Rgba = Rgba::rgb(0xAC, 0xC8, 0x37);
pub const SPRING_GREEN: Rgba = Rgba::rgb(0x00, 0xE0, 0x70);
pub const NAVY: Rgba = Rgba::rgb(0x00, 0x40, 0x71);
pub const TEAL: Rgba = Rgba::rgb(0x0B, 0x80, 0x83);
pub const GRAY: Rgba = Rgba::rgb(0x7F, 0x7F, 0x7F);
pub const GRAY_LIGHT: Rgba = Rgba::rgb(0x99, 0x99, 0x99);
pub const INK: Rgba = Rgba::rgb(0x22, 0x22, 0x22);

/// The categorical palette in its preferred order.
pub const PALETTE: [Rgba; 9] = [INDIGO, GOLD, JADE, CYAN, AZURE, LIME, SPRING_GREEN, NAVY, TEAL];

// Colours for comparing aligned sequences. Green says the aligned sequence has the reference's
// bases and red says it lacks them, which is how people read a comparison at a glance. Orange
// marks bases that differ and violet bases only the aligned sequence has. They lie outside the
// categorical palette on purpose: each palette colour already stands for a feature type, and a
// difference must never pass for a feature. Every kind is also drawn with its own shape and
// label, so colour is never the only cue.
/// Reference bases the aligned sequence has too.
pub const ALIGNED: Rgba = Rgba::rgb(0x2E, 0x9E, 0x5B);
/// A base that differs.
pub const CHANGED: Rgba = Rgba::rgb(0xE6, 0x7E, 0x22);
/// Bases only the aligned sequence has.
pub const INSERTED: Rgba = Rgba::rgb(0x8E, 0x44, 0xAD);
/// Reference bases the aligned sequence lacks.
pub const MISSING: Rgba = Rgba::rgb(0xD6, 0x27, 0x28);

// Colours for the pieces of a cloning protocol, as the user sketched them: green for the bases
// that go on into the product and red for the bases left over. Like the comparison colours they
// lie outside the categorical palette, which stands for feature types.
/// Bases that go on into the product.
pub const KEPT: Rgba = Rgba::rgb(0x2E, 0x9E, 0x5B);
/// Bases left over: the rest of a donor, or the part of a vector an insert replaces.
pub const LEFT_OVER: Rgba = Rgba::rgb(0xD6, 0x27, 0x28);

/// Colours for one appearance (light or dark).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Theme {
    pub dark: bool,
    pub background: Rgba,
    /// Main text.
    pub ink: Rgba,
    /// Secondary text (lengths, positions).
    pub muted: Rgba,
    /// Faint lines (grids, ticks).
    pub faint: Rgba,
    pub backbone: Rgba,
    pub enzyme: Rgba,
    pub primer: Rgba,
    pub selection: Rgba,
    pub accent: Rgba,
    /// Where an aligned sequence has the reference's bases.
    pub aligned: Rgba,
    /// Differences of each kind (see [`difference_color`]).
    pub changed: Rgba,
    pub inserted: Rgba,
    pub missing: Rgba,
    /// The pieces of a protocol that go on into the product, and those left over.
    pub kept: Rgba,
    pub left_over: Rgba,
}

impl Theme {
    pub fn light() -> Theme {
        Theme {
            dark: false,
            background: Rgba::rgb(0xFC, 0xFC, 0xFB),
            ink: INK,
            muted: Rgba::rgb(0x6B, 0x6B, 0x6B),
            faint: Rgba::rgb(0xDD, 0xDD, 0xDA),
            backbone: Rgba::rgb(0x8A, 0x8A, 0x8A),
            enzyme: Rgba::rgb(0x44, 0x44, 0x44),
            primer: TEAL,
            selection: INDIGO.with_alpha(56),
            accent: INDIGO,
            aligned: ALIGNED,
            changed: CHANGED,
            inserted: INSERTED,
            missing: MISSING,
            kept: KEPT,
            left_over: LEFT_OVER,
        }
    }

    pub fn dark() -> Theme {
        Theme {
            dark: true,
            background: Rgba::rgb(0x1B, 0x1C, 0x1F),
            ink: Rgba::rgb(0xE8, 0xE8, 0xE6),
            muted: Rgba::rgb(0x9A, 0x9A, 0x9A),
            faint: Rgba::rgb(0x3A, 0x3B, 0x3F),
            backbone: Rgba::rgb(0x8C, 0x8C, 0x8C),
            enzyme: Rgba::rgb(0xC8, 0xC8, 0xC8),
            primer: CYAN,
            selection: Rgba::rgb(0x8A, 0x96, 0xE0).with_alpha(64),
            accent: Rgba::rgb(0x8A, 0x96, 0xE0),
            // The difference colours, lightened to stand out on the dark background.
            aligned: Rgba::rgb(0x5C, 0xC8, 0x82),
            changed: Rgba::rgb(0xF2, 0x9B, 0x4E),
            inserted: Rgba::rgb(0xB7, 0x7D, 0xD6),
            missing: Rgba::rgb(0xF0, 0x5A, 0x5B),
            kept: Rgba::rgb(0x5C, 0xC8, 0x82),
            left_over: Rgba::rgb(0xF0, 0x5A, 0x5B),
        }
    }
}

/// The colour of a difference of this kind in this theme.
pub fn difference_color(kind: super::DifferenceKind, theme: &Theme) -> Rgba {
    match kind {
        super::DifferenceKind::Changed => theme.changed,
        super::DifferenceKind::Inserted => theme.inserted,
        super::DifferenceKind::Missing => theme.missing,
    }
}

/// A difference colour for text: a little darker on a light background, where the bright colour
/// alone would be hard to read.
pub fn difference_ink(kind: super::DifferenceKind, theme: &Theme) -> Rgba {
    let color = difference_color(kind, theme);
    if theme.dark { color } else { color.mix(theme.ink, 0.25) }
}

impl Default for Theme {
    fn default() -> Self {
        Theme::light()
    }
}

/// Default colour for a feature type, used when the file gives none.
pub fn feature_color(kind: &str) -> Rgba {
    match kind.to_ascii_lowercase().as_str() {
        "cds" | "gene" | "mrna" | "exon" => GOLD,
        "promoter" | "-35_signal" | "-10_signal" => JADE,
        "terminator" | "polya_signal" | "polya_site" => NAVY,
        "rep_origin" | "origin_of_replication" | "orit" => AZURE,
        "primer_bind" | "primer" => CYAN,
        "protein_bind" | "misc_binding" | "binding_site" => TEAL,
        "enhancer" | "regulatory" | "misc_signal" | "operator" => LIME,
        "rbs" | "ribosome_binding_site" | "5'utr" | "3'utr" => SPRING_GREEN,
        "ltr" | "repeat_region" | "mobile_element" | "misc_recomb" | "long_terminal_repeat" => INDIGO,
        "sig_peptide" | "mat_peptide" | "transit_peptide" => TEAL,
        _ => GRAY_LIGHT,
    }
}
