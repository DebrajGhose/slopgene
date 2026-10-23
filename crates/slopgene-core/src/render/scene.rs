//! A renderer-neutral drawing: a list of simple shapes, each optionally tagged with the thing
//! it depicts so an interactive view can tell what the pointer is over.

/// sRGB colour with straight (not premultiplied) alpha.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Rgba {
    pub r: u8,
    pub g: u8,
    pub b: u8,
    pub a: u8,
}

impl Rgba {
    pub const fn rgb(r: u8, g: u8, b: u8) -> Rgba {
        Rgba { r, g, b, a: 255 }
    }

    pub const fn new(r: u8, g: u8, b: u8, a: u8) -> Rgba {
        Rgba { r, g, b, a }
    }

    pub fn with_alpha(self, a: u8) -> Rgba {
        Rgba { a, ..self }
    }

    pub fn hex(self) -> String {
        format!("#{:02x}{:02x}{:02x}", self.r, self.g, self.b)
    }

    /// Mix toward `other` by `t` (0 = self, 1 = other).
    pub fn mix(self, other: Rgba, t: f32) -> Rgba {
        let lerp = |a: u8, b: u8| (a as f32 + (b as f32 - a as f32) * t).round().clamp(0.0, 255.0) as u8;
        Rgba { r: lerp(self.r, other.r), g: lerp(self.g, other.g), b: lerp(self.b, other.b), a: lerp(self.a, other.a) }
    }

    /// Relative luminance (0 dark .. 1 light), for picking readable text colours.
    pub fn luminance(self) -> f32 {
        let lin = |c: u8| {
            let c = c as f32 / 255.0;
            if c <= 0.04045 { c / 12.92 } else { ((c + 0.055) / 1.055).powf(2.4) }
        };
        0.2126 * lin(self.r) + 0.7152 * lin(self.g) + 0.0722 * lin(self.b)
    }
}

impl From<crate::model::Color> for Rgba {
    fn from(c: crate::model::Color) -> Rgba {
        Rgba::rgb(c.r, c.g, c.b)
    }
}

/// Horizontal alignment of text relative to its anchor point.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Anchor {
    Start,
    Middle,
    End,
}

/// One primitive. Coordinates are in pixels, y pointing down.
#[derive(Clone, Debug, PartialEq)]
pub enum Shape {
    /// An open stroked path (arcs are sampled into polylines). Ends are cut square.
    Polyline {
        points: Vec<[f32; 2]>,
        color: Rgba,
        width: f32,
    },
    /// A filled convex polygon.
    Polygon {
        points: Vec<[f32; 2]>,
        fill: Rgba,
        stroke: Option<(Rgba, f32)>,
    },
    Circle {
        center: [f32; 2],
        radius: f32,
        fill: Option<Rgba>,
        stroke: Option<(Rgba, f32)>,
    },
    /// Text with its baseline-centre (vertical middle) at `pos`.
    Text {
        pos: [f32; 2],
        text: String,
        size: f32,
        color: Rgba,
        anchor: Anchor,
        bold: bool,
    },
}

/// What a shape depicts.
#[derive(Clone, Debug, PartialEq, Default)]
pub enum Tag {
    #[default]
    None,
    Feature(u64),
    Primer(u64),
    /// Enzyme name and its top-strand cut position.
    Enzyme(String, usize),
    Backbone,
    Selection,
    /// An aligned sequence, by its index in [`super::MapOptions::alignments`].
    Track(usize),
    /// A difference: the track's index, then the mark's index in that track.
    Difference(usize, usize),
    /// A piece of a cut molecule, by its index in the ring's pieces.
    Piece(usize),
}

#[derive(Clone, Debug, PartialEq)]
pub struct Item {
    pub shape: Shape,
    pub tag: Tag,
}

/// A finished drawing.
#[derive(Clone, Debug, PartialEq)]
pub struct Scene {
    pub width: f32,
    pub height: f32,
    pub background: Rgba,
    pub items: Vec<Item>,
}

impl Scene {
    pub fn new(width: f32, height: f32, background: Rgba) -> Scene {
        Scene { width, height, background, items: Vec::new() }
    }

    pub fn push(&mut self, shape: Shape, tag: Tag) {
        self.items.push(Item { shape, tag });
    }
}

/// Measures text so layouts can avoid overlapping labels.
pub trait TextMeasure {
    fn text_width(&self, text: &str, size: f32, bold: bool) -> f32;
}

/// Width estimate for a proportional sans-serif font (about 0.55 em per character), used when
/// no real font metrics are available (SVG export).
pub struct ApproxMeasure;

impl TextMeasure for ApproxMeasure {
    fn text_width(&self, text: &str, size: f32, bold: bool) -> f32 {
        let per_char = if bold { 0.6 } else { 0.55 };
        text.chars().count() as f32 * size * per_char
    }
}

/// Points along a circular arc around `center`, from angle `a0` to `a1` (radians, clockwise
/// from 12 o'clock), roughly one point per 3 pixels of arc.
pub fn arc_points(center: [f32; 2], radius: f32, a0: f32, a1: f32) -> Vec<[f32; 2]> {
    let span = (a1 - a0).abs();
    let steps = ((span * radius / 3.0).ceil() as usize).clamp(2, 720);
    (0..=steps)
        .map(|ii| {
            let a = a0 + (a1 - a0) * ii as f32 / steps as f32;
            polar(center, radius, a)
        })
        .collect()
}

/// Point at `radius` and angle `a` (radians, clockwise from 12 o'clock).
pub fn polar(center: [f32; 2], radius: f32, a: f32) -> [f32; 2] {
    [center[0] + radius * a.sin(), center[1] - radius * a.cos()]
}
