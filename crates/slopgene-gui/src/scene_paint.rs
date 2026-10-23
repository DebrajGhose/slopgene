//! Paint a core [`Scene`] with egui, and find which tagged shape lies under the pointer.
//!
//! Maps are laid out once by `slopgene_core::render` (the same code that writes SVG files),
//! so the app and exported images match.

use crate::style;
use eframe::egui::{self, Align2, Color32, Painter, Pos2, Shape as EShape, Stroke, pos2};
use slopgene_core::render::{Anchor, Rgba, Scene, Shape, Tag, TextMeasure};

pub fn color(c: Rgba) -> Color32 {
    Color32::from_rgba_unmultiplied(c.r, c.g, c.b, c.a)
}

/// Measures text with the app's real fonts, so map labels are spaced exactly.
pub struct EguiMeasure<'a> {
    pub ctx: &'a egui::Context,
}

impl TextMeasure for EguiMeasure<'_> {
    fn text_width(&self, text: &str, size: f32, bold: bool) -> f32 {
        let font = if bold { style::bold_font(size) } else { style::ui_font(size) };
        self.ctx.fonts_mut(|f| f.layout_no_wrap(text.to_string(), font, Color32::WHITE).size().x)
    }
}

/// Draw the scene with its top-left corner at `origin`.
pub fn paint(painter: &Painter, scene: &Scene, origin: Pos2) {
    let at = |p: [f32; 2]| pos2(origin.x + p[0], origin.y + p[1]);
    for item in &scene.items {
        match &item.shape {
            Shape::Polyline { points, color: c, width } => {
                let pts: Vec<Pos2> = points.iter().map(|&p| at(p)).collect();
                painter.add(EShape::line(pts, Stroke::new(*width, color(*c))));
            }
            Shape::Polygon { points, fill, stroke } => {
                let pts: Vec<Pos2> = points.iter().map(|&p| at(p)).collect();
                let stroke = stroke.map(|(c, w)| Stroke::new(w, color(c))).unwrap_or(Stroke::NONE);
                painter.add(EShape::convex_polygon(pts, color(*fill), stroke));
            }
            Shape::Circle { center, radius, fill, stroke } => {
                if let Some(f) = fill {
                    painter.circle_filled(at(*center), *radius, color(*f));
                }
                if let Some((c, w)) = stroke {
                    painter.circle_stroke(at(*center), *radius, Stroke::new(*w, color(*c)));
                }
            }
            Shape::Text { pos, text, size, color: c, anchor, bold } => {
                let align = match anchor {
                    Anchor::Start => Align2::LEFT_CENTER,
                    Anchor::Middle => Align2::CENTER_CENTER,
                    Anchor::End => Align2::RIGHT_CENTER,
                };
                let font = if *bold { style::bold_font(*size) } else { style::ui_font(*size) };
                painter.text(at(*pos), align, text, font, color(*c));
            }
        }
    }
}

/// Distance from `p` to the segment `a`-`b`.
fn segment_distance(p: [f32; 2], a: [f32; 2], b: [f32; 2]) -> f32 {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let len2 = dx * dx + dy * dy;
    let t = if len2 == 0.0 { 0.0 } else { (((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / len2).clamp(0.0, 1.0) };
    let (cx, cy) = (a[0] + t * dx, a[1] + t * dy);
    ((p[0] - cx).powi(2) + (p[1] - cy).powi(2)).sqrt()
}

/// Whether `p` lies inside a convex polygon (either winding).
fn inside_convex(p: [f32; 2], points: &[[f32; 2]]) -> bool {
    let mut sign = 0.0f32;
    for ii in 0..points.len() {
        let a = points[ii];
        let b = points[(ii + 1) % points.len()];
        let cross = (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]);
        if cross.abs() < f32::EPSILON {
            continue;
        }
        if sign == 0.0 {
            sign = cross.signum();
        } else if cross.signum() != sign {
            return false;
        }
    }
    true
}

/// The tag of the topmost tagged shape under `p` (scene coordinates). Untagged shapes,
/// the backbone and the selection band are ignored.
pub fn hit(scene: &Scene, p: [f32; 2], measure: &dyn TextMeasure) -> Option<Tag> {
    for item in scene.items.iter().rev() {
        if matches!(item.tag, Tag::None | Tag::Backbone | Tag::Selection) {
            continue;
        }
        let over = match &item.shape {
            Shape::Polyline { points, width, .. } => {
                let reach = (width / 2.0).max(3.0) + 1.0;
                points.windows(2).any(|w| segment_distance(p, w[0], w[1]) <= reach)
            }
            Shape::Polygon { points, .. } => inside_convex(p, points),
            Shape::Circle { center, radius, .. } => ((p[0] - center[0]).powi(2) + (p[1] - center[1]).powi(2)).sqrt() <= *radius,
            Shape::Text { pos, text, size, anchor, bold, .. } => {
                let w = measure.text_width(text, *size, *bold);
                let x0 = match anchor {
                    Anchor::Start => pos[0],
                    Anchor::Middle => pos[0] - w / 2.0,
                    Anchor::End => pos[0] - w,
                };
                p[0] >= x0 - 2.0 && p[0] <= x0 + w + 2.0 && (p[1] - pos[1]).abs() <= size * 0.7
            }
        };
        if over {
            return Some(item.tag.clone());
        }
    }
    None
}
