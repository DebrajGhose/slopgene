//! Write a [`Scene`] as SVG.

use super::scene::{Anchor, Rgba, Scene, Shape};
use crate::xml::escape_attr;

const FONT_STACK: &str = "Inter, 'Segoe UI', -apple-system, 'Helvetica Neue', Arial, sans-serif";

fn color_attr(name: &str, c: Rgba) -> String {
    if c.a == 255 { format!(" {name}=\"{}\"", c.hex()) } else { format!(" {name}=\"{}\" {name}-opacity=\"{:.3}\"", c.hex(), c.a as f32 / 255.0) }
}

fn points_attr(points: &[[f32; 2]]) -> String {
    points.iter().map(|p| format!("{:.1},{:.1}", p[0], p[1])).collect::<Vec<_>>().join(" ")
}

/// Serialise the scene as a standalone SVG document.
pub fn to_svg(scene: &Scene) -> String {
    let mut out = String::new();
    out.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}\" height=\"{h}\" viewBox=\"0 0 {w} {h}\" font-family=\"{FONT_STACK}\">\n",
        w = scene.width.round(),
        h = scene.height.round()
    ));
    out.push_str(&format!("<rect width=\"100%\" height=\"100%\"{}/>\n", color_attr("fill", scene.background)));
    for item in &scene.items {
        match &item.shape {
            Shape::Polyline { points, color, width } => {
                out.push_str(&format!(
                    "<polyline points=\"{}\" fill=\"none\"{} stroke-width=\"{width:.2}\"/>\n",
                    points_attr(points),
                    color_attr("stroke", *color)
                ));
            }
            Shape::Polygon { points, fill, stroke } => {
                let stroke = stroke.map(|(c, w)| format!("{} stroke-width=\"{w:.2}\" stroke-linejoin=\"round\"", color_attr("stroke", c))).unwrap_or_default();
                out.push_str(&format!("<polygon points=\"{}\"{}{stroke}/>\n", points_attr(points), color_attr("fill", *fill)));
            }
            Shape::Circle { center, radius, fill, stroke } => {
                let fill = fill.map(|c| color_attr("fill", c)).unwrap_or_else(|| " fill=\"none\"".into());
                let stroke = stroke.map(|(c, w)| format!("{} stroke-width=\"{w:.2}\"", color_attr("stroke", c))).unwrap_or_default();
                out.push_str(&format!("<circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"{radius:.1}\"{fill}{stroke}/>\n", center[0], center[1]));
            }
            Shape::Text { pos, text, size, color, anchor, bold } => {
                let anchor = match anchor {
                    Anchor::Start => "start",
                    Anchor::Middle => "middle",
                    Anchor::End => "end",
                };
                let weight = if *bold { " font-weight=\"600\"" } else { "" };
                out.push_str(&format!(
                    "<text x=\"{:.1}\" y=\"{:.1}\" font-size=\"{size:.1}\" text-anchor=\"{anchor}\" dominant-baseline=\"central\"{weight}{}>{}</text>\n",
                    pos[0],
                    pos[1],
                    color_attr("fill", *color),
                    escape_attr(text)
                ));
            }
        }
    }
    out.push_str("</svg>\n");
    out
}
