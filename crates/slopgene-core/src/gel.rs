//! Virtual agarose gels.
//!
//! Migration follows the usual log-linear relation between fragment size and distance, over a
//! size range that shifts with agarose percentage (higher percentages resolve smaller
//! fragments). Supercoiled plasmid runs faster than its linear size and nicked circles slower.
//! Band brightness follows DNA mass, so in an equimolar digest larger fragments look brighter.

use crate::render::scene::{Anchor, Rgba, Scene, Shape, Tag};

/// A molecular-weight marker.
#[derive(Clone, Debug, PartialEq)]
pub struct Ladder {
    pub name: &'static str,
    /// Band sizes in bp, largest first.
    pub sizes: &'static [usize],
}

/// Common commercial ladders (band sizes from the manufacturers' charts).
pub static LADDERS: &[Ladder] = &[
    Ladder { name: "1kb", sizes: &[10002, 8001, 6001, 5001, 4001, 3001, 2000, 1500, 1000, 517, 500] },
    Ladder { name: "1kb-plus", sizes: &[10002, 8001, 6001, 5001, 4001, 3001, 2000, 1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100] },
    Ladder { name: "100bp", sizes: &[1517, 1200, 1000, 900, 800, 700, 600, 517, 500, 400, 300, 200, 100] },
    Ladder { name: "generuler-1kb", sizes: &[10000, 8000, 6000, 5000, 4000, 3500, 3000, 2500, 2000, 1500, 1000, 750, 500, 250] },
    Ladder { name: "generuler-1kb-plus", sizes: &[20000, 10000, 7000, 5000, 4000, 3000, 2000, 1500, 1000, 700, 500, 400, 300, 200, 75] },
];

impl Ladder {
    pub fn get(name: &str) -> Option<&'static Ladder> {
        let key = name.trim().to_ascii_lowercase().replace([' ', '_'], "-");
        LADDERS.iter().find(|l| l.name == key)
    }
}

/// The physical form of a DNA molecule, which changes how fast it runs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Form {
    Linear,
    Supercoiled,
    Nicked,
}

/// One band in a lane.
#[derive(Clone, Debug, PartialEq)]
pub struct Band {
    /// Length in bp.
    pub size: usize,
    pub form: Form,
    /// Relative number of molecules (1 for each fragment of a digest).
    pub molar: f32,
}

/// One lane of the gel.
#[derive(Clone, Debug, PartialEq)]
pub struct Lane {
    pub label: String,
    pub bands: Vec<Band>,
    /// Marker lanes get size labels.
    pub ladder: bool,
}

impl Lane {
    pub fn ladder(ladder: &Ladder) -> Lane {
        Lane { label: ladder.name.to_string(), bands: ladder.sizes.iter().map(|&s| Band { size: s, form: Form::Linear, molar: 1.0 }).collect(), ladder: true }
    }

    /// The fragments of a digest, one band each (identical sizes stack up).
    pub fn digest(label: impl Into<String>, sizes: &[usize]) -> Lane {
        Lane { label: label.into(), bands: sizes.iter().map(|&s| Band { size: s, form: Form::Linear, molar: 1.0 }).collect(), ladder: false }
    }

    /// An uncut plasmid: mostly supercoiled with some nicked circles.
    pub fn uncut_plasmid(label: impl Into<String>, size: usize) -> Lane {
        Lane {
            label: label.into(),
            bands: vec![Band { size, form: Form::Supercoiled, molar: 1.0 }, Band { size, form: Form::Nicked, molar: 0.3 }],
            ladder: false,
        }
    }
}

/// Relative migration (0 = at the well, 1 = bottom of the gel) of a band.
pub fn migration(band: &Band, agarose_percent: f64) -> f32 {
    // Apparent linear size: supercoils are compact and run ahead, nicked circles lag behind.
    let apparent = band.size as f64
        * match band.form {
            Form::Linear => 1.0,
            Form::Supercoiled => 0.65,
            Form::Nicked => 1.6,
        };
    let p = agarose_percent.clamp(0.3, 4.0);
    let min_bp = 250.0 / p.powf(1.5);
    let max_bp = 12_000.0 / p.powf(1.6);
    let t = (max_bp.log10() - apparent.max(1.0).log10()) / (max_bp.log10() - min_bp.log10());
    // Compress the extremes: very large fragments pile up near the well, small ones near the front.
    let t = t.clamp(-0.6, 1.6);
    let squashed = if t < 0.0 {
        t * 0.1
    } else if t > 1.0 {
        1.0 + (t - 1.0) * 0.15
    } else {
        t
    };
    (0.06 + 0.86 * squashed).clamp(0.02, 0.98) as f32
}

/// Options for drawing a gel.
#[derive(Clone, Debug)]
pub struct GelOptions {
    pub agarose_percent: f64,
    pub lane_width: f32,
    pub height: f32,
    pub font_size: f32,
}

impl Default for GelOptions {
    fn default() -> Self {
        GelOptions { agarose_percent: 1.0, lane_width: 56.0, height: 420.0, font_size: 11.0 }
    }
}

/// Draw the gel. The gel itself is dark with glowing bands, like a gel imager picture.
pub fn render_gel(lanes: &[Lane], options: &GelOptions) -> Scene {
    let font = options.font_size;
    let gap = 14.0;
    let left = 64.0;
    let top = font * 3.0;
    let gel_h = options.height;
    let width = left + lanes.len() as f32 * (options.lane_width + gap) + gap + 16.0;
    let height = top + gel_h + font * 2.0;
    let page = Rgba::rgb(0xFA, 0xFA, 0xFA);
    let gel = Rgba::rgb(0x1E, 0x21, 0x27);
    let glow = Rgba::rgb(0xF2, 0xF6, 0xFF);
    let ink = Rgba::rgb(0x22, 0x22, 0x22);
    let muted = Rgba::rgb(0x6B, 0x6B, 0x6B);
    let mut scene = Scene::new(width, height, page);

    let gel_x0 = left - 8.0;
    let gel_x1 = width - 16.0;
    scene.push(Shape::Polygon { points: vec![[gel_x0, top], [gel_x1, top], [gel_x1, top + gel_h], [gel_x0, top + gel_h]], fill: gel, stroke: None }, Tag::None);
    scene.push(
        Shape::Text {
            pos: [gel_x0 - 6.0, top - font * 0.9],
            text: format!("{}%", options.agarose_percent),
            size: font * 0.9,
            color: muted,
            anchor: Anchor::End,
            bold: false,
        },
        Tag::None,
    );

    for (ii, lane) in lanes.iter().enumerate() {
        let x0 = left + gap + ii as f32 * (options.lane_width + gap);
        let x1 = x0 + options.lane_width;
        let mid = (x0 + x1) / 2.0;
        // Lane label above the gel, and the loading well.
        scene.push(
            Shape::Text { pos: [mid, top - font * 0.9], text: lane.label.clone(), size: font, color: ink, anchor: Anchor::Middle, bold: false },
            Tag::None,
        );
        scene.push(
            Shape::Polygon {
                points: vec![[x0, top + 6.0], [x1, top + 6.0], [x1, top + 10.0], [x0, top + 10.0]],
                fill: Rgba::rgb(0x0E, 0x10, 0x13),
                stroke: None,
            },
            Tag::None,
        );
        // Bands: brightness follows mass (size x molar amount), relative to the lane's heaviest band.
        let mut last_label_y: Option<f32> = None;
        let max_mass = lane.bands.iter().map(|b| b.size as f32 * b.molar).fold(1.0f32, f32::max);
        for band in &lane.bands {
            let y = top + 14.0 + migration(band, options.agarose_percent) * (gel_h - 24.0);
            let mass = band.size as f32 * band.molar;
            let strength = if lane.ladder { 0.75 } else { 0.35 + 0.65 * (mass / max_mass).powf(0.6) };
            let alpha = (strength * 255.0) as u8;
            // A soft halo under a sharp core reads as a glowing band.
            for (pad, a) in [(2.5f32, alpha / 4), (0.0, alpha)] {
                let h = 2.2 + pad;
                scene.push(
                    Shape::Polygon {
                        points: vec![[x0 + 3.0 - pad, y - h], [x1 - 3.0 + pad, y - h], [x1 - 3.0 + pad, y + h], [x0 + 3.0 - pad, y + h]],
                        fill: glow.with_alpha(a),
                        stroke: None,
                    },
                    Tag::None,
                );
            }
            // Size labels for the first ladder, skipping any too close to the previous one.
            if lane.ladder && ii == 0 && last_label_y.is_none_or(|prev| y - prev >= font * 1.05) {
                last_label_y = Some(y);
                let label = if band.size >= 1000 { format!("{} kb", trim_kb(band.size)) } else { format!("{} bp", band.size) };
                scene.push(Shape::Text { pos: [gel_x0 - 6.0, y], text: label, size: font * 0.9, color: muted, anchor: Anchor::End, bold: false }, Tag::None);
            }
        }
    }
    scene
}

/// 1500 -> "1.5", 10002 -> "10".
fn trim_kb(bp: usize) -> String {
    let kb = (bp as f64 / 1000.0 * 10.0).round() / 10.0;
    if (kb - kb.round()).abs() < 1e-9 { format!("{}", kb as i64) } else { format!("{kb:.1}") }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bigger_runs_slower_and_percentages_shift() {
        let band = |size| Band { size, form: Form::Linear, molar: 1.0 };
        assert!(migration(&band(500), 1.0) > migration(&band(5000), 1.0));
        // A 500 bp fragment travels less far in a 2% gel than in a 0.8% gel.
        assert!(migration(&band(500), 2.0) < migration(&band(500), 0.8));
        let sc = Band { size: 3000, form: Form::Supercoiled, molar: 1.0 };
        assert!(migration(&sc, 1.0) > migration(&band(3000), 1.0));
        assert_eq!(trim_kb(1500), "1.5");
        assert_eq!(trim_kb(10002), "10");
    }

    #[test]
    fn gel_renders() {
        let lanes = vec![Lane::ladder(Ladder::get("1kb").unwrap()), Lane::digest("EcoRI", &[3000, 1200]), Lane::uncut_plasmid("uncut", 4200)];
        let scene = render_gel(&lanes, &GelOptions::default());
        let svg = crate::render::svg::to_svg(&scene);
        assert!(svg.contains("10 kb"));
        assert!(svg.contains("EcoRI"));
    }
}
