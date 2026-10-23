//! `slopgene render`: draw a map as SVG or PNG.

use crate::enzyme_cmd::choose_enzymes;
use crate::util::*;
use anyhow::{Result, anyhow};
use serde_json::json;
use slopgene_core::render::{self, ApproxMeasure, MapOptions, Scene, Theme};
use std::path::{Path, PathBuf};

/// Write a scene as SVG, or as PNG (rasterised at `scale`) when the path ends in `.png`.
pub fn save_scene(scene: &Scene, path: &Path, scale: f32) -> Result<()> {
    slopgene_script::files::save_scene(scene, path, scale).map_err(|e| anyhow!(e))
}

#[allow(clippy::too_many_arguments)]
pub fn run(
    path: &str,
    output: Option<&str>,
    width: f32,
    height: f32,
    enzymes: &str,
    dark: bool,
    scale: f32,
    font_size: f32,
    linear: bool,
    colors: &str,
) -> Result<Output> {
    let colors = render::ColorMode::parse(colors).ok_or_else(|| anyhow!("--colors must be palette or file"))?;
    let input = PathBuf::from(path);
    let mut doc = load(&input)?;
    if linear {
        doc.topology = slopgene_core::Topology::Linear;
    }
    let marks = if enzymes.eq_ignore_ascii_case("none") {
        Vec::new()
    } else {
        let db = enzyme_db();
        let chosen = choose_enzymes(db, enzymes)?;
        render::unique_cutter_marks(&doc, &chosen)
    };
    let options =
        MapOptions { width, height, theme: if dark { Theme::dark() } else { Theme::light() }, font_size, enzymes: marks, colors, ..MapOptions::default() };
    let scene = render::render_map(&doc, &options, &ApproxMeasure);
    let target = output.map(PathBuf::from).unwrap_or_else(|| input.with_extension("svg"));
    save_scene(&scene, &target, scale)?;
    let value = json!({ "output": path_text(&target), "width": scene.width, "height": scene.height, "enzymes_marked": options.enzymes.len() });
    Ok(Output::new(value, format!("wrote {}", path_text(&target))))
}
