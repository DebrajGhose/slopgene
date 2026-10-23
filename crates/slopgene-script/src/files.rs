//! File-system helpers shared by scripts, the command line and the desktop app: globbing,
//! listing sequence files, reading and writing text, and saving rendered scenes as SVG or PNG.
//!
//! Patterns are expanded here rather than by a shell, so they behave the same in PowerShell,
//! cmd, bash and zsh. Paths come back with forward slashes on every platform.

use rhai::{Array, Dynamic, Engine, EvalAltResult};
use slopgene_core::render::{Scene, svg::to_svg};
use std::path::Path;

type R<T> = Result<T, Box<EvalAltResult>>;

/// Extensions of the sequence files that folder listings pick up.
pub const SEQUENCE_EXTENSIONS: &[&str] = &["dna", "prot", "rna", "gb", "gbk", "genbank", "gbff", "fa", "fasta", "fna", "faa", "fas", "fsa", "json", "ape"];

/// Whether a path has one of the [`SEQUENCE_EXTENSIONS`], in any case.
pub fn is_sequence_file(path: &Path) -> bool {
    path.extension().and_then(|e| e.to_str()).is_some_and(|e| SEQUENCE_EXTENSIONS.contains(&e.to_ascii_lowercase().as_str()))
}

/// Whether a command-line argument is a glob pattern rather than a plain path.
pub fn has_glob_chars(text: &str) -> bool {
    text.contains(['*', '?', '['])
}

fn text(path: &Path) -> String {
    path.to_string_lossy().replace('\\', "/")
}

/// Files matching a glob pattern (`**` crosses folders), sorted. A pattern without wildcards
/// names one file, which is returned when it exists.
pub fn glob(pattern: &str) -> Result<Vec<String>, String> {
    let normalized = pattern.replace('\\', "/");
    // Walk from the longest leading part of the pattern that has no wildcards, and only as deep
    // as the pattern can reach.
    let parts: Vec<&str> = normalized.split('/').collect();
    let literal = parts.iter().take_while(|p| !has_glob_chars(p)).count();
    let base = if literal == 0 { ".".to_string() } else { parts[..literal].join("/") };
    let base = if base.is_empty() { "/".to_string() } else { base };
    if literal == parts.len() {
        return Ok(if Path::new(&normalized).exists() { vec![normalized] } else { Vec::new() });
    }
    let depth = if normalized.contains("**") { usize::MAX } else { parts.len() - literal };
    let matcher = globset::GlobBuilder::new(&normalized)
        .case_insensitive(cfg!(windows))
        .literal_separator(true)
        .build()
        .map_err(|e| format!("bad pattern '{pattern}': {e}"))?
        .compile_matcher();
    let mut out: Vec<String> = walkdir::WalkDir::new(&base)
        .max_depth(depth)
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_type().is_file())
        .map(|e| text(e.path()))
        .map(|p| p.strip_prefix("./").map(str::to_string).unwrap_or(p))
        .filter(|p| matcher.is_match(p))
        .collect();
    out.sort();
    Ok(out)
}

/// Sequence files in a folder (and its subfolders when `recursive`), sorted.
pub fn sequence_files(dir: &str, recursive: bool) -> Vec<String> {
    let mut out: Vec<String> = walkdir::WalkDir::new(dir)
        .max_depth(if recursive { usize::MAX } else { 1 })
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_type().is_file() && is_sequence_file(e.path()))
        .map(|e| text(e.path()))
        .collect();
    out.sort();
    out
}

/// Rasterise a scene to PNG bytes at `scale` pixels per scene unit. Text uses the system fonts.
pub fn scene_png(scene: &Scene, scale: f32) -> Result<Vec<u8>, String> {
    let svg = to_svg(scene);
    let mut options = resvg::usvg::Options::default();
    options.fontdb_mut().load_system_fonts();
    let tree = resvg::usvg::Tree::from_str(&svg, &options).map_err(|e| format!("could not rasterise the drawing: {e}"))?;
    let width = (scene.width * scale).ceil().max(1.0) as u32;
    let height = (scene.height * scale).ceil().max(1.0) as u32;
    let mut pixmap = resvg::tiny_skia::Pixmap::new(width, height).ok_or_else(|| format!("an image of {width}x{height} pixels is not possible"))?;
    resvg::render(&tree, resvg::tiny_skia::Transform::from_scale(scale, scale), &mut pixmap.as_mut());
    pixmap.encode_png().map_err(|e| format!("PNG encoding failed: {e}"))
}

/// Write a scene as SVG, or as PNG (rasterised at `scale`) when the path ends in `.png`.
pub fn save_scene(scene: &Scene, path: &Path, scale: f32) -> Result<(), String> {
    let is_png = path.extension().and_then(|e| e.to_str()).is_some_and(|e| e.eq_ignore_ascii_case("png"));
    let bytes = if is_png { scene_png(scene, scale)? } else { to_svg(scene).into_bytes() };
    slopgene_core::io::write_atomic(path, &bytes).map_err(|e| format!("{}: {e}", path.display()))
}

/// Encode an image of premultiplied RGBA pixels (row by row, 4 bytes each) as PNG.
pub fn rgba_png(width: u32, height: u32, rgba: Vec<u8>) -> Result<Vec<u8>, String> {
    let size = resvg::tiny_skia::IntSize::from_wh(width, height).ok_or("the image is empty")?;
    let pixmap = resvg::tiny_skia::Pixmap::from_vec(rgba, size).ok_or("the pixel data does not match the image size")?;
    pixmap.encode_png().map_err(|e| format!("PNG encoding failed: {e}"))
}

pub fn register(engine: &mut Engine) {
    engine.register_fn("glob", |pattern: &str| -> R<Array> { Ok(glob(pattern)?.into_iter().map(Dynamic::from).collect()) });
    engine.register_fn("files", |dir: &str| -> Array { sequence_files(dir, false).into_iter().map(Dynamic::from).collect() });
    engine.register_fn("files", |dir: &str, recursive: bool| -> Array { sequence_files(dir, recursive).into_iter().map(Dynamic::from).collect() });
    engine.register_fn("exists", |path: &str| Path::new(path).exists());
    engine.register_fn("file_name", |path: &str| Path::new(path).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default());
    engine.register_fn("file_stem", |path: &str| Path::new(path).file_stem().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default());
    engine.register_fn("read_text", |path: &str| -> R<String> { std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}").into()) });
    engine.register_fn("write_text", |path: &str, content: &str| -> R<()> {
        slopgene_core::io::write_atomic(Path::new(path), content.as_bytes()).map_err(|e| format!("{path}: {e}").into())
    });
    engine.register_fn("mkdir", |path: &str| -> R<()> { std::fs::create_dir_all(path).map_err(|e| format!("{path}: {e}").into()) });
}
