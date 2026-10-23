//! SlopGene desktop app.
//!
//! Usage: `slopgene-gui [FILES...] [--dark] [--view map|sequence|split|alignment] [--select A..B]
//! [--feature NAME] [--palette] [--console] [--digest ENZYMES] [--align READ] [--step N]
//! [--zoom] [--new-protocol [--enzymes NAMES]] [--size WxH] [--screenshot OUT.png [--after FRAMES] [--drag X1,Y1,X2,Y2]]`
//!
//! A file may be a protocol (`slopgene.protocol/1` JSON), which opens in a protocol tab.
//!
//! `--align READ` aligns READ to the first file and opens the alignment view.
//!
//! `--step N` shows step N (from 1) of the protocol that opens, and `--zoom` zooms into its
//! card. `--new-protocol` starts a cloning protocol from the files open, and `--enzymes` cuts its
//! insert out with those enzymes. The region marked with `--select` or `--feature` on the file
//! in front is the part it cuts out.
//!
//! `--screenshot` renders the window, saves it as PNG and quits; it is used to check the
//! interface without a person at the screen. `--drag` plays a held drag into it first, from one
//! point to another (in points from the window's top left), so the picture shows a drag on its way.

// No console window behind the app in release builds on Windows.
#![cfg_attr(all(windows, not(debug_assertions)), windows_subsystem = "windows")]

mod align_view;
mod app;
mod doc;
mod map_view;
mod panels;
mod protocol_view;
mod scene_paint;
mod seq_view;
mod style;

use app::{SlopGeneApp, StartOptions};
use doc::View;
use eframe::egui;
use std::path::PathBuf;

fn parse_args() -> (StartOptions, [f32; 2]) {
    let mut options = StartOptions { screenshot_after: 12, ..StartOptions::default() };
    let mut size = [1360.0, 860.0];
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--dark" => options.dark = true,
            "--palette" => options.palette = true,
            "--console" => options.console = true,
            "--view" => {
                options.view = match args.next().as_deref() {
                    Some("map") => Some(View::Map),
                    Some("sequence") | Some("seq") => Some(View::Sequence),
                    Some("split") => Some(View::Split),
                    Some("alignment") | Some("align") => Some(View::Alignment),
                    _ => None,
                }
            }
            "--select" => options.select = args.next(),
            "--feature" => options.feature = args.next(),
            "--digest" => options.digest = args.next(),
            "--align" => options.align = args.next().map(PathBuf::from),
            "--step" => options.step = args.next().and_then(|n| n.parse().ok()),
            "--zoom" => options.zoom = true,
            "--drag" => {
                // x1,y1,x2,y2 in points, from the top left of the window.
                let numbers: Vec<f32> = args.next().unwrap_or_default().split(',').filter_map(|n| n.trim().parse().ok()).collect();
                if let [x1, y1, x2, y2] = numbers[..] {
                    options.drag = Some((egui::pos2(x1, y1), egui::pos2(x2, y2)));
                }
            }
            "--new-protocol" => options.new_protocol = true,
            "--enzymes" => options.enzymes = args.next(),
            "--screenshot" => options.screenshot = args.next().map(PathBuf::from),
            "--after" => options.screenshot_after = args.next().and_then(|n| n.parse().ok()).unwrap_or(12),
            "--size" => {
                if let Some((w, h)) = args.next().as_deref().and_then(|s| s.split_once('x')) {
                    size = [w.parse().unwrap_or(size[0]), h.parse().unwrap_or(size[1])];
                }
            }
            "-h" | "--help" => {
                println!(
                    "slopgene-gui [FILES...] [--dark] [--view map|sequence|split|alignment] [--select A..B] [--feature NAME] [--palette] [--console] [--digest ENZYMES] [--align READ] [--step N] [--zoom] [--new-protocol [--enzymes NAMES]] [--size WxH] [--screenshot OUT.png [--after N] [--drag X1,Y1,X2,Y2]]"
                );
                std::process::exit(0);
            }
            other => options.files.push(PathBuf::from(other)),
        }
    }
    (options, size)
}

/// Where the app keeps its settings (recent files, dark mode, window size). The app was called
/// FukGene before, so on the first start under the new name the old settings are copied over.
fn settings_file() -> Option<PathBuf> {
    let file = eframe::storage_dir("SlopGene")?.join("app.ron");
    if !file.exists()
        && let Some(old) = eframe::storage_dir("FukGene").map(|dir| dir.join("app.ron"))
        && old.exists()
        && let Some(folder) = file.parent()
    {
        let _ = std::fs::create_dir_all(folder).and_then(|_| std::fs::copy(&old, &file));
    }
    Some(file)
}

fn main() -> eframe::Result<()> {
    let (start, size) = parse_args();
    // Linux desktops match a window to its launcher, and so to the launcher's icon, by this id.
    // The installers name the launcher slopgene-gui.desktop.
    let mut viewport = egui::ViewportBuilder::default()
        .with_title("SlopGene")
        .with_app_id("slopgene-gui")
        .with_inner_size(size)
        .with_min_inner_size([760.0, 480.0])
        .with_drag_and_drop(true);
    // The icon in the title bar, the taskbar and the window switcher. On macOS the app bundle
    // supplies the icon, and setting one here would replace it with this full-bleed square.
    if cfg!(not(target_os = "macos")) {
        let icon = eframe::icon_data::from_png_bytes(include_bytes!("../assets/icon/slopgene-256.png")).expect("the bundled icon is a PNG");
        viewport = viewport.with_icon(icon);
    }
    let native = eframe::NativeOptions {
        viewport,
        // Screenshot runs should not reuse or overwrite the saved window state.
        persist_window: start.screenshot.is_none(),
        // eframe would keep the settings in a folder named after the app id. They stay in the
        // folder named after the app, where they were before the app id was set.
        persistence_path: settings_file(),
        ..eframe::NativeOptions::default()
    };
    eframe::run_native("SlopGene", native, Box::new(move |cc| Ok(Box::new(SlopGeneApp::new(cc, start)))))
}
