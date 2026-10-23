//! Render the app icon, `assets/icon/slopgene.svg`, into the files the installers and the app use:
//!
//! - `slopgene-<size>.png` for the Linux packages (16 to 512 pixels) and the window (256);
//! - `slopgene.ico` for Windows (the program file, the Start menu and the installer);
//! - `slopgene.icns` for macOS, with the artwork inside the margin macOS icons keep.
//!
//! Run `cargo run -p slopgene-gui --example make_icons` after changing the SVG, and commit the
//! files it writes.

use resvg::tiny_skia::{Pixmap, Transform};
use resvg::usvg;
use std::path::Path;

/// Where the icon source lives and the rendered files go.
const DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/assets/icon");

/// PNG sizes for the Linux icon theme folders (`hicolor/<size>x<size>/apps`).
const PNG_SIZES: [u32; 7] = [16, 32, 48, 64, 128, 256, 512];

/// Image sizes inside the Windows icon. Windows picks among them for each place and screen scale.
const ICO_SIZES: [u32; 6] = [16, 24, 32, 48, 64, 256];

/// The macOS icon types and their pixel sizes. Each holds a PNG. The `ic11` to `ic14` and
/// `ic10` types are the double-density versions of the next size down.
const ICNS_TYPES: [(&[u8; 4], u32); 10] = [
    (b"icp4", 16),
    (b"ic11", 32),
    (b"icp5", 32),
    (b"ic12", 64),
    (b"ic07", 128),
    (b"ic13", 256),
    (b"ic08", 256),
    (b"ic14", 512),
    (b"ic09", 512),
    (b"ic10", 1024),
];

/// The share of the canvas the plate fills in a macOS icon: Apple's grid puts an 824-pixel plate
/// on a 1024-pixel canvas. The SVG's own plate is 896 pixels wide, so it shrinks by 824 / 896.
const MACOS_SCALE: f32 = 824.0 / 896.0;

/// Render the icon at `size` pixels. A `scale` below 1 shrinks the artwork about the centre and
/// leaves a transparent margin around it.
fn render(tree: &usvg::Tree, size: u32, scale: f32) -> Pixmap {
    let mut pixmap = Pixmap::new(size, size).expect("icon sizes are not zero");
    let fit = size as f32 / tree.size().width() * scale;
    let margin = size as f32 * (1.0 - scale) / 2.0;
    resvg::render(tree, Transform::from_scale(fit, fit).post_translate(margin, margin), &mut pixmap.as_mut());
    pixmap
}

fn png(pixmap: &Pixmap) -> Vec<u8> {
    pixmap.encode_png().expect("a rendered icon encodes as PNG")
}

/// A 32-bit device-independent bitmap, the classic format of the images inside an `.ico`:
/// a header, the pixels as bottom-up BGRA rows, then a 1-bit mask. The alpha channel does the
/// masking, so the mask is all zeros. The header's height counts the pixel rows and the mask rows.
fn dib(pixmap: &Pixmap) -> Vec<u8> {
    let (width, height) = (pixmap.width(), pixmap.height());
    let mask_row = width.div_ceil(32) * 4;
    let mut out = Vec::new();
    out.extend(40u32.to_le_bytes());
    out.extend((width as i32).to_le_bytes());
    out.extend((2 * height as i32).to_le_bytes());
    out.extend(1u16.to_le_bytes());
    out.extend(32u16.to_le_bytes());
    out.extend(0u32.to_le_bytes());
    out.extend((width * height * 4 + mask_row * height).to_le_bytes());
    out.extend([0u8; 16]);
    // ii = row, counted from the bottom; iii = column.
    for ii in (0..height).rev() {
        for iii in 0..width {
            let color = pixmap.pixel(iii, ii).expect("the pixel lies inside the image").demultiply();
            out.extend([color.blue(), color.green(), color.red(), color.alpha()]);
        }
    }
    out.extend(vec![0u8; (mask_row * height) as usize]);
    out
}

/// An `.ico` file: a directory of images, then the images. Sizes up to 64 pixels are stored as
/// bitmaps, which every Windows program reads; 256 pixels is stored as PNG, as Windows expects.
fn ico(tree: &usvg::Tree) -> Vec<u8> {
    let images: Vec<(u32, Vec<u8>)> = ICO_SIZES
        .iter()
        .map(|&size| {
            let pixmap = render(tree, size, 1.0);
            (size, if size >= 256 { png(&pixmap) } else { dib(&pixmap) })
        })
        .collect();
    let mut out = Vec::new();
    out.extend(0u16.to_le_bytes());
    out.extend(1u16.to_le_bytes());
    out.extend((images.len() as u16).to_le_bytes());
    let mut offset = 6 + 16 * images.len() as u32;
    for (size, data) in &images {
        // A width and height of 0 stand for 256 pixels.
        let side = if *size >= 256 { 0 } else { *size as u8 };
        out.extend([side, side, 0, 0]);
        out.extend(1u16.to_le_bytes());
        out.extend(32u16.to_le_bytes());
        out.extend((data.len() as u32).to_le_bytes());
        out.extend(offset.to_le_bytes());
        offset += data.len() as u32;
    }
    for (_, data) in &images {
        out.extend(data);
    }
    out
}

/// An `.icns` file: the magic number and total length, then each image as type, length and PNG.
/// Lengths are big-endian and include their own 8-byte headers.
fn icns(tree: &usvg::Tree) -> Vec<u8> {
    let mut body = Vec::new();
    for (kind, size) in ICNS_TYPES {
        let data = png(&render(tree, size, MACOS_SCALE));
        body.extend(kind);
        body.extend((8 + data.len() as u32).to_be_bytes());
        body.extend(data);
    }
    let mut out = Vec::new();
    out.extend(b"icns");
    out.extend((8 + body.len() as u32).to_be_bytes());
    out.extend(body);
    out
}

fn write(name: &str, data: &[u8]) {
    let path = Path::new(DIR).join(name);
    std::fs::write(&path, data).unwrap_or_else(|e| panic!("could not write {}: {e}", path.display()));
    println!("wrote {} ({} bytes)", path.display(), data.len());
}

fn main() {
    let source = Path::new(DIR).join("slopgene.svg");
    let svg = std::fs::read_to_string(&source).unwrap_or_else(|e| panic!("could not read {}: {e}", source.display()));
    let tree = usvg::Tree::from_str(&svg, &usvg::Options::default()).expect("the icon is valid SVG");
    for size in PNG_SIZES {
        write(&format!("slopgene-{size}.png"), &png(&render(&tree, size, 1.0)));
    }
    write("slopgene.ico", &ico(&tree));
    write("slopgene.icns", &icns(&tree));
}
