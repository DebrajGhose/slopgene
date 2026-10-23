//! Look and feel: fonts, colours and spacing.
//!
//! The aim is a quiet interface: one neutral background, hairline separators, a single accent
//! colour (indigo from the user's figure palette), Inter for the interface and JetBrains Mono for
//! sequence. Both fonts are bundled (SIL Open Font License, see `assets/fonts`), so the app looks
//! the same on every system.

use eframe::egui::epaint::text::VariationCoords;
use eframe::egui::{self, Color32, CornerRadius, FontData, FontDefinitions, FontFamily, FontId, FontTweak, Shadow, Stroke, TextStyle};
use std::sync::Arc;

static INTER: &[u8] = include_bytes!("../assets/fonts/Inter.ttf");
static JETBRAINS_MONO: &[u8] = include_bytes!("../assets/fonts/JetBrainsMono.ttf");

/// Font family for semibold interface text.
pub fn bold_family() -> FontFamily {
    FontFamily::Name("bold".into())
}

/// Font family for semibold sequence text.
pub fn mono_bold_family() -> FontFamily {
    FontFamily::Name("mono-bold".into())
}

pub fn ui_font(size: f32) -> FontId {
    FontId::new(size, FontFamily::Proportional)
}

pub fn bold_font(size: f32) -> FontId {
    FontId::new(size, bold_family())
}

pub fn mono_font(size: f32) -> FontId {
    FontId::new(size, FontFamily::Monospace)
}

/// A variable font fixed at the given weight (and optical size, where the font has one).
fn instance(bytes: &'static [u8], weight: f32, optical_size: Option<f32>) -> Arc<FontData> {
    let mut coords = VariationCoords::new([(b"wght", weight)]);
    if let Some(size) = optical_size {
        coords.push(b"opsz", size);
    }
    Arc::new(FontData::from_static(bytes).tweak(FontTweak { coords, ..FontTweak::default() }))
}

pub fn install_fonts(ctx: &egui::Context) {
    let mut fonts = FontDefinitions::default();
    fonts.font_data.insert("inter".into(), instance(INTER, 430.0, Some(14.0)));
    fonts.font_data.insert("inter-semibold".into(), instance(INTER, 600.0, Some(14.0)));
    fonts.font_data.insert("jbmono".into(), instance(JETBRAINS_MONO, 420.0, None));
    fonts.font_data.insert("jbmono-bold".into(), instance(JETBRAINS_MONO, 650.0, None));

    // Put the bundled fonts first and keep egui's own fonts behind them for symbols and emoji.
    let fallback_proportional = fonts.families.get(&FontFamily::Proportional).cloned().unwrap_or_default();
    let fallback_mono = fonts.families.get(&FontFamily::Monospace).cloned().unwrap_or_default();
    let with = |first: &str, rest: &[String]| -> Vec<String> {
        let mut v = vec![first.to_string()];
        v.extend(rest.iter().cloned());
        v
    };
    fonts.families.insert(FontFamily::Proportional, with("inter", &fallback_proportional));
    fonts.families.insert(FontFamily::Monospace, with("jbmono", &fallback_mono));
    fonts.families.insert(bold_family(), with("inter-semibold", &fallback_proportional));
    fonts.families.insert(mono_bold_family(), with("jbmono-bold", &fallback_mono));
    ctx.set_fonts(fonts);
}

/// The handful of colours the interface uses, for one appearance.
#[derive(Clone, Copy, Debug)]
pub struct Palette {
    pub dark: bool,
    /// Main background (canvas behind maps and sequence).
    pub canvas: Color32,
    /// Side panels and bars.
    pub panel: Color32,
    /// Raised surfaces: popups, text fields.
    pub raised: Color32,
    /// Hairline separators.
    pub line: Color32,
    pub ink: Color32,
    pub muted: Color32,
    pub faint: Color32,
    pub accent: Color32,
    /// Accent at low opacity, for selections.
    pub accent_soft: Color32,
    /// Caution: a step that waits for a choice, a file changed on disk.
    pub warn: Color32,
    /// A step that finished.
    pub ok: Color32,
    /// A step that failed.
    pub error: Color32,
}

impl Palette {
    pub fn light() -> Palette {
        Palette {
            dark: false,
            canvas: Color32::from_rgb(0xFC, 0xFC, 0xFB),
            panel: Color32::from_rgb(0xF6, 0xF6, 0xF4),
            raised: Color32::from_rgb(0xFF, 0xFF, 0xFF),
            line: Color32::from_rgb(0xE4, 0xE4, 0xE0),
            ink: Color32::from_rgb(0x22, 0x22, 0x22),
            muted: Color32::from_rgb(0x6B, 0x6B, 0x6B),
            faint: Color32::from_rgb(0xA8, 0xA8, 0xA4),
            accent: Color32::from_rgb(0x51, 0x5F, 0xAA),
            accent_soft: Color32::from_rgba_unmultiplied(0x51, 0x5F, 0xAA, 46),
            warn: Color32::from_rgb(0xB0, 0x6A, 0x00),
            ok: Color32::from_rgb(0x2F, 0x80, 0x55),
            error: Color32::from_rgb(0xB4, 0x3C, 0x36),
        }
    }

    pub fn dark() -> Palette {
        Palette {
            dark: true,
            canvas: Color32::from_rgb(0x1B, 0x1C, 0x1F),
            panel: Color32::from_rgb(0x22, 0x23, 0x27),
            raised: Color32::from_rgb(0x2A, 0x2B, 0x30),
            line: Color32::from_rgb(0x34, 0x35, 0x3A),
            ink: Color32::from_rgb(0xE8, 0xE8, 0xE6),
            muted: Color32::from_rgb(0x9A, 0x9A, 0x9A),
            faint: Color32::from_rgb(0x66, 0x67, 0x6C),
            accent: Color32::from_rgb(0x8A, 0x96, 0xE0),
            accent_soft: Color32::from_rgba_unmultiplied(0x8A, 0x96, 0xE0, 52),
            warn: Color32::from_rgb(0xE0, 0xA0, 0x40),
            ok: Color32::from_rgb(0x6C, 0xC0, 0x8E),
            error: Color32::from_rgb(0xE8, 0x7C, 0x72),
        }
    }

    pub fn get(dark: bool) -> Palette {
        if dark { Palette::dark() } else { Palette::light() }
    }
}

/// Apply the palette to egui's visuals and set sizes and spacing.
pub fn apply(ctx: &egui::Context, dark: bool) {
    let p = Palette::get(dark);
    let mut v = if dark { egui::Visuals::dark() } else { egui::Visuals::light() };
    v.panel_fill = p.panel;
    v.window_fill = p.raised;
    v.extreme_bg_color = p.raised;
    v.faint_bg_color = if dark { Color32::from_rgb(0x26, 0x27, 0x2B) } else { Color32::from_rgb(0xF0, 0xF0, 0xED) };
    v.code_bg_color = v.faint_bg_color;
    // Text fields need a fill that differs from both panels and windows to stay visible.
    v.text_edit_bg_color = Some(if dark { Color32::from_rgb(0x1E, 0x1F, 0x23) } else { Color32::from_rgb(0xF1, 0xF1, 0xEE) });
    v.window_stroke = Stroke::new(1.0, p.line);
    v.window_corner_radius = CornerRadius::same(10);
    v.menu_corner_radius = CornerRadius::same(8);
    v.window_shadow = Shadow { offset: [0, 6], blur: 24, spread: 0, color: Color32::from_black_alpha(if dark { 90 } else { 30 }) };
    v.popup_shadow = Shadow { offset: [0, 4], blur: 14, spread: 0, color: Color32::from_black_alpha(if dark { 80 } else { 26 }) };
    v.selection.bg_fill = p.accent_soft;
    v.selection.stroke = Stroke::new(1.0, p.accent);
    v.hyperlink_color = p.accent;
    v.warn_fg_color = p.warn;
    v.striped = false;

    let radius = CornerRadius::same(6);
    let w = &mut v.widgets;
    w.noninteractive.bg_fill = p.panel;
    w.noninteractive.weak_bg_fill = p.panel;
    w.noninteractive.bg_stroke = Stroke::new(1.0, p.line);
    w.noninteractive.fg_stroke = Stroke::new(1.0, p.ink);
    w.noninteractive.corner_radius = radius;
    let button = if dark { Color32::from_rgb(0x2E, 0x2F, 0x34) } else { Color32::from_rgb(0xEC, 0xEC, 0xE9) };
    let button_hover = if dark { Color32::from_rgb(0x38, 0x39, 0x3F) } else { Color32::from_rgb(0xE3, 0xE3, 0xDF) };
    let button_active = if dark { Color32::from_rgb(0x42, 0x43, 0x4A) } else { Color32::from_rgb(0xD9, 0xD9, 0xD5) };
    for (state, fill) in [(&mut w.inactive, button), (&mut w.hovered, button_hover), (&mut w.active, button_active), (&mut w.open, button_hover)] {
        state.bg_fill = fill;
        state.weak_bg_fill = fill;
        state.bg_stroke = Stroke::NONE;
        state.fg_stroke = Stroke::new(1.0, p.ink);
        state.corner_radius = radius;
        state.expansion = 0.0;
    }
    w.hovered.bg_stroke = Stroke::new(1.0, p.line);
    ctx.set_visuals(v);

    ctx.all_styles_mut(|style| {
        style.text_styles = [
            (TextStyle::Small, ui_font(11.0)),
            (TextStyle::Body, ui_font(13.5)),
            (TextStyle::Button, ui_font(13.5)),
            (TextStyle::Heading, bold_font(17.0)),
            (TextStyle::Monospace, mono_font(13.0)),
        ]
        .into();
        let s = &mut style.spacing;
        s.item_spacing = egui::vec2(8.0, 6.0);
        s.button_padding = egui::vec2(10.0, 4.0);
        s.interact_size.y = 26.0;
        s.window_margin = egui::Margin::same(14);
        s.menu_margin = egui::Margin::same(6);
        s.indent = 14.0;
        s.scroll.bar_width = 6.0;
        s.scroll.floating = true;
    });
}
