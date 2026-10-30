//! The inspector on the right: the selection, the selected feature's editor, and lists of
//! features, enzymes and primers.

use crate::app::{Action, ViewSettings};
use crate::doc::DocTab;
use crate::style::{self, Palette};
use eframe::egui::{self, Color32, CornerRadius, RichText, Sense, vec2};
use slopgene_core::primers::{TmConditions, analyze};
use slopgene_core::render::{feature_rgba, thousands};
use slopgene_core::translate::{GeneticCode, translate, translate_feature};
use slopgene_core::*;

/// Feature types offered in the editor (any other type typed in a file is kept as is).
pub const FEATURE_TYPES: &[&str] = &[
    "CDS",
    "gene",
    "promoter",
    "terminator",
    "rep_origin",
    "primer_bind",
    "protein_bind",
    "misc_feature",
    "enhancer",
    "polyA_signal",
    "LTR",
    "RBS",
    "sig_peptide",
    "regulatory",
    "repeat_region",
    "misc_RNA",
    "intron",
    "exon",
    "5'UTR",
    "3'UTR",
];

pub(crate) const ENZYME_SETS: &[&str] = &["common", "6+", "8+", "neb", "thermo", "golden-gate", "type-iis", "blunt", "all"];

pub(crate) fn section_title(ui: &mut egui::Ui, text: &str, palette: &Palette) {
    ui.add_space(4.0);
    ui.label(RichText::new(text.to_uppercase()).font(style::bold_font(10.5)).color(palette.muted));
    ui.add_space(2.0);
}

fn key_value(ui: &mut egui::Ui, key: &str, value: impl Into<String>, palette: &Palette) {
    ui.horizontal(|ui| {
        ui.label(RichText::new(key).color(palette.muted));
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label(value.into());
        });
    });
}

/// A small rounded colour swatch.
pub(crate) fn swatch(ui: &mut egui::Ui, color: Color32) {
    let (rect, _) = ui.allocate_exact_size(vec2(10.0, 10.0), Sense::hover());
    ui.painter().rect_filled(rect, CornerRadius::same(3), color);
}

pub fn inspector(ui: &mut egui::Ui, tab: &mut DocTab, settings: &ViewSettings, actions: &mut Vec<Action>) {
    let palette = Palette::get(settings.dark);
    egui::ScrollArea::vertical().id_salt(("inspector", tab.id)).auto_shrink([false, false]).show(ui, |ui| {
        ui.set_width(ui.available_width());
        selection_section(ui, tab, &palette, actions);
        if let Some(uid) = tab.selected_feature {
            ui.separator();
            feature_editor(ui, tab, uid, &palette, actions);
        }
        ui.separator();
        features_section(ui, tab, settings, &palette);
        ui.separator();
        enzymes_section(ui, tab, &palette);
        if !tab.doc.primers.is_empty() {
            ui.separator();
            primers_section(ui, tab, &palette);
        }
        ui.add_space(12.0);
    });
}

fn selection_section(ui: &mut egui::Ui, tab: &mut DocTab, palette: &Palette, actions: &mut Vec<Action>) {
    let len = tab.len();
    section_title(ui, "Selection", palette);
    let Some(sel) = tab.selection else {
        key_value(ui, "Cursor", format!("{} of {}", thousands(tab.cursor.min(len)), thousands(len)), palette);
        ui.label(RichText::new("Drag across the map or the sequence to select.").color(palette.faint).small());
        return;
    };
    let bases = sel.extract(&tab.doc.seq);
    key_value(ui, "Range", sel.format(len), palette);
    key_value(ui, "Length", format!("{} bp", thousands(bases.len())), palette);
    key_value(ui, "GC", format!("{:.1}%", alphabet::gc_fraction(&bases) * 100.0), palette);
    if (8..=80).contains(&bases.len()) {
        let report = analyze(&bases, &TmConditions::default());
        key_value(ui, "Tm", format!("{:.1} °C", report.tm), palette);
    }
    if bases.len() >= 3 {
        let protein = translate(&bases, GeneticCode::standard(), false);
        let shown: String = protein.iter().take(120).map(|&b| b as char).collect();
        ui.label(RichText::new("Translation (frame 1)").color(palette.muted).small());
        ui.label(RichText::new(if protein.len() > 120 { format!("{shown}…") } else { shown }).font(style::mono_font(11.5)));
    }
    ui.add_space(4.0);
    ui.horizontal_wrapped(|ui| {
        if ui
            .button("Clone")
            .on_hover_text("Start a cloning protocol that cuts these bases out. The enzymes that cut nearest them are listed to choose from.")
            .clicked()
        {
            actions.push(Action::NewProtocol);
        }
        if ui.button("Add feature").clicked() {
            actions.push(Action::AddFeatureFromSelection);
        }
        if ui.button("Copy").clicked() {
            actions.push(Action::CopySelection);
        }
        if ui.button("Copy rev. comp.").clicked() {
            actions.push(Action::CopyReverseComplement);
        }
        if (12..=80).contains(&bases.len()) && ui.button("Save as primer").clicked() {
            actions.push(Action::PrimerFromSelection);
        }
    });
}

fn feature_editor(ui: &mut egui::Ui, tab: &mut DocTab, uid: u64, palette: &Palette, actions: &mut Vec<Action>) {
    let Some(index) = tab.doc.feature_index(uid) else { return };
    let len = tab.len();
    section_title(ui, "Feature", palette);
    let feature = tab.doc.features[index].clone();

    // Name: edited in a buffer and applied as one undo step.
    let name_id = egui::Id::unique(("feature-name", uid));
    let mut name: String = ui.ctx().data(|d| d.get_temp(name_id)).unwrap_or_else(|| feature.name.clone());
    let response = ui.add(egui::TextEdit::singleline(&mut name).desired_width(f32::INFINITY).hint_text("name"));
    if response.changed() {
        ui.ctx().data_mut(|d| d.insert_temp(name_id, name.clone()));
        let value = name.clone();
        let _ = tab.edit_coalesced(&format!("name-{uid}"), |d| {
            if let Some(i) = d.feature_index(uid) {
                d.features[i].name = value;
            }
            Ok(())
        });
    }
    if response.lost_focus() {
        ui.ctx().data_mut(|d| d.remove::<String>(name_id));
    }

    ui.horizontal(|ui| {
        let mut kind = feature.kind.clone();
        egui::ComboBox::from_id_salt(("feature-type", uid)).selected_text(&kind).width(120.0).show_ui(ui, |ui| {
            for t in FEATURE_TYPES {
                ui.selectable_value(&mut kind, t.to_string(), *t);
            }
        });
        if kind != feature.kind {
            let _ = tab.edit(|d| {
                if let Some(i) = d.feature_index(uid) {
                    d.features[i].kind = kind;
                }
                Ok(())
            });
        }
        let mut strand = feature.strand;
        egui::ComboBox::from_id_salt(("feature-strand", uid)).selected_text(strand.name()).width(84.0).show_ui(ui, |ui| {
            for s in [Strand::Forward, Strand::Reverse, Strand::None, Strand::Both] {
                ui.selectable_value(&mut strand, s, s.name());
            }
        });
        if strand != feature.strand {
            let _ = tab.edit(|d| {
                if let Some(i) = d.feature_index(uid) {
                    d.features[i].strand = strand;
                }
                Ok(())
            });
        }
        let current = feature_rgba(&feature, slopgene_core::render::ColorMode::File);
        let mut rgb = [current.r, current.g, current.b];
        if ui.color_edit_button_srgb(&mut rgb).changed() {
            let _ = tab.edit_coalesced(&format!("color-{uid}"), |d| {
                if let Some(i) = d.feature_index(uid) {
                    d.features[i].set_color(Some(Color::rgb(rgb[0], rgb[1], rgb[2])));
                }
                Ok(())
            });
        }
    });
    let span = feature.span(len);
    key_value(ui, "Range", span.format(len), palette);
    key_value(ui, "Length", format!("{} bp", thousands(feature.len())), palette);
    if feature.segments.len() > 1 {
        key_value(ui, "Segments", feature.segments.len().to_string(), palette);
    }
    if feature.is_coding() {
        let protein = translate_feature(&feature, &tab.doc);
        let aa: Vec<u8> = protein.iter().copied().take_while(|&a| a != b'*').collect();
        key_value(ui, "Protein", format!("{} aa, {:.1} kDa", aa.len(), stats::protein_mass(&aa) / 1000.0), palette);
        let text: String = protein.iter().map(|&b| b as char).collect();
        ui.add(egui::Label::new(RichText::new(text).font(style::mono_font(11.0)).color(palette.muted)).wrap());
    }
    if let Some(note) = feature.qualifier_text("note").filter(|n| !n.is_empty()) {
        let plain = slopgene_core::io::snapgene::html_to_text(&note);
        ui.add(egui::Label::new(RichText::new(plain).small().color(palette.muted)).wrap());
    }
    ui.horizontal(|ui| {
        if ui.button("Copy sequence").clicked() {
            ui.ctx().copy_text(String::from_utf8_lossy(&feature.extract(&tab.doc.seq)).into_owned());
        }
        if ui.button("Delete").clicked() {
            actions.push(Action::DeleteSelectedFeature);
        }
    });
}

fn features_section(ui: &mut egui::Ui, tab: &mut DocTab, settings: &ViewSettings, palette: &Palette) {
    let len = tab.len();
    section_title(ui, &format!("Features · {}", tab.doc.features.len()), palette);
    let filter_id = egui::Id::unique(("feature-filter", tab.id));
    let mut filter: String = ui.ctx().data(|d| d.get_temp(filter_id)).unwrap_or_default();
    if tab.doc.features.len() > 8 && ui.add(egui::TextEdit::singleline(&mut filter).hint_text("filter").desired_width(f32::INFINITY)).changed() {
        ui.ctx().data_mut(|d| d.insert_temp(filter_id, filter.clone()));
    }
    let needle = filter.to_lowercase();
    let mut clicked = None;
    for feature in &tab.doc.features {
        if !needle.is_empty() && !feature.name.to_lowercase().contains(&needle) && !feature.kind.to_lowercase().contains(&needle) {
            continue;
        }
        let selected = tab.selected_feature == Some(feature.uid);
        let rgba = feature_rgba(feature, settings.colors);
        let response = ui
            .horizontal(|ui| {
                swatch(ui, Color32::from_rgb(rgba.r, rgba.g, rgba.b));
                let name = RichText::new(&feature.name).color(if selected { palette.accent } else { palette.ink });
                ui.add(egui::Label::new(name).truncate().sense(Sense::click()))
                    | ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        ui.label(RichText::new(feature.span(len).format(len)).small().color(palette.faint))
                    })
                    .inner
            })
            .inner;
        if response.clicked() {
            clicked = Some(feature.uid);
        }
    }
    if let Some(uid) = clicked {
        tab.select_feature(uid);
    }
}

fn enzymes_section(ui: &mut egui::Ui, tab: &mut DocTab, palette: &Palette) {
    section_title(ui, "Enzymes", palette);
    ui.horizontal(|ui| {
        let mut set = tab.enzyme_set.clone();
        egui::ComboBox::from_id_salt(("enzyme-set", tab.id)).selected_text(&set).width(110.0).show_ui(ui, |ui| {
            for s in ENZYME_SETS {
                ui.selectable_value(&mut set, s.to_string(), *s);
            }
        });
        if set != tab.enzyme_set {
            tab.enzyme_set = set;
        }
        ui.checkbox(&mut tab.unique_only, "unique only");
    });
    let unique_only = tab.unique_only;
    let counts: Vec<(String, String, Vec<usize>, Vec<usize>)> =
        tab.sites().1.iter().map(|c| (c.name.clone(), c.site.clone(), c.cuts.clone(), c.blocked.clone())).collect();
    let methylation = tab.doc.methylation.describe();
    let mut jump = None;
    let shown: Vec<&(String, String, Vec<usize>, Vec<usize>)> = counts.iter().filter(|c| !unique_only || c.2.len() == 1).collect();
    if shown.is_empty() {
        ui.label(RichText::new("No cuts with this set.").small().color(palette.faint));
    }
    for (name, site, cuts, blocked) in shown.iter().take(150) {
        let response = ui
            .horizontal(|ui| {
                let label = ui.add(egui::Label::new(RichText::new(name)).sense(Sense::click()));
                ui.label(RichText::new(site).small().font(style::mono_font(10.5)).color(palette.faint));
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    // The cut, or how many cuts, and how many sites methylation blocks.
                    let mut text = match cuts.len() {
                        0 => String::new(),
                        1 => cuts[0].to_string(),
                        n => format!("{n}×"),
                    };
                    if !blocked.is_empty() {
                        text = format!("{text}  {} blocked", blocked.len()).trim().to_string();
                    }
                    let count = ui.label(RichText::new(text).small().color(palette.muted));
                    if !blocked.is_empty() {
                        count.on_hover_text(format!("{} site(s) that methylation blocks ({methylation}), so they are not cut", blocked.len()));
                    }
                });
                label
            })
            .inner;
        if response.clicked() {
            jump = cuts.first().or(blocked.first()).copied();
        }
    }
    if let Some(pos) = jump {
        tab.select(None);
        tab.cursor = pos.min(tab.len());
        tab.scroll_to = Some(tab.cursor);
    }
}

fn primers_section(ui: &mut egui::Ui, tab: &mut DocTab, palette: &Palette) {
    section_title(ui, &format!("Primers · {}", tab.doc.primers.len()), palette);
    let len = tab.len();
    let mut select = None;
    for primer in &tab.doc.primers {
        let sites: Vec<&BindingSite> = primer.sites.iter().filter(|s| s.extra.attr("simplified") != Some("1")).collect();
        let response = ui
            .horizontal(|ui| {
                let label = ui.add(egui::Label::new(RichText::new(&primer.name)).truncate().sense(Sense::click()));
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let where_ = match sites.first() {
                        Some(s) => format!("{} {}", s.span.format(len), if s.strand.is_reverse() { "←" } else { "→" }),
                        None => "no site".into(),
                    };
                    ui.label(RichText::new(where_).small().color(palette.faint));
                });
                label
            })
            .inner;
        if response.clicked() {
            select = sites.first().map(|s| s.span);
        }
    }
    if let Some(span) = select {
        tab.select(Some(span));
        tab.scroll_to = Some(span.start);
    }
}
