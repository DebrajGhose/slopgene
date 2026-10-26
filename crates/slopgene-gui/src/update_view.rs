//! Updating the app from GitHub: the check for a newer version, the dialog that says what the
//! update would do and shows its progress, and the install itself, through `slopgene_update`.
//! The check and the download run on threads, and the dialog reads where they have got to each
//! frame.

use crate::style::{self, Palette};
use eframe::egui::{self, RichText};
use slopgene_update::{Asset, Done, Error, Installation, Os, Release, Version};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Duration;

/// Where the update has got to.
#[derive(Clone, Debug)]
pub enum Stage {
    /// Nothing asked yet.
    Idle,
    Checking,
    UpToDate,
    /// A newer release, with the file that updates this installation when the release has one.
    Available {
        release: Release,
        asset: Option<Asset>,
        installation: Installation,
    },
    Downloading {
        release: Release,
        done: u64,
        total: u64,
    },
    Installing {
        release: Release,
    },
    Finished {
        release: Release,
        done: Done,
    },
    Failed {
        error: String,
    },
}

/// The update check, the download and the dialog.
pub struct Updater {
    stage: Arc<Mutex<Stage>>,
    /// Set to stop a download.
    stop: Arc<AtomicBool>,
    /// Whether the dialog is open.
    pub open: bool,
    /// Whether the app looks for a newer version when it starts.
    pub at_start: bool,
    /// Set once the app should close so the new version can take its place.
    pub quit: bool,
}

impl Updater {
    pub fn new(at_start: bool) -> Updater {
        Updater { stage: Arc::new(Mutex::new(Stage::Idle)), stop: Arc::new(AtomicBool::new(false)), open: false, at_start, quit: false }
    }

    fn stage(&self) -> Stage {
        self.stage.lock().map(|s| s.clone()).unwrap_or(Stage::Idle)
    }

    fn set(stage: &Arc<Mutex<Stage>>, value: Stage) {
        if let Ok(mut current) = stage.lock() {
            *current = value;
        }
    }

    /// Ask GitHub for the newest release, on a thread.
    pub fn check(&mut self) {
        if matches!(self.stage(), Stage::Checking | Stage::Downloading { .. } | Stage::Installing { .. }) {
            return;
        }
        Self::set(&self.stage, Stage::Checking);
        let stage = self.stage.clone();
        std::thread::spawn(move || {
            let found = match slopgene_update::latest() {
                Ok(release) if release.version > Version::current() => {
                    let installation = Installation::detect();
                    let asset = slopgene_update::choose_asset(&release, &installation, Os::this()).cloned();
                    Stage::Available { release, asset, installation }
                }
                Ok(_) | Err(Error::NoRelease) => Stage::UpToDate,
                Err(e) => Stage::Failed { error: e.to_string() },
            };
            Self::set(&stage, found);
        });
    }

    /// The newer version found, for the button in the top bar.
    pub fn available(&self) -> Option<Version> {
        match self.stage() {
            Stage::Available { release, .. } | Stage::Downloading { release, .. } | Stage::Installing { release } | Stage::Finished { release, .. } => {
                Some(release.version)
            }
            _ => None,
        }
    }

    /// Open the dialog, looking for a newer version first when none is known yet.
    pub fn open_dialog(&mut self) {
        self.open = true;
        if matches!(self.stage(), Stage::Idle | Stage::UpToDate | Stage::Failed { .. }) {
            self.check();
        }
    }

    /// Download the release found and install it, on a thread.
    fn install(&mut self) {
        let Stage::Available { release, asset: Some(asset), installation } = self.stage() else { return };
        self.stop.store(false, Ordering::SeqCst);
        Self::set(&self.stage, Stage::Downloading { release: release.clone(), done: 0, total: asset.size });
        let stage = self.stage.clone();
        let stop = self.stop.clone();
        std::thread::spawn(move || {
            let told = stage.clone();
            let shown = release.clone();
            let result = slopgene_update::download(&release, &asset, &mut |done, total| {
                Self::set(&told, Stage::Downloading { release: shown.clone(), done, total });
                !stop.load(Ordering::SeqCst)
            })
            .and_then(|download| {
                Self::set(&stage, Stage::Installing { release: release.clone() });
                slopgene_update::install(&installation, &download, true)
            });
            let outcome = match result {
                Ok(done) => Stage::Finished { release, done },
                Err(Error::Cancelled) => Stage::Available { release, asset: Some(asset), installation },
                Err(e) => Stage::Failed { error: e.to_string() },
            };
            Self::set(&stage, outcome);
        });
    }

    /// Draw the dialog when it is open and act on its buttons. `unsaved` names the documents with
    /// unsaved changes, which the install would close.
    pub fn show(&mut self, ctx: &egui::Context, palette: &Palette, unsaved: &[String]) {
        let stage = self.stage();
        if matches!(stage, Stage::Checking | Stage::Downloading { .. } | Stage::Installing { .. }) {
            ctx.request_repaint_after(Duration::from_millis(100));
        }
        // Once the new version is on its way, the app closes for it.
        if let Stage::Finished { done, .. } = &stage
            && !matches!(done, Done::RunByHand { .. })
        {
            self.quit = true;
        }
        if !self.open {
            return;
        }
        let current = Version::current();
        let (mut close, mut install, mut again, mut stop) = (false, false, false, false);
        let modal = egui::Modal::new(egui::Id::unique("update")).show(ctx, |ui| {
            ui.set_width(460.0);
            ui.label(RichText::new("Update SlopGene").font(style::bold_font(16.0)).color(palette.ink));
            ui.add_space(8.0);
            match &stage {
                Stage::Idle | Stage::Checking => {
                    ui.horizontal(|ui| {
                        ui.spinner();
                        ui.label("Looking for a newer version…");
                    });
                }
                Stage::UpToDate => {
                    ui.label(format!("You have the newest version, {current}."));
                }
                Stage::Available { release, asset, installation } => {
                    ui.label(RichText::new(format!("SlopGene {} is available. You have {current}.", release.version)).color(palette.ink));
                    if !release.published.is_empty() {
                        ui.label(RichText::new(format!("Released on {}.", release.published)).small().color(palette.muted));
                    }
                    ui.add_space(6.0);
                    match asset {
                        Some(asset) => {
                            ui.label(format!(
                                "Download and install fetches {} ({}) and checks it against the release's checksums.",
                                asset.name,
                                megabytes(asset.size)
                            ));
                            ui.label(installation.describe());
                            if !unsaved.is_empty() {
                                ui.add_space(4.0);
                                ui.label(RichText::new(format!("Save your work first: {} unsaved.", unsaved.join(", "))).color(palette.warn));
                            }
                        }
                        None => {
                            ui.label(RichText::new(installation.describe()).color(palette.warn));
                            ui.label("The release page has every download.");
                        }
                    }
                    ui.add_space(4.0);
                    ui.hyperlink_to("Release page", &release.page);
                }
                Stage::Downloading { done, total, .. } => {
                    let fraction = if *total > 0 { *done as f32 / *total as f32 } else { 0.0 };
                    ui.label("Downloading…");
                    ui.add(egui::ProgressBar::new(fraction).text(format!("{} of {}", megabytes(*done), megabytes(*total))));
                }
                Stage::Installing { release } => {
                    ui.horizontal(|ui| {
                        ui.spinner();
                        ui.label(format!("Installing SlopGene {}…", release.version));
                    });
                }
                Stage::Finished { release, done: Done::RunByHand { command } } => {
                    ui.label(format!("SlopGene {} is downloaded. To install it, run this in a terminal:", release.version));
                    ui.horizontal(|ui| {
                        ui.code(command);
                        if ui.small_button("Copy").clicked() {
                            ui.ctx().copy_text(command.clone());
                        }
                    });
                }
                Stage::Finished { release, .. } => {
                    ui.label(format!("SlopGene closes now and starts again as version {}.", release.version));
                }
                Stage::Failed { error } => {
                    ui.label(RichText::new(format!("The update did not go through: {error}")).color(palette.warn));
                    ui.hyperlink_to("Release page", slopgene_update::RELEASES_PAGE);
                }
            }
            ui.add_space(12.0);
            ui.checkbox(&mut self.at_start, "Look for a newer version when SlopGene starts");
            ui.add_space(8.0);
            ui.horizontal(|ui| match &stage {
                Stage::Available { asset: Some(_), .. } => {
                    let ready = unsaved.is_empty();
                    if ui.add_enabled(ready, egui::Button::new("Download and install")).clicked() {
                        install = true;
                    }
                    if ui.button("Later").clicked() {
                        close = true;
                    }
                }
                Stage::Downloading { .. } => {
                    if ui.button("Stop").clicked() {
                        stop = true;
                    }
                }
                Stage::Installing { .. } | Stage::Finished { .. } if !matches!(&stage, Stage::Finished { done: Done::RunByHand { .. }, .. }) => {}
                Stage::Failed { .. } => {
                    if ui.button("Try again").clicked() {
                        again = true;
                    }
                    if ui.button("Close").clicked() {
                        close = true;
                    }
                }
                Stage::UpToDate => {
                    if ui.button("Check again").clicked() {
                        again = true;
                    }
                    if ui.button("Close").clicked() {
                        close = true;
                    }
                }
                _ => {
                    if ui.button("Close").clicked() {
                        close = true;
                    }
                }
            });
        });
        if install {
            self.install();
        }
        if again {
            self.check();
        }
        if stop {
            self.stop.store(true, Ordering::SeqCst);
        }
        if close || modal.should_close() {
            self.open = false;
        }
    }
}

fn megabytes(bytes: u64) -> String {
    format!("{:.1} MB", bytes as f64 / 1_000_000.0)
}
