//! Driving the window without a person at it, for pictures and animations of the app.
//!
//! `--play SCRIPT` moves the pointer, clicks, drags and presses keys as the script says, and
//! draws a cursor on the window so a viewer can follow. `--record OUT.gif` saves the frames into
//! an animated GIF while the script plays, which is how the README's animations are made. The
//! script and the recording run on egui's clock, so a recording plays back at the speed the
//! window ran, however long each frame took to draw and save.
//!
//! A script is a list of steps separated by semicolons or line breaks, with points in window
//! coordinates from the top left corner: `wait 1.5`, `move X,Y`, `click X,Y`, `drag X1,Y1,X2,Y2`,
//! `scroll X,Y,DY`, `key Escape` and `key ctrl+shift+n`. The pointer glides to each point rather
//! than jumping there, so the eye can follow it.

use eframe::egui::{self, Color32, Modifiers, Pos2, Stroke, pos2, vec2};
use slopgene_core::io;
use std::borrow::Cow;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::{SyncSender, sync_channel};
use std::thread::JoinHandle;

/// How long the button stays down in a click.
const HOLD: f64 = 0.1;
/// How long a step rests after its last event, so the window has drawn the result before the
/// next step begins.
const SETTLE: f64 = 0.2;
/// How long the pointer takes from one point of a drag to the other.
const DRAG_SECONDS: f64 = 0.7;
/// How long the wheel turns in a scroll.
const SCROLL_SECONDS: f64 = 0.5;
/// How long the ripple around a click lasts.
const RIPPLE: f64 = 0.45;
/// The shortest time between two frames of a recording, a little under a thirtieth of a second so
/// that every second frame of a 60 Hz window is kept.
const CAPTURE_INTERVAL: f64 = 1.0 / 30.0 - 0.002;
/// The first frame a recording keeps. The window takes its size over the first frames.
const FIRST_FRAME: u64 = 6;
/// The palette index that keeps the pixel of the frame before, so a frame carries only what
/// changed.
const TRANSPARENT: u8 = 255;

/// One step of a play script.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Step {
    /// Do nothing for this many seconds.
    Wait(f64),
    /// Glide the pointer to this point.
    Move(Pos2),
    /// Glide the pointer to this point and click there.
    Click(Pos2),
    /// Glide the pointer to the first point, press, move to the second point and let go.
    Drag(Pos2, Pos2),
    /// Glide the pointer to this point and turn the wheel there, so that what is under it
    /// scrolls down by this many points (up when negative).
    Scroll(Pos2, f32),
    /// Press a key and let it go.
    Key(egui::Key, Modifiers),
}

/// A play script and how far it has got.
#[derive(Clone, Debug, Default)]
pub struct Play {
    steps: Vec<Step>,
    /// The step being played.
    index: usize,
    /// When that step began, and where the pointer was then.
    began: Option<f64>,
    from: Pos2,
    /// How far the step has got (0 on its way, 1 pressed or scrolling, 2 let go) and when it
    /// got there.
    phase: u8,
    phase_at: f64,
    /// How far a scroll has turned the wheel so far.
    scrolled: f32,
    /// Where the pointer is, once the first frame has placed it.
    pointer: Option<Pos2>,
    /// Recent presses, for the ripple drawn around a click.
    presses: Vec<(Pos2, f64)>,
}

impl Play {
    /// Read a script. Steps are separated by semicolons or line breaks, and `#` starts a comment
    /// that runs to the end of the line.
    pub fn parse(text: &str) -> Result<Play, String> {
        let mut steps = Vec::new();
        for raw in text.split([';', '\n']) {
            let line = raw.split('#').next().unwrap_or_default().trim();
            if line.is_empty() {
                continue;
            }
            let (word, rest) = line.split_once(char::is_whitespace).unwrap_or((line, ""));
            let rest = rest.trim();
            // The numbers after the word, separated by commas.
            let numbers = || -> Result<Vec<f32>, String> {
                rest.split(',').map(|n| n.trim().parse::<f32>().map_err(|_| format!("cannot read the numbers in \"{line}\""))).collect()
            };
            let step = match word.to_ascii_lowercase().as_str() {
                "wait" => Step::Wait(rest.parse::<f64>().map_err(|_| format!("wait takes a number of seconds: \"{line}\""))?),
                "move" => match numbers()?[..] {
                    [x, y] => Step::Move(pos2(x, y)),
                    _ => return Err(format!("move takes X,Y: \"{line}\"")),
                },
                "click" => match numbers()?[..] {
                    [x, y] => Step::Click(pos2(x, y)),
                    _ => return Err(format!("click takes X,Y: \"{line}\"")),
                },
                "drag" => match numbers()?[..] {
                    [x1, y1, x2, y2] => Step::Drag(pos2(x1, y1), pos2(x2, y2)),
                    _ => return Err(format!("drag takes X1,Y1,X2,Y2: \"{line}\"")),
                },
                "scroll" => match numbers()?[..] {
                    [x, y, by] => Step::Scroll(pos2(x, y), by),
                    _ => return Err(format!("scroll takes X,Y,DY: \"{line}\"")),
                },
                "key" => parse_key(rest).ok_or_else(|| format!("cannot read the key in \"{line}\": write it as Escape, Enter, Down or ctrl+shift+n"))?,
                _ => return Err(format!("unknown step \"{word}\" in \"{line}\": the steps are wait, move, click, drag, scroll and key")),
            };
            steps.push(step);
        }
        if steps.is_empty() {
            return Err("the script has no steps".into());
        }
        Ok(Play { steps, ..Play::default() })
    }

    /// Whether every step has been played.
    pub fn finished(&self) -> bool {
        self.index >= self.steps.len()
    }

    /// Feed the events of this frame into the window's input. The first frame puts the pointer
    /// near the bottom of the window, and each later frame moves it along the step being played.
    pub fn tick(&mut self, raw: &mut egui::RawInput) {
        let Some(time) = raw.time else { return };
        if self.pointer.is_none() {
            let screen = raw.screen_rect.unwrap_or(egui::Rect::from_min_size(Pos2::ZERO, vec2(800.0, 600.0)));
            let start = pos2(screen.center().x, screen.max.y - 40.0);
            self.pointer = Some(start);
            raw.events.push(egui::Event::PointerMoved(start));
        }
        self.presses.retain(|(_, at)| time - at < RIPPLE);
        let Some(step) = self.steps.get(self.index).copied() else { return };
        if self.began.is_none() {
            self.began = Some(time);
            self.from = self.pointer.unwrap_or(Pos2::ZERO);
            self.phase = 0;
            self.phase_at = time;
        }
        let elapsed = time - self.began.unwrap_or(time);
        if self.run(step, time, elapsed, &mut raw.events) {
            self.index += 1;
            self.began = None;
        }
    }

    /// Play one frame of a step. Returns true when the step is over.
    fn run(&mut self, step: Step, time: f64, elapsed: f64, events: &mut Vec<egui::Event>) -> bool {
        match step {
            Step::Wait(seconds) => elapsed >= seconds,
            Step::Move(to) => {
                self.glide(to, elapsed, events);
                elapsed >= travel_time(self.from, to) + SETTLE
            }
            Step::Click(at) => match self.phase {
                0 => {
                    if self.glide(at, elapsed, events) {
                        self.press(at, time, events);
                    }
                    false
                }
                1 => {
                    if time - self.phase_at >= HOLD {
                        self.release(at, time, events);
                    }
                    false
                }
                _ => time - self.phase_at >= SETTLE,
            },
            Step::Drag(from, to) => match self.phase {
                0 => {
                    if self.glide(from, elapsed, events) {
                        self.press(from, time, events);
                    }
                    false
                }
                1 => {
                    // The button went down at `from`, and the pointer moves on to `to`.
                    let t = ((time - self.phase_at) / DRAG_SECONDS).min(1.0);
                    self.move_to(from.lerp(to, ease(t)), events);
                    if t >= 1.0 {
                        self.release(to, time, events);
                    }
                    false
                }
                _ => time - self.phase_at >= SETTLE,
            },
            Step::Scroll(at, by) => match self.phase {
                0 => {
                    if self.glide(at, elapsed, events) {
                        self.phase = 1;
                        self.phase_at = time;
                        self.scrolled = 0.0;
                    }
                    false
                }
                1 => {
                    // The wheel turns smoothly, and what is under the pointer has moved by `by`
                    // at the end. egui's wheel runs the other way: a negative delta scrolls down.
                    let t = ((time - self.phase_at) / SCROLL_SECONDS).min(1.0);
                    let turned = ease(t) * by;
                    let step = turned - self.scrolled;
                    if step != 0.0 {
                        events.push(egui::Event::MouseWheel {
                            unit: egui::MouseWheelUnit::Point,
                            delta: vec2(0.0, -step),
                            phase: egui::TouchPhase::Move,
                            source: egui::MouseWheelSource::Wheel,
                            modifiers: Modifiers::default(),
                        });
                        self.scrolled = turned;
                    }
                    if t >= 1.0 {
                        self.phase = 2;
                        self.phase_at = time;
                    }
                    false
                }
                _ => time - self.phase_at >= SETTLE,
            },
            Step::Key(key, modifiers) => match self.phase {
                0 => {
                    events.push(key_event(key, modifiers, true));
                    self.phase = 1;
                    self.phase_at = time;
                    false
                }
                1 => {
                    events.push(key_event(key, modifiers, false));
                    self.phase = 2;
                    self.phase_at = time;
                    false
                }
                _ => time - self.phase_at >= SETTLE,
            },
        }
    }

    /// Move the pointer along its way from where the step began to `to`. Returns true once it
    /// has arrived.
    fn glide(&mut self, to: Pos2, elapsed: f64, events: &mut Vec<egui::Event>) -> bool {
        let t = (elapsed / travel_time(self.from, to)).min(1.0);
        self.move_to(self.from.lerp(to, ease(t)), events);
        t >= 1.0
    }

    fn move_to(&mut self, at: Pos2, events: &mut Vec<egui::Event>) {
        if self.pointer != Some(at) {
            self.pointer = Some(at);
            events.push(egui::Event::PointerMoved(at));
        }
    }

    fn press(&mut self, at: Pos2, time: f64, events: &mut Vec<egui::Event>) {
        events.push(button_event(at, true));
        self.presses.push((at, time));
        self.phase = 1;
        self.phase_at = time;
    }

    fn release(&mut self, at: Pos2, time: f64, events: &mut Vec<egui::Event>) {
        events.push(button_event(at, false));
        self.phase = 2;
        self.phase_at = time;
    }

    /// Draw the cursor where the pointer is, with a ripple around each recent click.
    pub fn draw(&self, ctx: &egui::Context, accent: Color32) {
        if let Some(pointer) = self.pointer {
            let time = ctx.input(|i| i.time);
            draw_cursor(ctx, pointer, &self.presses, time, accent);
        }
    }
}

/// How long the pointer takes to a point: a quarter of a second, and longer for a long way.
fn travel_time(from: Pos2, to: Pos2) -> f64 {
    (0.25 + f64::from(from.distance(to)) / 1800.0).min(0.9)
}

/// Slow at both ends, so a movement starts and stops gently.
fn ease(t: f64) -> f32 {
    egui::emath::easing::cubic_in_out(t as f32)
}

/// A key name with modifiers, such as `Escape` or `ctrl+shift+n`.
fn parse_key(text: &str) -> Option<Step> {
    let mut modifiers = Modifiers::default();
    let mut key = None;
    let parts: Vec<&str> = text.split('+').map(str::trim).filter(|p| !p.is_empty()).collect();
    for (ii, part) in parts.iter().enumerate() {
        let last = ii + 1 == parts.len();
        match part.to_ascii_lowercase().as_str() {
            "ctrl" | "control" if !last => {
                modifiers.ctrl = true;
                modifiers.command = true;
            }
            "shift" if !last => modifiers.shift = true,
            "alt" if !last => modifiers.alt = true,
            "cmd" | "command" if !last => modifiers.command = true,
            _ => key = egui::Key::from_name(part),
        }
    }
    Some(Step::Key(key?, modifiers))
}

fn key_event(key: egui::Key, modifiers: Modifiers, pressed: bool) -> egui::Event {
    egui::Event::Key { key, physical_key: None, pressed, repeat: false, modifiers }
}

fn button_event(pos: Pos2, pressed: bool) -> egui::Event {
    egui::Event::PointerButton { pos, button: egui::PointerButton::Primary, pressed, modifiers: Modifiers::default() }
}

/// An arrow cursor at `at`, over a soft halo so it is easy to find, with a ring spreading from
/// each recent press.
pub fn draw_cursor(ctx: &egui::Context, at: Pos2, presses: &[(Pos2, f64)], time: f64, accent: Color32) {
    let painter = ctx.layer_painter(egui::LayerId::new(egui::Order::Tooltip, egui::Id::unique("play-cursor")));
    for (center, pressed_at) in presses {
        let age = ((time - pressed_at) / RIPPLE) as f32;
        if (0.0..1.0).contains(&age) {
            let color = Color32::from_rgba_unmultiplied(accent.r(), accent.g(), accent.b(), (220.0 * (1.0 - age)) as u8);
            painter.circle_stroke(*center, 8.0 + 26.0 * age, Stroke::new(2.5, color));
        }
    }
    painter.circle_filled(at, 13.0, Color32::from_rgba_unmultiplied(255, 204, 0, 70));
    // The outline of a standard arrow, scaled up a little, filled as five triangles (the arrow is
    // not convex) and then drawn round in white.
    let scale = 1.25;
    let outline: Vec<Pos2> = [(0.0, 0.0), (0.0, 16.0), (4.0, 12.0), (7.0, 18.0), (9.0, 17.0), (6.0, 11.0), (11.0, 11.0)]
        .iter()
        .map(|(x, y)| at + vec2(x * scale, y * scale))
        .collect();
    let mut mesh = egui::Mesh::default();
    for point in &outline {
        mesh.colored_vertex(*point, Color32::BLACK);
    }
    for [a, b, c] in [[0, 1, 2], [0, 2, 5], [0, 5, 6], [2, 3, 4], [2, 4, 5]] {
        mesh.add_triangle(a, b, c);
    }
    painter.add(egui::Shape::mesh(mesh));
    painter.add(egui::Shape::closed_line(outline, Stroke::new(1.2, Color32::WHITE)));
}

/// A frame of the window, with the time it was drawn.
pub struct Captured {
    pub time: f64,
    pub image: Arc<egui::ColorImage>,
}

/// What a recording came to.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Report {
    /// Frames written to the file. Frames that look the same as the one before cost nothing.
    pub frames: usize,
    /// How long the recording plays for.
    pub seconds: f64,
    /// Frames lost because the encoder could not keep up.
    pub lost: usize,
}

/// Frames of the window on their way into an animated GIF. The frames are asked for as the window
/// draws, and a thread encodes them as they arrive, so the window keeps its pace.
pub struct Recorder {
    path: PathBuf,
    /// Hands frames to the encoding thread. Dropping it ends the recording once the frames
    /// already asked for have arrived, since each request holds a copy until it is answered.
    sender: Option<SyncSender<Captured>>,
    worker: Option<JoinHandle<Result<Report, String>>>,
    /// When the last frame was asked for.
    last: Option<f64>,
    lost: Arc<AtomicUsize>,
    outcome: Option<Result<Report, String>>,
}

impl Recorder {
    pub fn new(path: PathBuf) -> Recorder {
        let (sender, receiver) = sync_channel::<Captured>(90);
        let lost = Arc::new(AtomicUsize::new(0));
        let file = path.clone();
        let lost_count = lost.clone();
        let worker = std::thread::spawn(move || -> Result<Report, String> {
            let (bytes, mut report) = encode(receiver.into_iter())?;
            io::write_atomic(&file, &bytes).map_err(|e| format!("{}: {e}", file.display()))?;
            report.lost = lost_count.load(Ordering::SeqCst);
            Ok(report)
        });
        Recorder { path, sender: Some(sender), worker: Some(worker), last: None, lost, outcome: None }
    }

    /// Ask for this frame when it is time for one, or with `stop` end the recording. The file is
    /// written once the frames already asked for have arrived.
    pub fn tick(&mut self, ctx: &egui::Context, time: f64, frame: u64, stop: bool) {
        if stop {
            self.sender = None;
        } else if frame >= FIRST_FRAME
            && self.last.is_none_or(|last| time - last >= CAPTURE_INTERVAL)
            && let Some(sender) = &self.sender
        {
            self.last = Some(time);
            let sender = sender.clone();
            let lost = self.lost.clone();
            ctx.request_screenshot(move |image| {
                if sender.try_send(Captured { time, image }).is_err() {
                    lost.fetch_add(1, Ordering::SeqCst);
                }
            });
        }
        if self.outcome.is_none()
            && self.worker.as_ref().is_some_and(JoinHandle::is_finished)
            && let Some(worker) = self.worker.take()
        {
            self.outcome = Some(worker.join().unwrap_or_else(|_| Err("the encoder stopped unexpectedly".into())));
        }
    }

    /// Whether the file has been written, or the recording failed.
    pub fn finished(&self) -> bool {
        self.outcome.is_some()
    }

    /// What became of the recording, in a line for the terminal.
    pub fn summary(&self) -> String {
        match &self.outcome {
            Some(Ok(report)) => {
                let lost = if report.lost > 0 { format!(", {} frames lost to a slow encoder", report.lost) } else { String::new() };
                format!("recorded {} frames over {:.1} s into {}{lost}", report.frames, report.seconds, self.path.display())
            }
            Some(Err(e)) => format!("recording failed: {e}"),
            None => "recording".into(),
        }
    }
}

/// A frame ready for the file: the part of the window that changed, as palette indices, with
/// the pixels that did not change marked transparent.
struct Prepared {
    left: u16,
    top: u16,
    width: u16,
    height: u16,
    /// 256 colours, index 255 being the transparent one.
    palette: Vec<u8>,
    indices: Vec<u8>,
    first: bool,
}

/// A frame waiting for the next one, which decides how long it is shown.
struct Pending {
    time: f64,
    pixels: Vec<u32>,
    prepared: Prepared,
}

/// Encode frames into a GIF. The first frame sets the size, and frames of another size are left
/// out. A frame that looks the same as the one before only lengthens that frame's time on
/// screen, and any other frame carries just the rectangle that changed.
pub fn encode(frames: impl Iterator<Item = Captured>) -> Result<(Vec<u8>, Report), String> {
    let failed = |e: gif::EncodingError| format!("GIF encoding failed: {e}");
    let mut encoder: Option<gif::Encoder<Vec<u8>>> = None;
    let mut size = [0usize; 2];
    // The last frame written, the frame waiting to be written, and the clock.
    let mut written: Vec<u32> = Vec::new();
    let mut pending: Option<Pending> = None;
    let (mut first_time, mut last_time) = (0.0, 0.0);
    let mut carry = 0.0;
    let mut count = 0;
    for captured in frames {
        let pixels: Vec<u32> = captured.image.pixels.iter().map(|p| (u32::from(p.r()) << 16) | (u32::from(p.g()) << 8) | u32::from(p.b())).collect();
        if encoder.is_none() {
            size = captured.image.size;
            if size[0] == 0 || size[1] == 0 || size[0] > usize::from(u16::MAX) || size[1] > usize::from(u16::MAX) {
                return Err(format!("a frame of {} by {} pixels cannot go into a GIF", size[0], size[1]));
            }
            let mut fresh = gif::Encoder::new(Vec::new(), size[0] as u16, size[1] as u16, &[]).map_err(failed)?;
            fresh.set_repeat(gif::Repeat::Infinite).map_err(failed)?;
            encoder = Some(fresh);
            first_time = captured.time;
        }
        if captured.image.size != size || pixels.len() != size[0] * size[1] {
            continue;
        }
        last_time = captured.time;
        if let Some(waiting) = &pending
            && waiting.pixels == pixels
        {
            continue;
        }
        if let Some(waiting) = pending.take()
            && let Some(encoder) = encoder.as_mut()
        {
            let delay = centiseconds(captured.time - waiting.time, &mut carry);
            write_frame(encoder, &waiting.prepared, delay).map_err(failed)?;
            written = waiting.pixels;
            count += 1;
        }
        let before = if written.is_empty() { None } else { Some(written.as_slice()) };
        if let Some(prepared) = prepare(before, &pixels, size[0], size[1]) {
            pending = Some(Pending { time: captured.time, pixels, prepared });
        }
    }
    let Some(mut encoder) = encoder else { return Err("no frames were recorded".into()) };
    if let Some(waiting) = pending.take() {
        // The last frame stays for the time the recording left it, and at least a second.
        let delay = centiseconds((last_time - waiting.time).max(1.0), &mut carry);
        write_frame(&mut encoder, &waiting.prepared, delay).map_err(failed)?;
        count += 1;
    }
    let bytes = encoder.into_inner().map_err(failed)?;
    Ok((bytes, Report { frames: count, seconds: last_time - first_time, lost: 0 }))
}

/// A time in hundredths of a second, carrying the rounding error into the next one so the
/// recording keeps time over many frames. Browsers show a frame for at least two hundredths.
fn centiseconds(seconds: f64, carry: &mut f64) -> u16 {
    let exact = seconds * 100.0 + *carry;
    let rounded = exact.round().max(2.0);
    *carry = exact - rounded;
    rounded.min(f64::from(u16::MAX)) as u16
}

fn write_frame(encoder: &mut gif::Encoder<Vec<u8>>, prepared: &Prepared, delay: u16) -> Result<(), gif::EncodingError> {
    let frame = gif::Frame {
        delay,
        dispose: gif::DisposalMethod::Keep,
        transparent: if prepared.first { None } else { Some(TRANSPARENT) },
        top: prepared.top,
        left: prepared.left,
        width: prepared.width,
        height: prepared.height,
        palette: Some(prepared.palette.clone()),
        buffer: Cow::Borrowed(&prepared.indices),
        ..gif::Frame::default()
    };
    encoder.write_frame(&frame)
}

/// The rectangle of `after` that differs from `before`, as palette indices, with the pixels
/// inside it that did not change marked transparent. None when nothing changed. Without a
/// `before`, the whole frame.
fn prepare(before: Option<&[u32]>, after: &[u32], width: usize, height: usize) -> Option<Prepared> {
    let changed = |ii: usize| before.is_none_or(|b| b[ii] != after[ii]);
    // The rows and columns that changed.
    let (mut left, mut top, mut right, mut bottom) = (usize::MAX, usize::MAX, 0, 0);
    for ii in 0..height {
        for iii in 0..width {
            if changed(ii * width + iii) {
                left = left.min(iii);
                right = right.max(iii);
                top = top.min(ii);
                bottom = bottom.max(ii);
            }
        }
    }
    if left == usize::MAX {
        return None;
    }
    // The colours of the changed pixels, by how often each occurs.
    let mut counts: HashMap<u32, u32> = HashMap::new();
    for ii in top..=bottom {
        for iii in left..=right {
            let at = ii * width + iii;
            if changed(at) {
                *counts.entry(after[at]).or_insert(0) += 1;
            }
        }
    }
    let (palette, index) = quantize(&counts);
    let mut indices = Vec::with_capacity((right - left + 1) * (bottom - top + 1));
    for ii in top..=bottom {
        for iii in left..=right {
            let at = ii * width + iii;
            indices.push(if changed(at) { index.get(&after[at]).copied().unwrap_or(0) } else { TRANSPARENT });
        }
    }
    Some(Prepared {
        left: left as u16,
        top: top as u16,
        width: (right - left + 1) as u16,
        height: (bottom - top + 1) as u16,
        palette,
        indices,
        first: before.is_none(),
    })
}

/// Fold the colours into at most 255, by median cut: the box of colours with the widest spread
/// along a channel is split at its median, by weight, until there are 255 boxes or every box
/// holds one colour. Each colour then maps to the weighted mean of its box. Returns the palette
/// as 256 RGB triples (the last one spare, for transparency) and the index of every colour.
fn quantize(counts: &HashMap<u32, u32>) -> (Vec<u8>, HashMap<u32, u8>) {
    let channel = |rgb: u32, ii: usize| ((rgb >> (16 - 8 * ii)) & 0xFF) as i32;
    let mut boxes: Vec<Vec<(u32, u32)>> = vec![counts.iter().map(|(&rgb, &n)| (rgb, n)).collect()];
    boxes[0].sort_unstable();
    while boxes.len() < 255 {
        // The box with the widest channel, and that channel.
        let mut widest: Option<(usize, usize, i32)> = None;
        for (ii, colors) in boxes.iter().enumerate() {
            if colors.len() < 2 {
                continue;
            }
            for iii in 0..3 {
                let (low, high) = colors.iter().fold((255, 0), |(low, high), &(rgb, _)| (low.min(channel(rgb, iii)), high.max(channel(rgb, iii))));
                if widest.is_none_or(|(_, _, spread)| high - low > spread) {
                    widest = Some((ii, iii, high - low));
                }
            }
        }
        let Some((at, by, _)) = widest else { break };
        let mut colors = boxes.swap_remove(at);
        colors.sort_by_key(|&(rgb, _)| channel(rgb, by));
        // Split where half the weight lies, keeping at least one colour on each side.
        let total: u64 = colors.iter().map(|&(_, n)| u64::from(n)).sum();
        let mut seen = 0;
        let mut split = colors.len() - 1;
        for (ii, &(_, n)) in colors.iter().enumerate() {
            seen += u64::from(n);
            if seen * 2 >= total {
                split = ii + 1;
                break;
            }
        }
        let split = split.clamp(1, colors.len() - 1);
        let rest = colors.split_off(split);
        boxes.push(colors);
        boxes.push(rest);
    }
    let mut palette = Vec::with_capacity(256 * 3);
    let mut index = HashMap::with_capacity(counts.len());
    for (ii, colors) in boxes.iter().enumerate() {
        let weight: u64 = colors.iter().map(|&(_, n)| u64::from(n)).sum::<u64>().max(1);
        for iii in 0..3 {
            let sum: u64 = colors.iter().map(|&(rgb, n)| u64::from(n) * channel(rgb, iii) as u64).sum();
            palette.push(sum.div_ceil(weight).min(255) as u8);
        }
        for &(rgb, _) in colors {
            index.insert(rgb, ii as u8);
        }
    }
    palette.resize(256 * 3, 0);
    (palette, index)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_script_reads_into_steps() {
        let play =
            Play::parse("wait 1.5; click 10,20 # the gene\nDRAG 1,2, 3,4; key Escape; key ctrl+shift+n; # a comment\n move 5,6; scroll 7,8,-90").unwrap();
        let ctrl_shift = Modifiers { ctrl: true, command: true, shift: true, ..Modifiers::default() };
        assert_eq!(
            play.steps,
            vec![
                Step::Wait(1.5),
                Step::Click(pos2(10.0, 20.0)),
                Step::Drag(pos2(1.0, 2.0), pos2(3.0, 4.0)),
                Step::Key(egui::Key::Escape, Modifiers::default()),
                Step::Key(egui::Key::N, ctrl_shift),
                Step::Move(pos2(5.0, 6.0)),
                Step::Scroll(pos2(7.0, 8.0), -90.0),
            ]
        );
        assert!(Play::parse("").is_err());
        assert!(Play::parse("click 1").unwrap_err().contains("X,Y"));
        assert!(Play::parse("jump 1,2").unwrap_err().contains("unknown step"));
        assert!(Play::parse("key Hyperspace").unwrap_err().contains("key"));
    }

    /// Run a script frame by frame at 60 Hz and gather the events it feeds the window.
    fn events_of(script: &str, seconds: f64) -> (Vec<egui::Event>, bool) {
        let mut play = Play::parse(script).unwrap();
        let mut events = Vec::new();
        let frames = (seconds * 60.0) as usize;
        for ii in 0..frames {
            let mut raw = egui::RawInput {
                time: Some(ii as f64 / 60.0),
                screen_rect: Some(egui::Rect::from_min_size(Pos2::ZERO, vec2(800.0, 600.0))),
                ..Default::default()
            };
            play.tick(&mut raw);
            events.extend(raw.events);
        }
        (events, play.finished())
    }

    #[test]
    fn a_click_glides_there_presses_and_lets_go() {
        let (events, finished) = events_of("click 100,50", 2.0);
        assert!(finished);
        let presses: Vec<(Pos2, bool)> =
            events.iter().filter_map(|e| if let egui::Event::PointerButton { pos, pressed, .. } = e { Some((*pos, *pressed)) } else { None }).collect();
        assert_eq!(presses, vec![(pos2(100.0, 50.0), true), (pos2(100.0, 50.0), false)]);
        // The pointer moved several times on its way, and arrived before the press.
        let moves: Vec<Pos2> = events.iter().filter_map(|e| if let egui::Event::PointerMoved(p) = e { Some(*p) } else { None }).collect();
        assert!(moves.len() > 5);
        assert_eq!(moves.last(), Some(&pos2(100.0, 50.0)));
        let press_at = events.iter().position(|e| matches!(e, egui::Event::PointerButton { pressed: true, .. })).unwrap();
        assert!(matches!(events[press_at - 1], egui::Event::PointerMoved(p) if p == pos2(100.0, 50.0)));
    }

    #[test]
    fn a_drag_holds_the_button_while_it_moves_and_a_key_goes_down_and_up() {
        let (events, finished) = events_of("drag 10,10,200,100; key Enter", 3.0);
        assert!(finished);
        let down = events.iter().position(|e| matches!(e, egui::Event::PointerButton { pressed: true, .. })).unwrap();
        let up = events.iter().position(|e| matches!(e, egui::Event::PointerButton { pressed: false, .. })).unwrap();
        let between = events[down..up].iter().filter(|e| matches!(e, egui::Event::PointerMoved(_))).count();
        assert!(between > 10, "{between} moves between press and release");
        assert!(matches!(events[up], egui::Event::PointerButton { pos, .. } if pos == pos2(200.0, 100.0)));
        let keys: Vec<bool> =
            events.iter().filter_map(|e| if let egui::Event::Key { key: egui::Key::Enter, pressed, .. } = e { Some(*pressed) } else { None }).collect();
        assert_eq!(keys, vec![true, false]);
    }

    #[test]
    fn a_scroll_turns_the_wheel_by_the_whole_distance() {
        let (events, finished) = events_of("scroll 50,50,300", 2.0);
        assert!(finished);
        let turned: f32 = events.iter().filter_map(|e| if let egui::Event::MouseWheel { delta, .. } = e { Some(delta.y) } else { None }).sum();
        assert!((turned + 300.0).abs() < 0.01, "the wheel turned {turned}");
        let first_turn = events.iter().position(|e| matches!(e, egui::Event::MouseWheel { .. })).unwrap();
        assert!(matches!(events[first_turn - 1], egui::Event::PointerMoved(p) if p == pos2(50.0, 50.0)));
    }

    #[test]
    fn a_wait_takes_its_time() {
        let (_, finished) = events_of("wait 1", 0.9);
        assert!(!finished);
        let (_, finished) = events_of("wait 1", 1.1);
        assert!(finished);
    }

    #[test]
    fn many_colours_fold_into_a_palette_and_few_stay_as_they_are() {
        let mut counts = HashMap::new();
        for ii in 0..2000u32 {
            counts.insert((ii * 7919) & 0xFF_FF_FF, 1 + ii % 5);
        }
        let (palette, index) = quantize(&counts);
        assert_eq!(palette.len(), 256 * 3);
        assert_eq!(index.len(), 2000);
        assert!(index.values().all(|&ii| ii < 255));
        let few: HashMap<u32, u32> = [(0xFF0000, 10), (0x00FF00, 20), (0x0000FF, 30)].into_iter().collect();
        let (palette, index) = quantize(&few);
        for (rgb, ii) in index {
            let at = 3 * usize::from(ii);
            assert_eq!(u32::from(palette[at]) << 16 | u32::from(palette[at + 1]) << 8 | u32::from(palette[at + 2]), rgb);
        }
    }

    fn frame(time: f64, size: [usize; 2], color: Color32, spot: Option<(usize, usize, Color32)>) -> Captured {
        let mut image = egui::ColorImage::filled(size, color);
        if let Some((x, y, c)) = spot {
            image.pixels[y * size[0] + x] = c;
        }
        Captured { time, image: Arc::new(image) }
    }

    #[test]
    fn a_recording_keeps_only_the_frames_that_change_and_their_times() {
        let size = [40, 30];
        let frames = vec![
            frame(0.0, size, Color32::WHITE, None),
            frame(0.5, size, Color32::WHITE, None),
            frame(1.0, size, Color32::WHITE, Some((10, 5, Color32::RED))),
            frame(1.2, size, Color32::WHITE, Some((10, 5, Color32::RED))),
            frame(1.5, [10, 10], Color32::BLACK, None),
            frame(2.0, size, Color32::WHITE, Some((30, 20, Color32::BLUE))),
            frame(2.5, size, Color32::WHITE, Some((30, 20, Color32::BLUE))),
        ];
        let (bytes, report) = encode(frames.into_iter()).unwrap();
        assert_eq!(report.frames, 3);
        assert!((report.seconds - 2.5).abs() < 1e-9);
        let mut options = gif::DecodeOptions::new();
        options.set_color_output(gif::ColorOutput::Indexed);
        let mut decoder = options.read_info(bytes.as_slice()).unwrap();
        assert_eq!((decoder.width(), decoder.height()), (40, 30));
        let mut seen = Vec::new();
        while let Some(frame) = decoder.read_next_frame().unwrap() {
            seen.push((frame.delay, frame.left, frame.top, frame.width, frame.height, frame.transparent));
        }
        // The first frame is whole and stays a second. The second carries one pixel and stays
        // another second. The last one stays at least a second too.
        assert_eq!(seen[0], (100, 0, 0, 40, 30, None));
        assert_eq!(seen[1], (100, 10, 5, 1, 1, Some(TRANSPARENT)));
        assert_eq!(seen[2].0, 100);
        assert_eq!((seen[2].1, seen[2].2, seen[2].3, seen[2].4), (10, 5, 21, 16));
    }

    #[test]
    fn a_changed_rectangle_keeps_the_pixels_around_the_change() {
        let before = vec![0u32; 16];
        let mut after = before.clone();
        after[5] = 0xFF0000;
        after[10] = 0x00FF00;
        let prepared = prepare(Some(&before), &after, 4, 4).unwrap();
        assert_eq!((prepared.left, prepared.top, prepared.width, prepared.height), (1, 1, 2, 2));
        // Only the two changed pixels have colours; the other two keep what was there.
        assert_eq!(prepared.indices.iter().filter(|&&ii| ii == TRANSPARENT).count(), 2);
        assert!(prepare(Some(&after), &after, 4, 4).is_none());
        assert!(prepare(None, &after, 4, 4).unwrap().first);
    }

    #[test]
    fn hundredths_carry_their_rounding() {
        let mut carry = 0.0;
        let total: u32 = (0..100).map(|_| u32::from(centiseconds(0.033, &mut carry))).sum();
        assert_eq!(total, 330);
        assert_eq!(centiseconds(0.001, &mut carry), 2);
    }
}
