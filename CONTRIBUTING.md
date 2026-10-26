# Contributing to SlopGene

SlopGene is a plasmid and sequence editor in Rust: a core library, a command line, a Rhai scripting engine and an egui desktop app.
[docs/architecture.md](docs/architecture.md) explains how the crates fit together.
This guide covers building and testing, the rules the code follows, and how a release is made.

## Build and test

| Task | Command |
|---|---|
| Build everything | `cargo build --release` |
| All tests | `cargo test --workspace` |
| Format | `cargo fmt --all` (settings in `rustfmt.toml`) |
| Lints, as CI runs them | `cargo clippy --workspace --all-targets -- -D warnings` |
| Windows installer, locally | `cargo build --release -p slopgene-gui -p slopgene-cli`, then `cargo packager --release -p slopgene-gui --formats nsis` (cargo-packager 0.11.8) |
| App icon files, after editing `assets/icon/slopgene.svg` | `cargo run -p slopgene-gui --example make_icons` |

The minimum Rust version is 1.95, which egui 0.37 sets.
The README's [Building from source](README.md#building-from-source) section lists what each system needs besides Rust.
CI (`.github/workflows/ci.yml`) builds, lints and tests on Linux, macOS and Windows.
A pull request should pass the format check, the lints and the tests.

### Tests that use SnapGene

Three tests compare SlopGene with SnapGene.
Each of them passes without checking anything when SnapGene's files are absent.

- `crates/slopgene-core/tests/corpus.rs` reads and rewrites every SnapGene and GenBank file it finds and checks that the bytes come back unchanged.
  It reads the folders that `SLOPGENE_CORPUS` names (separated by `;` on Windows and `:` elsewhere), or SnapGene's installation folder when that variable is unset.
  Run it after any change to `io/snapgene.rs`, for example with `SLOPGENE_CORPUS='/path/to/maps' cargo test -p slopgene-core --test corpus -- --nocapture`.
  Every SnapGene file should come back identical byte for byte.
- `crates/slopgene-core/tests/snapgene_parity.rs` repeats the PCR, restriction cloning and Golden Gate workflows of SnapGene's sample project.
  It finds the sample project in SnapGene's installation folder, or where `SLOPGENE_SNAPGENE_SAMPLES` points.
- `crates/slopgene-core/tests/snapgene_opens.rs` has SnapGene's own command-line converter read `.dna` files that SlopGene edited, created and converted.
  It launches SnapGene, so it runs only when `SLOPGENE_SNAPGENE_EXE` names the program, as in `SLOPGENE_SNAPGENE_EXE='C:\Program Files\SnapGene\SnapGene.exe'`.
  Run it after any change to `io/snapgene.rs`.

Sequencing verification has synthetic structural cases (deletions, insertions, duplications, inversions, a dimer) in the tests of `crates/slopgene-core/src/align.rs`.

### Checking the app without a person at the screen

To check the app without a person at the screen, render a screenshot and look at it:
`slopgene-gui file.dna --view split --size 1360x860 --screenshot shot.png`.
`--palette` opens the command palette for the picture.
`--align read.gbk` opens the alignment view with that read, and `--select A..B` then scrolls its aligned bases to base A.
`--view alignment` opens the alignment view with nothing aligned.
A protocol file (`slopgene.protocol/1` JSON) opens in a protocol tab, and `--step N` shows its step N.
`--zoom` zooms into that step's card on the third frame.
`--after 45` then shows the detailed view, and `--after 6` catches the zoom on its way.
`--new-protocol` starts a cloning protocol from the files open, with the region marked by `--feature` or `--select` as the part to cut out.
`--enzymes NcoI,EagI` picks the insert's enzymes for it.
`--drag X1,Y1,X2,Y2` presses at the first point, moves to the second and holds, so the picture shows a drag on its way.
The drag ends just before the screenshot, so a zoom has finished first.
Its points are window coordinates, read off an earlier screenshot of the same window.
`--play SCRIPT` plays pointer and key actions into the window first, with a cursor drawn where the pointer is, and the screenshot is taken once the script has played.
A script is a list of steps separated by semicolons, such as `--play "wait 1; click 300,200; drag 100,100,400,100; key Escape"`, and SCRIPT may name a file that holds it.
The steps are `wait` (seconds), `move`, `click`, `drag`, `scroll` (`scroll X,Y,DY` turns the wheel at a point) and `key` (`Escape`, `Enter`, `Down`, `ctrl+shift+n`), and the pointer glides to each point.
`--record out.gif` saves the frames into an animated GIF while the script plays, in place of a screenshot, and quits when the script ends.
A frame that looks the same as the one before costs nothing, and the recording plays back at the speed the window ran.
`scripts/readme-animations.sh` records the README's animations this way, from the maps in `examples/maps`.

## Rules

- **Clean room.** SlopGene learns SnapGene's file format only from SnapGene documents and public parsers. Never disassemble or decompile SnapGene's programs, and never scan them for strings. Never commit SnapGene's data (feature library, vector maps, sample files) or anyone's private sequence files. Read them locally as test data only.
- **Core first.** Every calculation lives in `slopgene-core`. The command line, scripts and app only call it, so all three give the same answers.
- **Coordinates.** The core uses 0-based half-open `Span`s, which may run through the origin. Every surface people or scripts see uses 1-based inclusive positions. Convert only with `Span::from_1based`, `Span::to_1based` and `span::parse_range`.
- **Commands.** A command returns an `Output` that holds both its JSON and its text. Errors become `{"ok": false, "error": ...}` in JSON modes.
- **Writes.** Files are written through `io::write_atomic`.
- **Script functions.** A new function needs three things: its registration in `slopgene-script`, a row in the `API` table of `crates/slopgene-script/src/api.rs`, and the same signature in `docs/scripting.md`. Tests fail when the three disagree.
- **Agent guide.** `docs/agent-guide.md` is compiled into `slopgene guide`. Update it when a command changes.

## Style

- Loop indices are double letters that grow with nesting: `ii`, then `iii`, then `iiii`. Never use a bare `i`.
- Comments explain what each block does and why, so a reader can follow the code from the comments alone.
- Prose in the docs, commit messages and pull requests follows a few rules. Each sentence makes one claim. A paragraph opens with the sentence that sums it up. Commas only separate the items of a list, and semicolons are rare. There are no em dashes. Filler such as "it is worth noting" goes.
- Map colours come from the palette in `crates/slopgene-core/src/render/theme.rs`. The app's interface colours are in `crates/slopgene-gui/src/style.rs`.

## Making a release

The version lives in `[workspace.package]` of `Cargo.toml`.
To make a release, set the new version there, commit it, and push a tag that names it:

```
git tag v0.1.0
git push origin v0.1.0
```

`.github/workflows/release.yml` then builds the Windows, macOS and Linux installers on GitHub's machines.
It stops when the tag and the version in `Cargo.toml` differ.
The installers, the plain archives and their checksums go into a GitHub release, published at once as the latest release.
The README's download buttons and the app's update check read that release, so pushing the tag is the whole release.
Starting the workflow by hand from the Actions tab builds the same files as downloadable artifacts and makes no release.
The comments in `crates/slopgene-gui/Cargo.toml` show how to build each installer on your own computer.

## Licence

SlopGene is licensed under either the MIT licence ([LICENSE-MIT](LICENSE-MIT)) or the Apache License 2.0 ([LICENSE-APACHE](LICENSE-APACHE)), at your option.
Unless you state otherwise, a contribution you submit for inclusion in SlopGene is licensed in the same two ways, without any other terms or conditions.
