# Architecture

SlopGene is one Rust library with three thin programs around it.
The library owns the document model, the file formats, every algorithm and the layout of maps.
The command line, the scripting engine and the desktop app only translate between that library and their users.
A digest, a PCR or a map therefore gives the same result in every program.

```
                 slopgene-core
   model · formats · algorithms · map layout (Scene)
        ▲                 ▲                  ▲
        │                 │                  │
 slopgene-script   slopgene-cli      slopgene-gui
   (Rhai engine)    (`slopgene`)     (`slopgene-gui`, egui)
        ▲                 │                  │
        └─────────────────┴──────────────────┘
          both programs run scripts through it
```

## Crates

| Crate | Depends on | Purpose |
|---|---|---|
| `slopgene-core` | serde, lzma-rs, thiserror | Document model, file formats, algorithms, map layout |
| `slopgene-script` | core, rhai, resvg | The Rhai engine with SlopGene's functions, and PNG export |
| `slopgene-cli` | core, script, clap, rayon | The `slopgene` command |
| `slopgene-gui` | core, script, eframe (egui with wgpu), rfd | The `slopgene-gui` desktop app |

The core has no graphics, terminal or network dependency.
It compiles anywhere Rust does, which keeps the algorithms testable without a display.

## The document model

`model.rs` defines a `Document`: one DNA, RNA or protein molecule with its annotations.
It holds the name, the molecule type, the topology, the bases, the features, the primers, the end structure of a linear fragment, methylation settings and file metadata.
A feature is a list of segments, so it can be split (an intron, or a feature that runs through the origin).
Each feature and segment keeps the SnapGene attributes SlopGene does not model, in their original order.

### Coordinates

Inside the library, positions are 0-based and half-open.
A `Span` (`span.rs`) runs from `start` to `end`.
When `end` exceeds the sequence length, the span continues through the origin of a circular molecule.
`Span::pieces` turns such a span into its one or two plain stretches.

Every surface a person or script sees uses 1-based inclusive positions, like SnapGene and GenBank.
The conversion happens in one place (`Span::from_1based`, `Span::to_1based`, `span::parse_range`).
Code outside the core never converts positions by hand.

## File formats

`io/mod.rs` detects a file's format from its first bytes and dispatches to a reader.
Writers go through `io::write_atomic`, which writes a temporary file beside the target and renames it.
A program that reads the file at the same moment sees either the old file or the new one.

| Module | Formats |
|---|---|
| `io/snapgene.rs` | SnapGene `.dna`, `.rna` and `.prot` |
| `io/genbank.rs` | GenBank flat files, single or multi-record |
| `io/fasta.rs` | FASTA, single or multi-record |
| `json.rs` | SlopGene document JSON (`slopgene.document/1`) |

### SnapGene fidelity

A SnapGene file is a sequence of typed packets.
Several packets hold XML, and one holds the edit history as xz-compressed XML.
`io/snapgene.rs` names each packet type and holds the rules for reading and writing it.

The reader keeps every packet's original bytes next to the parsed model.
The writer compares each packet's content with the model.
A packet whose content did not change is written from its original bytes.
A file that was only read is therefore written back byte for byte, and the corpus test checks this on every SnapGene file it can find.

After an edit to the bases, some packets describe a sequence that no longer exists.
These are the history, cached analyses and alignment data.
The writer drops them, and SnapGene rebuilds what it needs when it opens the file.
The list lives in `STALE_AFTER_SEQUENCE_EDIT` in `io/snapgene.rs`.

## Editing

`edit.rs` implements every change to the bases: insert, delete, replace, invert, reverse-complement, set the origin and change the topology.
Each edit moves, clips or splits the features and primers so they keep covering the same bases.
A replacement of the same length is a point mutation and moves nothing.
After an edit, the caller recomputes where primers bind (`primers::refresh_binding_sites`).

## Algorithms

| Module | What it computes |
|---|---|
| `alphabet.rs` | IUPAC codes, complements, cleaning of pasted sequence |
| `enzymes.rs` | The enzyme database, enzyme sets, site finding and cut positions |
| `digest.rs` | Fragments with exact sticky ends, ligation and circularisation |
| `primers.rs` | Nearest-neighbour Tm, binding sites with mismatches, dimer and hairpin checks, primer design |
| `pcr.rs` | PCR products, including 5' tails and the template's features |
| `cloning.rs` | Restriction cloning, Gibson / HiFi / In-Fusion assembly, Golden Gate |
| `excise.rs` | The enzymes that cut nearest a marked region on each side, and what a chosen pair cuts out, checked on the vector |
| `protocol.rs` | Cloning protocols: steps that open, digest, amplify by PCR, pick, dephosphorylate, flip and ligate, run in order with their choices and problems |
| `translate.rs` | Translation with the NCBI genetic codes, open reading frames |
| `search.rs` | DNA motifs with mismatches on both strands, protein motifs in six frames |
| `align.rs` | Read-to-reference alignment, the protein effect of each difference, and the aligned bases laid out in lines |
| `annotate.rs` | Feature detection from a library of trusted annotated maps |
| `gel.rs` | Band positions on a virtual agarose gel |
| `stats.rs` | Molecular weight, isoelectric point, extinction coefficient |

A few choices shape these modules.

- **Restriction enzymes** come from REBASE and are compiled into the program from `data/enzymes.tsv`. The environment variable `SLOPGENE_ENZYMES` points at a replacement file, which `slopgene enzymes import-rebase` writes.
- **Site finding** indexes the enzymes by the first four bases of their site in both orientations. Each position of the sequence is compared only with the enzymes whose site can start there. A test compares the result with an exhaustive scan on random sequences.
- **Sticky ends** are stored on the fragment documents themselves (`Document::ends`). Ligation checks that the overhangs pair and merges them, so cloning products contain each base once.
- **Alignment** follows the seed, chain and align method of SnapGene's "Align to Reference" and minimap2. Shared 15-mers are merged into runs of exact matches, the runs are chained with minimap2's gap cost, and the stretches between runs are aligned with affine gaps (Gotoh). A chain carries across an indel of any size, so changes of 50 bp or more become structural events (deletion, insertion, duplication, inversion, replacement) rather than a collapse in identity. Whole-plasmid reads are aligned twice, the second time opened inside their longest exact match, so a change across the read's arbitrary start is found too. The alignment is kept column by column (`Alignment::aligned`). `Alignment::blocks` lays the columns out in lines, with each structural change in a block of its own. The app's alignment view and `slopgene align --show` both draw those blocks. `Alignment::track` turns an alignment into an `AlignmentTrack` for `MapOptions::alignments`. A track holds the stretch the read covers, cut into segments (aligned, turned round, missing) by the structural changes, and a mark for each difference. `Alignment::tally` counts the reference bases the read pairs, lacks and turns round, and the bases only the read has. Circular maps draw tracks as rings outside the primers, and linear maps as lines under the features. An aligned segment is a solid band and a missing one a thin line, so a ring reads as the other sequence, broken where it lacks the map's bases. The comparison colours (green aligned, red missing, orange changed, violet inserted) sit in `render/theme.rs` beside the palette and stay out of it on purpose: every palette colour already marks a feature type.
- **Excisions** start from a region a person marks, not from enzyme names. `excise::nearest_cuts` lists the enzymes of a set that cut nearest the region on each side, each at its cut nearest the region, with the end that cut leaves and how often the enzyme cuts. Enzymes that cut inside the region are set apart. Nothing is paired or ranked beyond distance, because the person picks one enzyme on each side. `excise::excision` then digests with the chosen pair, and the fragment that holds the region gives the exact ends, length and bases carried along. With a vector, it also checks that the pair cuts the vector once each (one enzyme once or twice) and that the backbone's ends pair with the insert's. `excise::can_join` says whether an end can join a given piece, which the app shows as a tick beside an enzyme. The app's bench, `slopgene excise` and the script functions `nearest_cuts` and `excision` all call these, so they agree.
- **Protocols** are data (`Protocol`, saved as `slopgene.protocol/1` JSON), and `protocol::run` runs them. Each step names its inputs by the ids of earlier steps, so a protocol is a chain of steps that can branch and join. A run gives every step its molecules, findings (cut sites, the fragment kept, junctions), notes and state. A digest may name its part (`around`). The run then finds that part on the molecule, says which fragment holds it, and warns when an enzyme cuts inside it. Every fragment keeps its `Origin`, the stretch of the molecule it came from. A ligation records the stretch of the product each input became (`Findings::parts`), so a picture can draw where each piece went. A step that several options fit stops with `NeedsChoice` and lists the options with their products. A digest without enzymes stops there too, with no options, until it gets them. An `amplify` step copies a part by PCR: it designs primers for the part with `primer_design::design_pairs`, with the sites of its tail enzymes on their 5' tails (`site_tail`), makes the product with `pcr::pcr`, and records the primers (`Findings::primers`), the pairs designed (`Findings::primer_pairs`) and the tails' lengths. The primers are designed anew on each run until a pair is kept in the step (`set_primers`), after which every run checks them against the part. A PCR product carries no 5' phosphates, and a digest keeps an uncut end's phosphate state, so only the cut ends of the piece between the tails carry them. A PCR and the digest that trims its tails share their enzymes: `set_enzymes` on either changes both, a PCR's new tails lead the vector as an insert's cut does, and a vector's new enzymes lead the PCRs whose tails had the same enzymes, since a tail can carry any site. `copy_out` and `cut_out` switch a lane between a cut and a PCR. Choosing writes the choice into the step's parameters, so the file records it. A failure marks the steps that depend on it as skipped rather than stopping the run. Ligation checks each junction for pairing overhangs and a 5' phosphate on at least one end. It also lists what else its tube makes (`Findings::byproducts`): every other way round of the pieces, and for a circle the backbone (the first piece) closed with fewer of the other pieces or with them in another order, up to three inserts. Each `Byproduct` holds the backbone, since only such a circle gives colonies, and one that holds every piece in order carries the `Pin` that makes it the product. `protocol::story` tells each step in plain words, and `protocol::summary` tells the whole protocol, so the app, the command line and JSON describe it alike. The `Loader` trait finds the maps: `Folder` reads files beside the protocol file, and the app prefers an open tab of the same file. The edits a person makes on the app's bench are methods of `Protocol` (`add_source`, `keep_params`, `put_piece`, `join_into`, `add_to_join`, `take_out_of_join`, `remove_source`, `dephosphorylate_input`, `set_orientation`, `set_part`, `set_enzymes`, `set_tails`, `set_primers`, `copy_out`, `cut_out` and `leave_out`). `Change` names each of them as data, and `Protocol::make` makes one and says what it did. `change_file` reads a protocol file, makes a `Change`, writes the file back and runs it, and `write_new` writes a new restriction-cloning file. The command line's `protocol` commands and the script functions `protocol_*` go through these, so all three surfaces change a protocol in the same way.
- **Read verification** ranks every reference by the k-mers it shares with the read, counting both strands. Ranking costs far less than aligning. Only the best-ranked reference is then aligned.

## Map layout

`render/` turns a document into a `Scene`: a list of shapes (paths, arcs, text) in pixel coordinates.
Shapes carry a tag naming the feature, enzyme or primer they draw.
`render/circular.rs` lays out plasmid maps and `render/linear.rs` lays out linear maps.
Label placement spreads crowded labels into columns and draws leader lines to their features.
A circular map can turn a chosen base to 12 o'clock (`MapOptions::top`).
`render/ring.rs` draws the small rings of the protocol bench.
A ring is a thick circle in the colours of its stretches: green where bases go on into the product, red where they are left over, and grey where nothing is decided.
A cut ring breaks into its pieces, and a piece can lift away from the centre.
A ring's radius is a fixed share of its canvas, so the app can grow a ring into a full map of the same size.

Two painters consume a `Scene`.
`render/svg.rs` writes SVG, which the script crate rasterises to PNG with resvg.
The desktop app paints the same `Scene` with egui and uses the tags for hovering and clicking.
An exported map is therefore the map the app shows.

## The command line

`slopgene-cli` defines its commands with clap in `main.rs`.
Each command module returns an `Output`, which holds a JSON value and a text rendering of the same result.
The global flags `--json` and `--jsonl` choose which one is printed.
Errors become `{"ok": false, "error": ...}` in JSON modes, so a script always receives something it can parse.

File arguments accept files, folders and glob patterns.
SlopGene expands patterns itself, because the Windows shells do not.
Commands over many files load and analyse them in parallel with rayon, and print results in input order.

## Scripting

`slopgene-script` wraps a `Document` in a `Doc` type and registers about 80 functions with a Rhai engine.
Positions in scripts are 1-based and inclusive, like the command line.
The table in `api.rs` documents every function, and `slopgene run --api` prints it.
Tests check that the table matches the functions the engine registers and that [scripting.md](scripting.md) lists every signature.

## The desktop app

`slopgene-gui` is an eframe application drawn with egui on wgpu.
It bundles its fonts (Inter and JetBrains Mono), so it looks the same on every system.

| Module | Role |
|---|---|
| `app.rs` | Tabs, views, the command palette, dialogs, the digest window, the script console |
| `protocol_view.rs` | Protocol tabs: the bench of small rings (one column per map, its cut under it, the join below), the animated zoom into a card's detailed view, and the panel that tells each step in plain words with its choices, such as the enzymes nearest a part on each side. The bench is edited by hand: a drag across a zoomed map marks the part, a click on an enzyme there picks or drops a cut, a kept piece dragged onto a plasmid or a join goes into it, and each card's menu holds the other changes. A join's card widens to show the other plasmids its tube makes in grey beside the product, and its zoom lists them in a strip with "Use this instead". A PCR has a card of its own between its plasmid and its cut, which draws the product as a bar with the copy in green and the tails in grey, and whose panel designs the primers and picks the tail enzymes. Under every zoomed picture a read-only sequence view (`seq_view`) shows the molecule's bases with the same stretch selected, and a stretch selected there shows on the picture and becomes the part |
| `align_view.rs` | The alignment view: the map with the aligned sequences, the aligned bases with features and coloured differences, and the inspector with verdicts, differences and the method |
| `doc.rs` | One open document: file, selection, undo history, cached sites, reload from disk |
| `map_view.rs` | The circular or linear map, painted from the core's `Scene` |
| `seq_view.rs` | The sequence view: both strands, ruler, features, enzyme cuts, translations, typing |
| `panels.rs` | The inspector: selection details, the feature editor, lists of features, enzymes and primers |
| `scene_paint.rs` | Painting a `Scene` with egui, and hit-testing its tagged shapes |
| `record.rs` | Playing a script of pointer and key actions into the window with a cursor drawn (`--play`), and recording the frames into an animated GIF (`--record`), for pictures of the app and the README's animations |
| `style.rs` | Fonts, colours and spacing |

Each document keeps up to 200 undo steps as snapshots of the whole `Document`.
A revision counter increases with every change, and cached analyses rebuild when it moves.
The app checks the file's modification time about once a second.
When the file changed and the tab has no unsaved edits, the tab reloads.
An agent that edits a file with `slopgene edit` therefore shows its changes in the open window.

## Packaging

The installers come from cargo-packager, configured under `[package.metadata.packager]` in `crates/slopgene-gui/Cargo.toml`.
Each installer holds both programs, with the desktop app as the main one.
`.github/workflows/release.yml` builds the installers on Windows, macOS and Linux machines when a version tag is pushed.

| File | Role |
|---|---|
| `crates/slopgene-gui/assets/icon/slopgene.svg` | The app icon, drawn as a small plasmid map |
| `crates/slopgene-gui/examples/make_icons.rs` | Renders the icon into the PNG, `.ico` and `.icns` files beside it |
| `crates/slopgene-gui/build.rs` and `slopgene-gui.rc` | Put the icon into `slopgene-gui.exe` when building on Windows |
| `crates/slopgene-gui/packaging/slopgene-gui.desktop` | The template of the Linux launcher |

The macOS app is universal.
The workflow joins the Apple silicon and Intel builds with `lipo`, signs the bundle ad hoc and makes the disk image with `hdiutil`.

## Testing

| Test | What it checks |
|---|---|
| Unit tests in each core module | Algorithms on small constructed molecules |
| `tests/corpus.rs` | Every SnapGene file found reads, rewrites byte for byte, and survives an edit |
| `tests/snapgene_parity.rs` | PCR, restriction cloning and Golden Gate products match those in SnapGene's sample project |
| `slopgene-script` tests | Script behaviour, and agreement of the API table, the engine and the documentation |
| `tests/snapgene_opens.rs` | SnapGene itself reads `.dna` files SlopGene edited, created and converted (opt-in: set `SLOPGENE_SNAPGENE_EXE`) |

The corpus and parity tests use SnapGene's installed files and the user's own files.
Those files are not in the repository.
The tests pass without checking anything when the files are absent.
Set `SLOPGENE_CORPUS` and `SLOPGENE_SNAPGENE_SAMPLES` to run them elsewhere.
