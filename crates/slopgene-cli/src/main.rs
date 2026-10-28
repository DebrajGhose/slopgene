//! `slopgene`: inspect, edit and analyse DNA, RNA and protein files from the command line.
//!
//! Every command prints readable text by default, `--json` for one JSON document, or `--jsonl`
//! for one JSON object per line (for streaming over many files). Coordinates are 1-based and
//! inclusive, like SnapGene and GenBank; on circular sequences a range may run through the
//! origin (`2600..50`). File arguments accept folders and glob patterns (`lab/**/*.dna`).

mod cloning_cmd;
mod design_cmd;
mod edit_cmd;
mod enzyme_cmd;
mod inspect;
mod protocol_cmd;
mod render_cmd;
mod update_cmd;
mod util;
mod verify_cmd;

use clap::{Args, Parser, Subcommand};
use edit_cmd::EditOp;
use slopgene_core::primers::TmConditions;
use slopgene_core::protocol::{Change, Orientation};
use util::{Mode, Output};

#[derive(Parser)]
#[command(
    name = "slopgene",
    version,
    about = "Inspect, edit and analyse sequence files (SnapGene .dna, GenBank, FASTA, JSON).",
    long_about = "Inspect, edit and analyse sequence files (SnapGene .dna, GenBank, FASTA, JSON).\n\n\
        Coordinates are 1-based and inclusive (like SnapGene and GenBank). On circular sequences a range \
        may run through the origin, for example 2600..50. File arguments accept folders and glob \
        patterns such as \"plasmids/**/*.dna\".\n\nAdd --json to any command for machine-readable output, \
        or --jsonl for one JSON object per line. Run `slopgene guide` for a guide written for scripts and AI agents."
)]
struct Cli {
    /// Print one JSON document instead of text.
    #[arg(long, global = true)]
    json: bool,
    /// Print one JSON object per line (for many files).
    #[arg(long, global = true, conflicts_with = "json")]
    jsonl: bool,
    #[command(subcommand)]
    command: Command,
}

#[derive(Args)]
struct Files {
    /// Files, folders or glob patterns.
    #[arg(required = true)]
    paths: Vec<String>,
    /// Also search subfolders of folder arguments.
    #[arg(short, long)]
    recursive: bool,
}

#[derive(Subcommand)]
enum Command {
    /// One-line summary per file: length, topology, GC, feature and primer counts.
    Info(Files),
    /// Everything in one file: metadata, features and primers (full document with --json).
    Show {
        path: String,
        /// Leave the sequence out of the JSON.
        #[arg(long)]
        no_sequence: bool,
    },
    /// Print the sequence, or part of it.
    Seq {
        path: String,
        /// Region, e.g. 100..250 (1-based, inclusive).
        #[arg(long)]
        range: Option<String>,
        /// A feature by name (or #N), read along its own strand.
        #[arg(long)]
        feature: Option<String>,
        /// Reverse-complement the result.
        #[arg(long)]
        revcomp: bool,
        /// Print as FASTA.
        #[arg(long)]
        fasta: bool,
        /// Translate instead of printing bases.
        #[arg(long)]
        translate: bool,
        /// NCBI genetic code for --translate.
        #[arg(long, default_value_t = 1)]
        table: u32,
        /// Reading frame for --translate (1-3).
        #[arg(long, default_value_t = 1)]
        frame: usize,
    },
    /// List features.
    Features {
        #[command(flatten)]
        files: Files,
        /// Only this feature type (CDS, promoter, ...).
        #[arg(long = "type")]
        kind: Option<String>,
        /// Only features whose name contains this text.
        #[arg(long)]
        name: Option<String>,
    },
    /// List primers with binding sites, Tm and warnings.
    Primers {
        #[command(flatten)]
        files: Files,
        /// Recompute binding sites instead of using the stored ones.
        #[arg(long)]
        refresh: bool,
    },
    /// Search for a DNA motif (IUPAC codes, mismatches) or a protein motif.
    Find {
        /// The motif.
        pattern: String,
        #[command(flatten)]
        files: Files,
        /// Allowed mismatches (DNA only).
        #[arg(short, long, default_value_t = 0)]
        mismatches: usize,
        /// Search translated sequence in all six frames.
        #[arg(long)]
        protein: bool,
        /// Search the top strand only.
        #[arg(long)]
        one_strand: bool,
        /// NCBI genetic code for --protein.
        #[arg(long, default_value_t = 1)]
        table: u32,
    },
    /// Find open reading frames.
    Orfs {
        path: String,
        /// Minimum length in codons (SnapGene's default is 75).
        #[arg(long, default_value_t = 75)]
        min_codons: usize,
        /// Start codons: atg, alternative or any (stop to stop).
        #[arg(long, default_value = "atg")]
        starts: String,
        /// NCBI genetic code.
        #[arg(long, default_value_t = 1)]
        table: u32,
        /// Top strand only.
        #[arg(long)]
        one_strand: bool,
        /// Include each ORF's protein in the output.
        #[arg(long)]
        protein: bool,
    },
    /// Translate a feature, a region or a sequence, with mass and pI.
    Translate {
        path: Option<String>,
        /// Translate this sequence instead of a file.
        #[arg(long)]
        seq: Option<String>,
        #[arg(long)]
        range: Option<String>,
        #[arg(long)]
        feature: Option<String>,
        #[arg(long, default_value_t = 1)]
        table: u32,
        #[arg(long, default_value_t = 1)]
        frame: usize,
        /// Stop at the first stop codon.
        #[arg(long)]
        to_stop: bool,
    },
    /// Edit a file in place (or write the result elsewhere with -o).
    Edit {
        path: String,
        /// Write the result here instead of overwriting the input.
        #[arg(short, long, global = true)]
        output: Option<String>,
        /// Report what would change without writing anything.
        #[arg(long, global = true)]
        dry_run: bool,
        #[command(subcommand)]
        op: EditCommand,
    },
    /// Convert between formats (dna, genbank, fasta, json).
    Convert {
        #[command(flatten)]
        files: Files,
        /// Target format: dna, genbank, fasta or json.
        #[arg(long)]
        to: Option<String>,
        /// Output file (one input only).
        #[arg(short, long)]
        output: Option<String>,
        /// Output folder (default: next to each input).
        #[arg(long)]
        out_dir: Option<String>,
        /// Overwrite existing files.
        #[arg(long)]
        force: bool,
    },
    /// Restriction enzymes: list, look up, sets, import from REBASE.
    Enzymes {
        #[command(subcommand)]
        action: EnzymeCommand,
    },
    /// Find restriction sites.
    Sites {
        #[command(flatten)]
        files: Files,
        /// Enzyme names (EcoRI,BamHI) or a set: common, 6+, 8+, neb, thermo, golden-gate, type-iis, blunt, all.
        #[arg(short, long, default_value = "common")]
        enzymes: String,
        /// Only enzymes that cut exactly once.
        #[arg(long)]
        unique: bool,
        #[arg(long)]
        min_cuts: Option<usize>,
        #[arg(long)]
        max_cuts: Option<usize>,
        /// Also list enzymes that do not cut.
        #[arg(long)]
        noncutters: bool,
    },
    /// Digest with restriction enzymes: fragments with exact ends, optional fragment files and gel.
    Digest {
        path: String,
        /// Enzyme names, e.g. EcoRI,BamHI.
        #[arg(short, long)]
        enzymes: String,
        /// Write each fragment as a .dna file into this folder.
        #[arg(long)]
        out_dir: Option<String>,
        /// Draw a virtual agarose gel (.svg or .png).
        #[arg(long)]
        gel: Option<String>,
        /// Ladder for the gel: 1kb, 1kb-plus, 100bp, generuler-1kb, generuler-1kb-plus.
        #[arg(long, default_value = "1kb")]
        ladder: String,
        /// Agarose percentage for the gel.
        #[arg(long, default_value_t = 1.0)]
        agarose: f64,
        /// Add an uncut-plasmid lane to the gel.
        #[arg(long)]
        uncut: bool,
    },
    /// Melting temperature and quality checks for oligos, as Primer3 computes them. The
    /// conditions default to Primer3's.
    Tm {
        /// One or more sequences.
        #[arg(required = true)]
        seqs: Vec<String>,
        /// Monovalent cations (mM).
        #[arg(long, default_value_t = 50.0)]
        na: f64,
        /// Mg2+ (mM).
        #[arg(long, default_value_t = 1.5)]
        mg: f64,
        /// dNTPs (mM).
        #[arg(long, default_value_t = 0.6)]
        dntp: f64,
        /// Oligo concentration (nM).
        #[arg(long, default_value_t = 50.0)]
        oligo: f64,
    },
    /// Simulate a PCR. Primers can be sequences or names of primers stored in the template.
    Pcr {
        template: String,
        #[arg(long)]
        fwd: String,
        #[arg(long)]
        rev: String,
        /// Write the product here (.dna, .gb, .fasta, .json).
        #[arg(short, long)]
        output: Option<String>,
        /// Minimum Tm of the annealed part.
        #[arg(long, default_value_t = 40.0)]
        min_tm: f64,
        /// Require perfect annealing (no mismatches).
        #[arg(long)]
        strict: bool,
    },
    /// Design PCR primers that amplify a region, with optional 5' tails.
    DesignPrimers {
        path: String,
        /// Region to amplify, e.g. 120..980.
        #[arg(long)]
        range: String,
        /// Target Tm of the annealing part.
        #[arg(long, default_value_t = 60.0)]
        tm: f64,
        /// 5' tail for the forward primer (restriction site, overlap ...).
        #[arg(long, default_value = "")]
        fwd_tail: String,
        /// 5' tail for the reverse primer.
        #[arg(long, default_value = "")]
        rev_tail: String,
        /// Store both primers in the file.
        #[arg(long)]
        add: bool,
        #[arg(short, long)]
        output: Option<String>,
    },
    /// Draw a plasmid or linear map (.svg, or .png).
    Render {
        path: String,
        #[arg(short, long)]
        output: Option<String>,
        #[arg(long, default_value_t = 900.0)]
        width: f32,
        #[arg(long, default_value_t = 700.0)]
        height: f32,
        /// Mark unique cutters from this set or list ("none" to hide enzymes).
        #[arg(long, default_value = "common")]
        enzymes: String,
        #[arg(long)]
        dark: bool,
        /// Pixel density for PNG output.
        #[arg(long, default_value_t = 2.0)]
        scale: f32,
        #[arg(long, default_value_t = 11.0)]
        font_size: f32,
        /// Draw a linear map even for circular sequences.
        #[arg(long)]
        linear: bool,
        /// Feature colours: palette (one calm colour per type) or file (the colours stored in the file).
        #[arg(long, default_value = "palette")]
        colors: String,
    },
    /// Run a Rhai script (see `slopgene run --api` and docs/scripting.md).
    Run {
        /// The script file.
        script: Option<String>,
        /// Run this code instead of a file.
        #[arg(short = 'e', long = "eval")]
        code: Option<String>,
        /// Values for the script's ARGS map, as key=value (repeatable).
        #[arg(long = "arg")]
        args: Vec<String>,
        /// List every scripting function and exit.
        #[arg(long)]
        api: bool,
    },
    /// Align one read to one reference and list the differences (with protein effects).
    Align {
        /// The reference map.
        reference: String,
        /// The read: a whole-plasmid result (.gbk, .fasta, .dna) or a Sanger read.
        read: String,
        /// Also show the aligned bases, the reference over the read. The JSON gets both rows.
        #[arg(long)]
        show: bool,
        /// Columns per line of aligned bases, with --show.
        #[arg(long, default_value_t = 60)]
        width: usize,
        /// Also draw the reference map with the read as a ring and its differences marked, as SVG
        /// (or PNG when the name ends in .png).
        #[arg(long, value_name = "FILE")]
        map: Option<String>,
    },
    /// Match each read to its best reference among many, align, and report the verdicts.
    Verify {
        /// Sequencing results (files, folders or patterns).
        #[arg(required = true)]
        reads: Vec<String>,
        /// Reference maps (files, folders or patterns).
        #[arg(long, required = true, num_args = 1..)]
        refs: Vec<String>,
        /// Also search subfolders of folder arguments.
        #[arg(short, long)]
        recursive: bool,
        /// How many candidate references to list per read.
        #[arg(long, default_value_t = 3)]
        top: usize,
    },
    /// Restriction cloning: cut vector and insert with the same enzymes and ligate.
    Clone {
        /// The vector to open (its largest fragment is kept as the backbone).
        vector: String,
        /// The molecule that carries the insert.
        insert: String,
        /// Enzymes, e.g. EcoRI,BamHI.
        #[arg(short, long)]
        enzymes: String,
        /// Write the product (two orientations get _1 and _2).
        #[arg(short, long)]
        output: Option<String>,
    },
    /// Homology assembly (Gibson, NEBuilder HiFi, In-Fusion) of fragments given in order.
    Assemble {
        /// The fragments in assembly order. Each one's end must overlap the next one's start.
        #[arg(required = true)]
        fragments: Vec<String>,
        /// Minimum overlap between neighbouring fragments.
        #[arg(long, default_value_t = 15)]
        min_overlap: usize,
        /// Do not close the product into a circle.
        #[arg(long)]
        linear: bool,
        /// Write the product to this file.
        #[arg(short, long)]
        output: Option<String>,
    },
    /// Golden Gate assembly: the vector first, then the parts (any order).
    GoldenGate {
        /// The destination vector, then the parts in any order.
        #[arg(required = true)]
        parts: Vec<String>,
        /// The Type IIS enzyme.
        #[arg(short, long, default_value = "BsaI")]
        enzyme: String,
        /// Write the product to this file.
        #[arg(short, long)]
        output: Option<String>,
    },
    /// Annotate files with features from a library of trusted annotated files.
    Annotate {
        #[command(flatten)]
        files: Files,
        /// Annotated files to take features from (files, folders or patterns).
        #[arg(long, required = true, num_args = 1..)]
        library: Vec<String>,
        /// Minimum identity over the whole feature, in percent.
        #[arg(long, default_value_t = 95.0)]
        min_identity: f64,
        /// Ignore library features shorter than this.
        #[arg(long, default_value_t = 30)]
        min_length: usize,
        /// Write the annotated file here (one target).
        #[arg(short, long)]
        output: Option<String>,
        /// Add the features to the targets themselves.
        #[arg(long)]
        in_place: bool,
        /// Also save the library as one annotated sequence file, for reuse.
        #[arg(long)]
        save_library: Option<String>,
    },
    /// The enzymes that cut nearest a region of a map on each side, or what a chosen pair cuts out.
    Excise {
        /// The map that carries the insert.
        path: String,
        /// The region to cut out, e.g. 1201..1950 (1-based, inclusive).
        #[arg(long, required_unless_present = "feature", conflicts_with = "feature")]
        range: Option<String>,
        /// Or a feature by name (or #N) whose bases are the region.
        #[arg(long)]
        feature: Option<String>,
        /// The vector the insert goes into: each enzyme says how often it cuts it, and a chosen pair
        /// is checked on it.
        #[arg(long)]
        vector: Option<String>,
        /// The enzymes to choose from: a set name or a list.
        #[arg(short, long, default_value = "common")]
        enzymes: String,
        /// How many enzymes to list on each side of the region.
        #[arg(long, default_value_t = slopgene_core::excise::DEFAULT_PER_SIDE)]
        count: usize,
        /// The enzyme that cuts before the region. With --right, shows the piece the pair cuts out.
        #[arg(long, requires = "right")]
        left: Option<String>,
        /// The enzyme that cuts after the region.
        #[arg(long, requires = "left")]
        right: Option<String>,
    },
    /// Cloning protocols: plan a construction step by step, run it and check every junction.
    Protocol {
        #[command(subcommand)]
        action: ProtocolCommand,
    },
    /// Open files in the desktop app.
    View {
        /// The files to open, each in its own tab.
        files: Vec<String>,
    },
    /// Print a guide for scripts and AI agents.
    Guide,
    /// Say whether a newer version of SlopGene is released on GitHub, and install it with --install.
    Update {
        /// Download the newer version and install it the way this copy was installed.
        #[arg(long)]
        install: bool,
    },
}

#[derive(Subcommand)]
enum ProtocolCommand {
    /// Run a protocol file: each step's molecules, cut sites, junctions and problems.
    Run {
        /// The protocol file (.json).
        path: String,
        /// Settle a step that waits for a choice, as STEP=OPTION (insert=2). Repeat for more steps.
        #[arg(long)]
        choose: Vec<String>,
        /// Write the choices into the protocol file.
        #[arg(long)]
        save: bool,
        /// Write the molecules each step makes into this folder.
        #[arg(long)]
        out_dir: Option<String>,
    },
    /// Write a restriction-cloning protocol for a vector and an insert, then run it.
    New {
        /// The protocol file to write (.json).
        path: String,
        /// The vector to open.
        #[arg(long)]
        vector: String,
        /// The map that carries the insert.
        #[arg(long)]
        insert: String,
        /// Enzymes that cut the insert out (and open the vector), e.g. EcoRI,BamHI. Without them the
        /// cuts wait for enzymes, chosen in the app or with `slopgene excise`.
        #[arg(short, long)]
        enzymes: Option<String>,
        /// Other enzymes for the vector, when its sites differ from the insert's (BglII for BamHI).
        #[arg(long)]
        vector_enzymes: Option<String>,
        /// The part of the insert's map to cut out: a feature's name or a range such as 389..1108.
        #[arg(long, alias = "insert-feature")]
        cut_out: Option<String>,
        /// The part of the vector the insert replaces: a feature's name or a range.
        #[arg(long)]
        replace: Option<String>,
        /// A name for the protocol.
        #[arg(long)]
        name: Option<String>,
        /// Leave the backbone's 5' phosphates on (no dephosphorylation step).
        #[arg(long)]
        keep_phosphates: bool,
        /// Copy the insert out of its map by PCR, with the enzymes' sites on the primers' tails,
        /// rather than cut it out. The cut of the PCR product then trims the tails.
        #[arg(long)]
        pcr: bool,
        /// With --pcr, the melting temperature the primer design aims for, for the part of each
        /// primer that anneals (40 to 80, 60 unless set).
        #[arg(long, requires = "pcr")]
        tm: Option<f64>,
        /// Replace the protocol file if it exists.
        #[arg(long)]
        force: bool,
    },
    /// Add a map to cut a piece from: open it, cut it around a part and keep the piece that holds
    /// the part. The file is changed and run.
    Add {
        /// The protocol file (.json).
        path: String,
        /// The map to add.
        map: String,
        /// The part to cut out: a feature's name or a range such as 389..1108.
        #[arg(long)]
        cut_out: Option<String>,
        /// Enzymes that cut it out, e.g. EcoRI,BamHI. Without them the cut waits for enzymes.
        #[arg(short, long)]
        enzymes: Option<String>,
        /// Copy the part out by PCR, with the enzymes' sites on the primers' tails, rather than
        /// cut it out.
        #[arg(long)]
        pcr: bool,
    },
    /// Change a cut: the part it is about, its enzymes or both. A vector cut with the insert's
    /// enzymes takes the new ones too. For a PCR step, the part is the one it copies and the
    /// enzymes go on its primers' tails. The file is changed and run.
    Cut {
        /// The protocol file (.json).
        path: String,
        /// The cut's step, such as cut_donor.
        step: String,
        /// The part the cut is about: a feature's name or a range such as 389..1108.
        #[arg(long)]
        part: Option<String>,
        /// Forget the part the cut is about.
        #[arg(long, conflicts_with = "part")]
        no_part: bool,
        /// The enzymes, e.g. EcoRI,BamHI. An empty list (-e "") takes them off.
        #[arg(short, long)]
        enzymes: Option<String>,
    },
    /// Keep another piece of a cut, by its number along the cut map, as a click on that piece
    /// does in the app. The file is changed and run.
    Keep {
        /// The protocol file (.json).
        path: String,
        /// The pick that keeps a piece of a cut, such as insert.
        pick: String,
        /// The piece to keep, counted from 1 as `protocol run` lists them.
        #[arg(long)]
        piece: usize,
    },
    /// Put a kept piece into a map or a join. A map is cut with the piece's enzymes and its
    /// largest piece joined to it. The file is changed and run.
    Put {
        /// The protocol file (.json).
        path: String,
        /// The step that keeps the piece, such as insert.
        piece: String,
        /// The map's open step (such as vector) or the join (such as product) that takes it.
        #[arg(long)]
        into: String,
        /// Leave the backbone's 5' phosphates on when a map is cut to take the piece.
        #[arg(long)]
        keep_phosphates: bool,
    },
    /// Take a piece out of a join. The file is changed and run.
    TakeOut {
        /// The protocol file (.json).
        path: String,
        /// The piece the join takes, such as insert_2.
        piece: String,
        /// The join, such as product.
        #[arg(long)]
        join: String,
    },
    /// Take a map out of the protocol with every step that works on it. The file is changed and
    /// run.
    RemoveMap {
        /// The protocol file (.json).
        path: String,
        /// The map's open step, such as donor_2.
        step: String,
    },
    /// Leave a step with one input out, such as a dephosphorylation, and pass its input on. The
    /// file is changed and run.
    LeaveOut {
        /// The protocol file (.json).
        path: String,
        /// The step to leave out, such as cip.
        step: String,
    },
    /// Take the 5' phosphates off a piece on its way into a join, so it cannot close on itself.
    /// The file is changed and run.
    Dephosphorylate {
        /// The protocol file (.json).
        path: String,
        /// The piece the join takes, such as backbone.
        piece: String,
        /// The join, such as product.
        #[arg(long)]
        join: String,
    },
    /// Say which way a piece goes into a join, as "Use this instead" on one of the join's other
    /// plasmids does in the app. The file is changed and run.
    Turn {
        /// The protocol file (.json).
        path: String,
        /// The piece the join takes, such as insert. The first piece sets the direction.
        piece: String,
        /// The join, such as product.
        #[arg(long)]
        join: String,
        /// forward (as given), reverse (turned round) or either (the join then waits for a choice
        /// when both ways fit).
        #[arg(long)]
        way: String,
    },
    /// Copy a map's part out by PCR instead of cutting it out: a PCR step goes between the map
    /// and its cut, with the cut's enzymes on the primers' tails, and the cut then trims the
    /// tails. The file is changed and run.
    CopyOut {
        /// The protocol file (.json).
        path: String,
        /// Any step of the map's lane: the map (such as donor), its cut or its pick.
        step: String,
    },
    /// Cut a map's part out instead of copying it by PCR, the reverse of copy-out. The file is
    /// changed and run.
    CutOut {
        /// The protocol file (.json).
        path: String,
        /// Any step of the map's lane: the map (such as donor), its PCR, its cut or its pick.
        step: String,
    },
    /// The primers of a PCR step: keep one of the pairs the run designs, use two of your own, or
    /// design them anew on each run. The file is changed and run.
    Primers {
        /// The protocol file (.json).
        path: String,
        /// The PCR step, such as pcr_donor.
        step: String,
        /// Keep this pair of those `protocol run` lists, counted from 1 (1 is the best).
        #[arg(long, conflicts_with_all = ["fwd", "rev", "clear"])]
        pair: Option<usize>,
        /// The forward primer, 5' to 3' with its tail.
        #[arg(long, requires = "rev", conflicts_with = "clear")]
        fwd: Option<String>,
        /// The reverse primer, 5' to 3' with its tail.
        #[arg(long, requires = "fwd")]
        rev: Option<String>,
        /// Forget the kept primers, so the run designs them anew.
        #[arg(long)]
        clear: bool,
        /// The melting temperature the design aims for, for the part of each primer that anneals
        /// (40 to 80, 60 unless set). The primers are designed anew at it, unless --pair keeps one.
        #[arg(long, conflicts_with_all = ["fwd", "rev"])]
        tm: Option<f64>,
    },
    /// Draw the protocol as the app's bench: each map as a small ring, its cut or PCR under it with
    /// the pieces coloured, and the join with the new plasmid, as SVG or PNG.
    Render {
        /// The protocol file (.json).
        path: String,
        /// The picture to write (.svg or .png); by default the protocol file's name with .svg.
        #[arg(short, long)]
        output: Option<String>,
        #[arg(long)]
        dark: bool,
        /// Pixel density for PNG output.
        #[arg(long, default_value_t = 2.0)]
        scale: f32,
    },
}

#[derive(Subcommand)]
enum EditCommand {
    /// Insert bases so they start at position --at (use length+1 to append).
    Insert {
        #[arg(long)]
        at: usize,
        #[arg(long)]
        seq: String,
    },
    /// Delete a range.
    Delete {
        #[arg(long)]
        range: String,
    },
    /// Replace a range with new bases (same length = point mutation, nothing moves).
    Replace {
        #[arg(long)]
        range: String,
        #[arg(long)]
        seq: String,
    },
    /// Reverse-complement the whole molecule.
    Revcomp,
    /// Reverse-complement a region in place (flip an insert).
    Invert {
        #[arg(long)]
        range: String,
    },
    /// Make base --at the new base 1 (circular only).
    Origin {
        #[arg(long)]
        at: usize,
    },
    /// Make the sequence linear or circular.
    Topology { topology: String },
    /// Set the description.
    Describe { text: String },
    /// Add a feature.
    AddFeature {
        #[arg(long)]
        name: String,
        #[arg(long = "type", default_value = "misc_feature")]
        kind: String,
        #[arg(long)]
        range: String,
        #[arg(long, default_value = "forward")]
        strand: String,
        #[arg(long)]
        color: Option<String>,
        /// Qualifier as key=value (repeatable).
        #[arg(long = "qualifier")]
        qualifiers: Vec<String>,
    },
    /// Remove a feature by name (or #N).
    RemoveFeature {
        name: String,
        /// Remove every feature with this name.
        #[arg(long)]
        all: bool,
    },
    /// Change a feature's name, type, range, strand or colour.
    UpdateFeature {
        name: String,
        #[arg(long)]
        new_name: Option<String>,
        #[arg(long = "type")]
        kind: Option<String>,
        #[arg(long)]
        range: Option<String>,
        #[arg(long)]
        strand: Option<String>,
        #[arg(long)]
        color: Option<String>,
    },
    /// Add a primer; its binding sites are found automatically.
    AddPrimer {
        #[arg(long)]
        name: String,
        #[arg(long)]
        seq: String,
        #[arg(long)]
        description: Option<String>,
    },
    /// Remove a primer by name.
    RemovePrimer { name: String },
    /// Recompute every primer's binding sites.
    RefreshPrimers,
}

#[derive(Subcommand)]
enum EnzymeCommand {
    /// List enzymes (optionally one set, or matching a name or site).
    List {
        #[arg(long)]
        set: Option<String>,
        #[arg(long)]
        search: Option<String>,
    },
    /// Details for named enzymes.
    Info {
        #[arg(required = true)]
        names: Vec<String>,
    },
    /// The named enzyme sets.
    Sets,
    /// Build an enzyme table from REBASE's EMBOSS files (use it via SLOPGENE_ENZYMES).
    ImportRebase {
        #[arg(long)]
        emboss_e: String,
        #[arg(long)]
        emboss_r: String,
        #[arg(short, long)]
        output: String,
    },
}

const GUIDE: &str = include_str!("../../../docs/agent-guide.md");

fn run(command: Command) -> anyhow::Result<Output> {
    match command {
        Command::Info(f) => inspect::info(&f.paths, f.recursive),
        Command::Show { path, no_sequence } => inspect::show(&path, !no_sequence),
        Command::Seq { path, range, feature, revcomp, fasta, translate, table, frame } => {
            inspect::seq(&path, range.as_deref(), feature.as_deref(), revcomp, fasta, translate, table, frame)
        }
        Command::Features { files, kind, name } => inspect::features(&files.paths, files.recursive, kind.as_deref(), name.as_deref()),
        Command::Primers { files, refresh } => inspect::primers(&files.paths, files.recursive, refresh),
        Command::Find { pattern, files, mismatches, protein, one_strand, table } => {
            inspect::find(&files.paths, files.recursive, &pattern, mismatches, protein, one_strand, table)
        }
        Command::Orfs { path, min_codons, starts, table, one_strand, protein } => inspect::orfs(&path, min_codons, &starts, table, one_strand, protein),
        Command::Translate { path, seq, range, feature, table, frame, to_stop } => {
            inspect::translate_cmd(path.as_deref(), seq.as_deref(), range.as_deref(), feature.as_deref(), table, frame, to_stop)
        }
        Command::Edit { path, output, dry_run, op } => {
            let op = match op {
                EditCommand::Insert { at, seq } => EditOp::Insert { at, seq },
                EditCommand::Delete { range } => EditOp::Delete { range },
                EditCommand::Replace { range, seq } => EditOp::Replace { range, seq },
                EditCommand::Revcomp => EditOp::Revcomp,
                EditCommand::Invert { range } => EditOp::Invert { range },
                EditCommand::Origin { at } => EditOp::Origin { at },
                EditCommand::Topology { topology } => EditOp::Topology { topology },
                EditCommand::Describe { text } => EditOp::Describe { text },
                EditCommand::AddFeature { name, kind, range, strand, color, qualifiers } => EditOp::AddFeature { name, kind, range, strand, color, qualifiers },
                EditCommand::RemoveFeature { name, all } => EditOp::RemoveFeature { name, all },
                EditCommand::UpdateFeature { name, new_name, kind, range, strand, color } => {
                    EditOp::UpdateFeature { name, new_name, kind, range, strand, color }
                }
                EditCommand::AddPrimer { name, seq, description } => EditOp::AddPrimer { name, seq, description },
                EditCommand::RemovePrimer { name } => EditOp::RemovePrimer { name },
                EditCommand::RefreshPrimers => EditOp::RefreshPrimers,
            };
            edit_cmd::run(&path, op, output.as_deref(), dry_run)
        }
        Command::Convert { files, to, output, out_dir, force } => {
            design_cmd::convert(&files.paths, files.recursive, to.as_deref(), output.as_deref(), out_dir.as_deref(), force)
        }
        Command::Enzymes { action } => match action {
            EnzymeCommand::List { set, search } => enzyme_cmd::list(set.as_deref(), search.as_deref()),
            EnzymeCommand::Info { names } => enzyme_cmd::info(&names),
            EnzymeCommand::Sets => enzyme_cmd::sets(),
            EnzymeCommand::ImportRebase { emboss_e, emboss_r, output } => enzyme_cmd::import_rebase(&emboss_e, &emboss_r, &output),
        },
        Command::Sites { files, enzymes, unique, min_cuts, max_cuts, noncutters } => {
            enzyme_cmd::sites(&files.paths, files.recursive, &enzymes, unique, min_cuts, max_cuts, noncutters)
        }
        Command::Digest { path, enzymes, out_dir, gel, ladder, agarose, uncut } => {
            enzyme_cmd::digest_cmd(&path, &enzymes, out_dir.as_deref(), gel.as_deref(), &ladder, agarose, uncut)
        }
        Command::Tm { seqs, na, mg, dntp, oligo } => {
            design_cmd::tm(&seqs, TmConditions { monovalent_mm: na, magnesium_mm: mg, dntp_mm: dntp, oligo_nm: oligo })
        }
        Command::Pcr { template, fwd, rev, output, min_tm, strict } => design_cmd::pcr_cmd(&template, &fwd, &rev, output.as_deref(), min_tm, strict),
        Command::DesignPrimers { path, range, tm, fwd_tail, rev_tail, add, output } => {
            design_cmd::design_primers(&path, &range, tm, &fwd_tail, &rev_tail, add, output.as_deref())
        }
        Command::Render { path, output, width, height, enzymes, dark, scale, font_size, linear, colors } => {
            render_cmd::run(&path, output.as_deref(), width, height, &enzymes, dark, scale, font_size, linear, &colors)
        }
        Command::Guide => Ok(Output::new(serde_json::json!({ "guide": GUIDE }), GUIDE)),
        Command::Update { install } => update_cmd::update(install),
        Command::Align { reference, read, show, width, map } => verify_cmd::align(&reference, &read, show.then_some(width), map.as_deref()),
        Command::Clone { vector, insert, enzymes, output } => cloning_cmd::clone_cmd(&vector, &insert, &enzymes, output.as_deref()),
        Command::Assemble { fragments, min_overlap, linear, output } => cloning_cmd::assemble(&fragments, min_overlap, linear, output.as_deref()),
        Command::GoldenGate { parts, enzyme, output } => cloning_cmd::golden_gate_cmd(&parts, &enzyme, output.as_deref()),
        Command::View { files } => cloning_cmd::view(&files),
        Command::Excise { path, range, feature, vector, enzymes, count, left, right } => cloning_cmd::excise(&cloning_cmd::ExcisePlan {
            path: &path,
            range: range.as_deref(),
            feature: feature.as_deref(),
            vector: vector.as_deref(),
            enzymes: &enzymes,
            count,
            left: left.as_deref(),
            right: right.as_deref(),
        }),
        Command::Protocol { action } => match action {
            ProtocolCommand::Run { path, choose, save, out_dir } => protocol_cmd::run(&path, &choose, save, out_dir.as_deref()),
            ProtocolCommand::New { path, vector, insert, enzymes, vector_enzymes, cut_out, replace, name, keep_phosphates, pcr, tm, force } => {
                protocol_cmd::new_cloning(&protocol_cmd::CloningPlan {
                    path: &path,
                    tm,
                    vector: &vector,
                    insert: &insert,
                    enzymes: enzymes.as_deref(),
                    vector_enzymes: vector_enzymes.as_deref(),
                    cut_out: cut_out.as_deref(),
                    replace: replace.as_deref(),
                    name: name.as_deref(),
                    keep_phosphates,
                    pcr,
                    force,
                })
            }
            // The changes the app's bench makes, each a `protocol::Change` made to the file.
            ProtocolCommand::Add { path, map, cut_out, enzymes, pcr } => {
                let enzymes = protocol_cmd::enzyme_names(enzymes.as_deref().unwrap_or(""))?;
                protocol_cmd::change(&path, Change::AddMap { map: map.into(), cut_out, enzymes, pcr })
            }
            ProtocolCommand::Cut { path, step, part, no_part, enzymes } => {
                let part = if no_part { Some(None) } else { part.map(Some) };
                let enzymes = enzymes.as_deref().map(protocol_cmd::enzyme_names).transpose()?;
                protocol_cmd::change(&path, Change::Cut { step, part, enzymes })
            }
            ProtocolCommand::Keep { path, pick, piece } => protocol_cmd::change(&path, Change::Keep { pick, piece }),
            ProtocolCommand::Put { path, piece, into, keep_phosphates } => {
                protocol_cmd::change(&path, Change::Put { piece, into, dephosphorylate: !keep_phosphates })
            }
            ProtocolCommand::TakeOut { path, piece, join } => protocol_cmd::change(&path, Change::TakeOut { piece, join }),
            ProtocolCommand::RemoveMap { path, step } => protocol_cmd::change(&path, Change::RemoveMap { step }),
            ProtocolCommand::LeaveOut { path, step } => protocol_cmd::change(&path, Change::LeaveOut { step }),
            ProtocolCommand::Dephosphorylate { path, piece, join } => protocol_cmd::change(&path, Change::Dephosphorylate { piece, join }),
            ProtocolCommand::Turn { path, piece, join, way } => {
                let way = Orientation::parse(&way).ok_or_else(|| anyhow::anyhow!("--way takes forward, reverse or either (got '{way}')"))?;
                protocol_cmd::change(&path, Change::Turn { piece, join, way })
            }
            ProtocolCommand::CopyOut { path, step } => protocol_cmd::change(&path, Change::CopyOut { step }),
            ProtocolCommand::CutOut { path, step } => protocol_cmd::change(&path, Change::CutOut { step }),
            ProtocolCommand::Primers { path, step, pair, fwd, rev, clear, tm } => {
                if pair.is_none() && fwd.is_none() && !clear && tm.is_none() {
                    anyhow::bail!("say what to do with the primers: --pair N, --fwd and --rev, --tm T, or --clear");
                }
                protocol_cmd::change(&path, Change::Primers { step, pair, forward: fwd, reverse: rev, tm })
            }
            ProtocolCommand::Render { path, output, dark, scale } => protocol_cmd::render(&path, output.as_deref(), dark, scale),
        },
        Command::Annotate { files, library, min_identity, min_length, output, in_place, save_library } => {
            cloning_cmd::annotate(&files.paths, files.recursive, &library, min_identity, min_length, output.as_deref(), in_place, save_library.as_deref())
        }
        Command::Verify { reads, refs, recursive, top } => verify_cmd::verify(&reads, &refs, recursive, top),
        Command::Run { .. } => unreachable!("handled in main"),
    }
}

/// Run a script. Text mode streams `print` output as it happens; JSON modes collect it and
/// report it with the script's final value.
fn run_script(script: Option<String>, code: Option<String>, args: Vec<String>, api: bool, mode: Mode) -> anyhow::Result<Option<Output>> {
    if api {
        return Ok(Some(Output::new(slopgene_script::api_json(), slopgene_script::api_text())));
    }
    let source = match (code, script) {
        (Some(c), _) => c,
        (None, Some(path)) => std::fs::read_to_string(&path).map_err(|e| anyhow::anyhow!("{path}: {e}"))?,
        (None, None) => anyhow::bail!("give a script file or -e CODE (or --api for the function list)"),
    };
    let mut pairs = Vec::new();
    for a in args {
        let (k, v) = a.split_once('=').ok_or_else(|| anyhow::anyhow!("--arg must be key=value (got '{a}')"))?;
        pairs.push((k.to_string(), v.to_string()));
    }
    let buffer = (mode != Mode::Text).then(|| std::sync::Arc::new(std::sync::Mutex::new(String::new())));
    let options = slopgene_script::RunOptions { args: pairs, output: buffer.clone(), document: None };
    let result = slopgene_script::run(&source, options).map_err(|e| anyhow::anyhow!(e))?;
    match buffer {
        None => {
            if !result.is_unit() {
                slopgene_script::stdout_line(&result.to_string());
            }
            Ok(None)
        }
        Some(b) => {
            let printed = b.lock().map(|s| s.clone()).unwrap_or_default();
            let value: serde_json::Value = slopgene_script::to_json_value(&result);
            Ok(Some(Output::new(serde_json::json!({ "ok": true, "output": printed, "result": value }), printed)))
        }
    }
}

fn main() {
    let cli = Cli::parse();
    let mode = if cli.json {
        Mode::Json
    } else if cli.jsonl {
        Mode::JsonLines
    } else {
        Mode::Text
    };
    let outcome = match cli.command {
        Command::Run { script, code, args, api } => run_script(script, code, args, api, mode),
        other => run(other).map(Some),
    };
    match outcome {
        Ok(Some(output)) => output.print(mode),
        Ok(None) => {}
        Err(e) => {
            // Errors are JSON too in JSON modes, so a script always gets something to parse.
            let message = format!("{e:#}");
            match mode {
                Mode::Text => eprintln!("error: {message}"),
                _ => slopgene_script::stdout_line(&serde_json::json!({ "ok": false, "error": message }).to_string()),
            }
            std::process::exit(1);
        }
    }
}
