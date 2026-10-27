//! The scripting API as a script writer sees it: every function and property, one line each.
//!
//! `slopgene run --api` prints this table, and docs/scripting.md lists the same signatures.
//! The tests below check that the table names exactly the functions and properties the engine
//! registers, and that the documentation contains every signature, so the three cannot drift.

/// One documented function or property.
pub struct ApiEntry {
    /// Section heading, such as "Files and folders".
    pub section: &'static str,
    /// The call as a script writes it, such as `open(path) -> Doc` or `doc.subseq(start, end)`.
    pub signature: &'static str,
    /// What it does, in one sentence.
    pub summary: &'static str,
}

const fn entry(section: &'static str, signature: &'static str, summary: &'static str) -> ApiEntry {
    ApiEntry { section, signature, summary }
}

const FILES: &str = "Files and folders";
const OPENING: &str = "Opening and creating documents";
const PROPERTIES: &str = "Document properties";
const READING: &str = "Reading parts of a document";
const EDITING: &str = "Editing in place";
const SAVING: &str = "Saving and drawing";
const ENZYMES: &str = "Restriction enzymes";
const ANALYSIS: &str = "Analysis";
const CLONING: &str = "Cloning";
const PROTOCOLS: &str = "Cloning protocols";
const SEQUENCES: &str = "Sequences and oligos";
const UTILITIES: &str = "Utilities";

/// Every function and property, in the order the documentation presents them.
pub const API: &[ApiEntry] = &[
    entry(FILES, "glob(pattern) -> [string]", "Paths that match a pattern such as \"maps/**/*.dna\", sorted."),
    entry(FILES, "files(dir) -> [string]", "Sequence files in a folder, sorted."),
    entry(FILES, "files(dir, recursive) -> [string]", "Sequence files in a folder and, when recursive is true, in its subfolders."),
    entry(FILES, "exists(path) -> bool", "Whether a file or folder exists."),
    entry(FILES, "file_name(path) -> string", "The last part of a path, such as \"pUC19.dna\"."),
    entry(FILES, "file_stem(path) -> string", "The file name without its extension, such as \"pUC19\"."),
    entry(FILES, "read_text(path) -> string", "The contents of a text file."),
    entry(FILES, "write_text(path, text)", "Write text to a file, replacing the file in one step."),
    entry(FILES, "mkdir(path)", "Create a folder and any missing parent folders."),
    entry(OPENING, "open(path) -> Doc", "Open a SnapGene, GenBank, FASTA or JSON file (the first record when it holds several)."),
    entry(OPENING, "open_all(path) -> [Doc]", "Every record in a file."),
    entry(OPENING, "new_dna(name, seq) -> Doc", "A new linear DNA document."),
    entry(OPENING, "new_dna(name, seq, topology) -> Doc", "A new DNA document that is \"linear\" or \"circular\"."),
    entry(OPENING, "from_json(text) -> Doc", "A document from SlopGene document JSON, the format `to_json` writes."),
    entry(OPENING, "doc.copy() -> Doc", "A copy of the document (assigning a document to a new variable also copies it)."),
    entry(PROPERTIES, "doc.name", "The name, which can be assigned."),
    entry(PROPERTIES, "doc.description", "The description, which can be assigned."),
    entry(PROPERTIES, "doc.length", "The number of bases or residues."),
    entry(PROPERTIES, "doc.len() -> int", "The number of bases or residues, the same as `doc.length`."),
    entry(PROPERTIES, "doc.seq", "The sequence as a string."),
    entry(PROPERTIES, "doc.topology", "\"linear\" or \"circular\". Assigning either value changes the topology."),
    entry(PROPERTIES, "doc.molecule", "\"dna\", \"rna\" or \"protein\"."),
    entry(PROPERTIES, "doc.is_circular", "Whether the molecule is circular."),
    entry(PROPERTIES, "doc.gc", "GC content in percent."),
    entry(PROPERTIES, "doc.path", "The file the document was opened from or last saved to, or \"\" for a new document."),
    entry(PROPERTIES, "doc.features", "Every feature as a map with name, type, strand, start, end, length, color, segments, qualifiers and index."),
    entry(PROPERTIES, "doc.primers", "Every primer as a map with name, sequence, length, description and sites."),
    entry(PROPERTIES, "doc.to_map() -> map", "The whole document as a map with the content of `to_json`."),
    entry(PROPERTIES, "doc.to_json() -> string", "The whole document as SlopGene document JSON."),
    entry(PROPERTIES, "doc.to_string() -> string", "The name, length and topology on one line, which is what `print(doc)` shows."),
    entry(READING, "doc.subseq(start, end) -> string", "The bases from start to end, through the origin when start > end on a circular molecule."),
    entry(READING, "doc.feature(name) -> map", "One feature, by name in any case or by \"#N\" for the Nth feature."),
    entry(READING, "doc.has_feature(name) -> bool", "Whether any feature has this name."),
    entry(READING, "doc.features_of_type(type) -> [map]", "The features of one type, such as \"CDS\"."),
    entry(READING, "doc.feature_seq(name) -> string", "A feature's bases, read 5' to 3' along its strand."),
    entry(READING, "doc.translate_feature(name) -> string", "The protein a CDS feature encodes."),
    entry(READING, "doc.extract(start, end) -> Doc", "A region as a new linear document, with the features that overlap it clipped to the region."),
    entry(EDITING, "doc.insert(pos, bases)", "Insert bases so that the first new base is at pos, from 1 to length + 1."),
    entry(EDITING, "doc.delete(start, end)", "Delete the bases from start to end."),
    entry(EDITING, "doc.replace(start, end, bases)", "Replace the bases from start to end with new bases."),
    entry(EDITING, "doc.invert(start, end)", "Reverse-complement a region in place, flipping the features wholly inside it."),
    entry(EDITING, "doc.reverse_complement()", "Reverse-complement the whole molecule, flipping every feature and primer with it."),
    entry(EDITING, "doc.set_origin(pos)", "Make pos the first base of a circular molecule."),
    entry(EDITING, "doc.add_feature(name, type, start, end, strand) -> map", "Add a feature whose strand is \"forward\", \"reverse\", \"both\" or \"none\"."),
    entry(EDITING, "doc.add_feature(spec) -> map", "Add a feature from a map with name, type, start, end, strand, color and qualifiers."),
    entry(EDITING, "doc.update_feature(name, changes)", "Change name, type, strand, color or range (\"100..250\"), and set any other key as a qualifier."),
    entry(EDITING, "doc.remove_feature(name) -> int", "Remove every feature with this name and return how many were removed."),
    entry(EDITING, "doc.add_primer(name, seq) -> map", "Add a primer and find where it binds."),
    entry(EDITING, "doc.remove_primer(name) -> int", "Remove every primer with this name and return how many were removed."),
    entry(EDITING, "doc.refresh_primers()", "Recompute where every primer binds."),
    entry(SAVING, "doc.save() -> string", "Write the document back to the file it came from and return the path."),
    entry(SAVING, "doc.save(path) -> string", "Write to a file in the format its extension names (.dna, .gbk, .fasta or .json)."),
    entry(SAVING, "doc.save(path, format) -> string", "Write in a named format: \"dna\", \"genbank\", \"fasta\" or \"json\"."),
    entry(SAVING, "doc.render(path) -> string", "Draw the map as SVG, or as PNG when the path ends in .png."),
    entry(SAVING, "doc.render(path, options) -> string", "Draw the map with options for width, height, enzymes, dark and colors."),
    entry(ENZYMES, "doc.sites(enzymes) -> [map]", "Every cutting site as a map with enzyme, start, end, strand and cut, in position order."),
    entry(ENZYMES, "doc.cut_counts(enzymes) -> map", "How many times each enzyme cuts, keyed by enzyme name."),
    entry(ENZYMES, "doc.unique_cutters(enzymes) -> [map]", "The enzymes that cut exactly once, as maps with enzyme and cut."),
    entry(
        ENZYMES,
        "doc.digest(enzymes) -> [map]",
        "The fragments, each with length, start, end, the enzyme and end type on each side, and the fragment itself as doc.",
    ),
    entry(ENZYMES, "enzyme(name) -> map", "One enzyme with its name, site, site_with_cuts and suppliers."),
    entry(ENZYMES, "enzymes(spec) -> [string]", "The enzyme names of a set or of a comma-separated list."),
    entry(ANALYSIS, "doc.find(pattern) -> [map]", "Matches of a DNA motif on both strands, as maps with start, end, strand and mismatches."),
    entry(ANALYSIS, "doc.find(pattern, mismatches) -> [map]", "Matches of a DNA motif with up to this many mismatched bases."),
    entry(ANALYSIS, "doc.find_protein(motif) -> [map]", "Matches of a protein motif in all six reading frames."),
    entry(ANALYSIS, "doc.orfs() -> [map]", "Open reading frames of at least 75 codons, as maps with start, end, strand, frame and codons."),
    entry(ANALYSIS, "doc.orfs(min_codons) -> [map]", "Open reading frames of at least this many codons."),
    entry(
        ANALYSIS,
        "align(reference, read) -> map",
        "The alignment of a read to a reference: verdict, identity, coverage, structural changes and each difference.",
    ),
    entry(ANALYSIS, "doc.annotate(library) -> [map]", "Add the library features found at 95% identity or more, and return the added features."),
    entry(ANALYSIS, "doc.annotate(library, min_identity) -> [map]", "Add the library features found at this identity in percent or more."),
    entry(
        CLONING,
        "doc.nearest_cuts(start, end, enzymes) -> map",
        "The enzymes that cut nearest the region on each side, nearest first (before and after), and those that cut inside it.",
    ),
    entry(
        CLONING,
        "doc.nearest_cuts(start, end, enzymes, options) -> map",
        "The same, with options vector (a document, so each enzyme says how often it cuts it) and count (enzymes on each side).",
    ),
    entry(
        CLONING,
        "doc.excision(start, end, left, right) -> map",
        "What cutting with left before the region and right after it gives: the piece that holds the region, its ends and direction.",
    ),
    entry(CLONING, "doc.excision(start, end, left, right, options) -> map", "The same, with the option vector, on which the pair is checked too."),
    entry(CLONING, "doc.pcr(fwd, rev) -> Doc", "The PCR product of two primers, given as sequences or as names of primers in the document."),
    entry(CLONING, "ligate(left, right) -> Doc", "Join two linear fragments whose facing ends are compatible."),
    entry(CLONING, "circularize(fragment) -> Doc", "Close a linear fragment whose two ends are compatible."),
    entry(
        CLONING,
        "restriction_clone(vector, insert, enzymes) -> [Doc]",
        "Cut both molecules with the same enzymes and ligate, giving one product per orientation that fits.",
    ),
    entry(CLONING, "assemble(fragments) -> Doc", "Join fragments in order through end overlaps of at least 15 bp and close the circle."),
    entry(CLONING, "assemble(fragments, min_overlap) -> Doc", "Join fragments in order through end overlaps of at least this length."),
    entry(CLONING, "golden_gate(parts, enzyme) -> Doc", "Golden Gate assembly with a type IIS enzyme, the vector first and the parts in any order."),
    entry(PROTOCOLS, "protocol_run(path) -> map", "Run a protocol file and report every step, as `slopgene protocol run --json` does."),
    entry(
        PROTOCOLS,
        "protocol_new(path, vector, insert) -> map",
        "Write a protocol that cuts a piece out of the insert's map and puts it into the vector, then run it.",
    ),
    entry(
        PROTOCOLS,
        "protocol_new(path, vector, insert, options) -> map",
        "The same, with the options enzymes, vector_enzymes, cut_out, replace, name, keep_phosphates, pcr (copy the insert out by PCR), tm (the melting temperature its primer design aims for) and force.",
    ),
    entry(PROTOCOLS, "protocol_add(path, map) -> map", "Add a map to cut a piece from. Its cut waits for a part and enzymes."),
    entry(PROTOCOLS, "protocol_add(path, map, options) -> map", "The same, with the options cut_out (the part), enzymes and pcr (copy the part out by PCR)."),
    entry(
        PROTOCOLS,
        "protocol_cut(path, step, options) -> map",
        "Change a cut's part, its enzymes or both, with the options part and enzymes (\"\" clears either). For a PCR step, the part it copies and the enzymes on its primers' tails.",
    ),
    entry(PROTOCOLS, "protocol_keep(path, pick, piece) -> map", "Make a pick keep another piece of its cut, counted from 1 along the cut map."),
    entry(PROTOCOLS, "protocol_put(path, piece, into) -> map", "Put a kept piece into a map, which is cut with the piece's enzymes, or into a join."),
    entry(PROTOCOLS, "protocol_take_out(path, piece, join) -> map", "Take a piece out of a join."),
    entry(PROTOCOLS, "protocol_remove_map(path, step) -> map", "Take a map out of the protocol with every step that works on it."),
    entry(PROTOCOLS, "protocol_leave_out(path, step) -> map", "Leave out a step that has one input, and pass that input on."),
    entry(PROTOCOLS, "protocol_dephosphorylate(path, piece, join) -> map", "Take the 5' phosphates off a piece on its way into a join."),
    entry(PROTOCOLS, "protocol_turn(path, piece, join, way) -> map", "Say which way a piece goes into a join: \"forward\", \"reverse\" or \"either\"."),
    entry(
        PROTOCOLS,
        "protocol_copy_out(path, step) -> map",
        "Copy a map's part out by PCR instead of cutting it out: a PCR step goes between the map and its cut, with the cut's enzymes on the primers' tails.",
    ),
    entry(PROTOCOLS, "protocol_cut_out(path, step) -> map", "Cut a map's part out instead of copying it by PCR, the reverse of protocol_copy_out."),
    entry(
        PROTOCOLS,
        "protocol_render(path, out) -> map",
        "Draw the protocol as the app's bench, each map a small ring with its cut or PCR under it and the join below, and write it as SVG or PNG by the extension of out.",
    ),
    entry(PROTOCOLS, "protocol_render(path, out, options) -> map", "The same, with the options dark and scale (the pixel density of a PNG)."),
    entry(
        PROTOCOLS,
        "protocol_primers(path, step, options) -> map",
        "The primers of a PCR step: pair keeps pair N of those the run designs, forward and reverse give your own, tm aims the design at a melting temperature, and clear: true designs them anew.",
    ),
    entry(SEQUENCES, "revcomp(seq) -> string", "The reverse complement."),
    entry(SEQUENCES, "complement(seq) -> string", "The complement, in the same direction."),
    entry(SEQUENCES, "gc(seq) -> float", "GC content in percent."),
    entry(SEQUENCES, "translate(seq) -> string", "The protein from the first base, with the standard genetic code."),
    entry(SEQUENCES, "translate(seq, table) -> string", "The protein from the first base, with the NCBI genetic code of this table number."),
    entry(SEQUENCES, "tm(seq) -> float", "The nearest-neighbour melting temperature at SnapGene's default conditions."),
    entry(SEQUENCES, "primer_report(seq) -> map", "Length, tm, gc_percent, self_dimer, three_prime_dimer, hairpin_stem and warnings for an oligo."),
    entry(UTILITIES, "today() -> string", "Today's date as YYYY-MM-DD."),
    entry(UTILITIES, "join(items, separator) -> string", "The items of an array joined into one string."),
    entry(UTILITIES, "to_json(value) -> string", "Any value as JSON text."),
    entry(UTILITIES, "parse_json(text) -> any", "JSON text as a map, an array or a simple value."),
];

/// Width of the signature column in [`api_text`]. Longer signatures put their summary on the
/// next line, so the listing stays readable in an 80-to-120-column terminal.
const SIGNATURE_COLUMN: usize = 40;

/// The table as text, one section after another, for `slopgene run --api`.
pub fn api_text() -> String {
    let mut text = String::new();
    let mut section = "";
    for e in API {
        if e.section != section {
            if !text.is_empty() {
                text.push('\n');
            }
            text.push_str(e.section);
            text.push('\n');
            section = e.section;
        }
        if e.signature.len() <= SIGNATURE_COLUMN {
            text.push_str(&format!("  {:<SIGNATURE_COLUMN$}  {}\n", e.signature, e.summary));
        } else {
            text.push_str(&format!("  {}\n  {:<SIGNATURE_COLUMN$}  {}\n", e.signature, "", e.summary));
        }
    }
    text
}

/// The table as JSON, for `slopgene run --api --json`.
pub fn api_json() -> serde_json::Value {
    API.iter().map(|e| serde_json::json!({ "section": e.section, "signature": e.signature, "summary": e.summary })).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    /// Registered functions that are deliberately not documented (Rhai calls `to_debug` for
    /// `debug(doc)`, and it says the same as `to_string`).
    const UNDOCUMENTED: &[&str] = &["to_debug"];

    /// The function or property name of a signature, and whether it is a property.
    fn name_of(signature: &str) -> (&str, bool) {
        let call = signature.strip_prefix("doc.").unwrap_or(signature);
        match call.find('(') {
            Some(paren) => (&call[..paren], false),
            None => (call.split_whitespace().next().unwrap_or(call), true),
        }
    }

    /// Function names and property names the engine registers.
    fn registered() -> (BTreeSet<String>, BTreeSet<String>) {
        let mut functions = BTreeSet::new();
        let mut properties = BTreeSet::new();
        for sig in crate::engine(None).gen_fn_signatures(false) {
            let name = sig.split('(').next().unwrap_or("").to_string();
            if let Some(p) = name.strip_prefix("get$") {
                properties.insert(p.to_string());
            } else if !name.starts_with("set$") && !UNDOCUMENTED.contains(&name.as_str()) {
                functions.insert(name);
            }
        }
        (functions, properties)
    }

    #[test]
    fn table_lists_exactly_what_the_engine_registers() {
        let (functions, properties) = registered();
        let mut documented_functions = BTreeSet::new();
        let mut documented_properties = BTreeSet::new();
        for e in API {
            let (name, is_property) = name_of(e.signature);
            if is_property {
                documented_properties.insert(name.to_string());
            } else {
                documented_functions.insert(name.to_string());
            }
        }
        assert_eq!(documented_functions, functions, "functions in the API table versus the engine");
        assert_eq!(documented_properties, properties, "properties in the API table versus the engine");
    }

    #[test]
    fn documentation_contains_every_signature() {
        let docs = include_str!("../../../docs/scripting.md");
        let missing: Vec<&str> = API.iter().map(|e| e.signature).filter(|s| !docs.contains(&format!("`{s}`"))).collect();
        assert!(missing.is_empty(), "docs/scripting.md lacks: {missing:?}");
    }

    #[test]
    fn text_has_every_section() {
        let text = api_text();
        for section in [FILES, OPENING, PROPERTIES, READING, EDITING, SAVING, ENZYMES, ANALYSIS, CLONING, PROTOCOLS, SEQUENCES, UTILITIES] {
            assert!(text.contains(section), "{section}");
        }
    }
}
