//! SnapGene's native files: `.dna` (DNA), `.rna` (RNA), `.prot` (protein).
//!
//! A SnapGene file is a list of packets. Each packet is one type byte, a 4-byte big-endian
//! length, then that many bytes of data. The first packet is the 0x09 "cookie" that identifies
//! the file. Annotations (features, primers, notes, settings) are short XML documents; the
//! sequence is plain ASCII after a flags byte. The layout was worked out from SnapGene documents
//! and public parsers such as Biopython's, never from SnapGene's program code.
//!
//! The reader fills in the format-neutral [`Document`] and keeps every packet in
//! [`SnapGeneRaw`]. The writer regenerates the packets the model covers and copies all others
//! through unchanged. Reading and writing an unedited file gives back the same bytes. After a
//! sequence edit it drops the packets derived from the old sequence and its history, and
//! empties the restriction-site cache so SnapGene rebuilds it.

use crate::error::{Error, Result, parse_err};
use crate::model::*;
use crate::span::Span;
use crate::xml::{self, Element, Node};

/// Packet type bytes, each with its meaning as found in SnapGene documents.
pub mod packet {
    pub const DNA: u8 = 0x00;
    /// Only found nested inside history nodes.
    pub const COMPRESSED_DNA: u8 = 0x01;
    /// 4-mer index of the sequence (derived).
    pub const KMER_INDEX: u8 = 0x02;
    /// Restriction-site cache (derived).
    pub const ENZYME_CACHE_3: u8 = 0x03;
    pub const PRIMERS: u8 = 0x05;
    pub const NOTES: u8 = 0x06;
    pub const HISTORY_TREE: u8 = 0x07;
    pub const PROPERTIES: u8 = 0x08;
    pub const COOKIE: u8 = 0x09;
    pub const FEATURES: u8 = 0x0A;
    pub const HISTORY_NODE: u8 = 0x0B;
    /// Enzyme display settings (the chosen enzyme set); not derived.
    pub const ENZYME_DISPLAY: u8 = 0x0D;
    pub const CUSTOM_ENZYME_SETS: u8 = 0x0E;
    /// Container for an attached sequence or trace.
    pub const SEQUENCE_CONTAINER: u8 = 0x10;
    pub const ALIGNABLE_SEQUENCES: u8 = 0x11;
    pub const PROTEIN: u8 = 0x15;
    /// BAM alignment of attached traces (derived).
    pub const TRACE_ALIGNMENT: u8 = 0x1B;
    pub const ENZYME_VISIBILITIES: u8 = 0x1C;
    pub const HISTORY_DELTA: u8 = 0x1D;
    pub const HISTORY_CONTENT: u8 = 0x1E;
    pub const RNA: u8 = 0x20;
    /// RNA secondary structure (derived).
    pub const SECONDARY_STRUCTURE: u8 = 0x22;
    /// Small xz-compressed JSON, probably saved searches.
    pub const SAVED_SEARCHES: u8 = 0x23;
}

/// Cookie sequence-type codes.
pub mod seq_type {
    pub const DNA: u16 = 1;
    pub const PROTEIN: u16 = 2;
    pub const DNA_ALIGNMENT: u16 = 3;
    pub const PROTEIN_ALIGNMENT: u16 = 4;
    pub const RNA: u16 = 7;
    pub const GEL: u16 = 8;
}

/// Versions written into the cookie of files SlopGene creates from scratch. They match
/// what current SnapGene releases write (export 15, import 19).
const DEFAULT_EXPORT_VERSION: u16 = 15;
const DEFAULT_IMPORT_VERSION: u16 = 19;

/// Packets derived from the sequence (or describing its history) that a writer drops after the
/// sequence was edited: the 4-mer index, history tree and nodes, the trace alignment, strand
/// colours, numbering and the secondary structure. The
/// restriction-site cache (0x03) is not dropped but emptied, which makes SnapGene rebuild it.
const STALE_AFTER_SEQUENCE_EDIT: &[u8] = &[0x02, 0x07, 0x0B, 0x14, 0x16, 0x1B, 0x1D, 0x1E, 0x22];

/// DNA packet flag bits.
const FLAG_CIRCULAR: u8 = 0x01;
const FLAG_DOUBLE_STRANDED: u8 = 0x02;
const FLAG_DAM: u8 = 0x04;
const FLAG_DCM: u8 = 0x08;
const FLAG_ECOKI: u8 = 0x10;

/// Whether the bytes start with a SnapGene cookie packet.
pub fn is_snapgene(bytes: &[u8]) -> bool {
    bytes.len() >= 13 && bytes[0] == packet::COOKIE && &bytes[5..13] == b"SnapGene"
}

/// Split a file into its packets.
pub fn split_packets(bytes: &[u8]) -> Result<Vec<RawPacket>> {
    let mut packets = Vec::new();
    let mut pos = 0usize;
    while pos < bytes.len() {
        if bytes.len() - pos < 5 {
            return Err(parse_err(format!("truncated packet header at byte {pos}")));
        }
        let kind = bytes[pos];
        let len = u32::from_be_bytes([bytes[pos + 1], bytes[pos + 2], bytes[pos + 3], bytes[pos + 4]]) as usize;
        let start = pos + 5;
        let end = start
            .checked_add(len)
            .filter(|&e| e <= bytes.len())
            .ok_or_else(|| parse_err(format!("packet 0x{kind:02x} at byte {pos} runs past the end of the file")))?;
        packets.push(RawPacket { kind, data: bytes[start..end].to_vec() });
        pos = end;
    }
    Ok(packets)
}

/// Read a SnapGene document.
pub fn read(bytes: &[u8]) -> Result<Document> {
    if !is_snapgene(bytes) {
        return Err(parse_err("not a SnapGene file (missing 'SnapGene' cookie)"));
    }
    let packets = split_packets(bytes)?;

    // Cookie: "SnapGene", then sequence type, export version and import version (u16 BE each).
    let cookie = &packets[0].data;
    if cookie.len() < 14 {
        return Err(parse_err("SnapGene cookie packet is too short"));
    }
    let be16 = |b: &[u8]| u16::from_be_bytes([b[0], b[1]]);
    let mut raw =
        SnapGeneRaw { seq_type: be16(&cookie[8..10]), export_version: be16(&cookie[10..12]), import_version: be16(&cookie[12..14]), ..SnapGeneRaw::default() };
    match raw.seq_type {
        seq_type::DNA | seq_type::PROTEIN | seq_type::RNA => {}
        seq_type::DNA_ALIGNMENT | seq_type::PROTEIN_ALIGNMENT => {
            return Err(Error::Unsupported("SnapGene alignment files (.dnaln/.praln) are not sequence documents".into()));
        }
        seq_type::GEL => return Err(Error::Unsupported("SnapGene gel files (.gel) are not sequence documents".into())),
        other => return Err(Error::Unsupported(format!("SnapGene document type {other}"))),
    }

    let mut doc = Document::default();
    let mut have_sequence = false;

    // First pass: the sequence, because feature coordinates need its length and topology.
    for pkt in &packets {
        let (molecule, is_seq) = match pkt.kind {
            packet::DNA => (Molecule::Dna, true),
            packet::RNA => (Molecule::Rna, true),
            packet::PROTEIN => (Molecule::Protein, true),
            _ => (Molecule::Dna, false),
        };
        if !is_seq {
            continue;
        }
        if have_sequence {
            return Err(parse_err("file contains more than one sequence packet"));
        }
        let Some((&flags, seq)) = pkt.data.split_first() else {
            return Err(parse_err("empty sequence packet"));
        };
        doc.molecule = molecule;
        doc.seq = seq.to_vec();
        doc.topology = if flags & FLAG_CIRCULAR != 0 { Topology::Circular } else { Topology::Linear };
        doc.strandedness = if molecule == Molecule::Dna && flags & FLAG_DOUBLE_STRANDED == 0 {
            Strandedness::Single
        } else if molecule == Molecule::Dna {
            Strandedness::Double
        } else {
            Strandedness::Single
        };
        doc.methylation = Methylation { dam: flags & FLAG_DAM != 0, dcm: flags & FLAG_DCM != 0, ecoki: flags & FLAG_ECOKI != 0 };
        have_sequence = true;
    }
    if !have_sequence {
        if packets.iter().any(|p| p.kind == packet::COMPRESSED_DNA) {
            return Err(Error::Unsupported("SnapGene compressed-DNA packet without a plain sequence".into()));
        }
        return Err(parse_err("SnapGene file has no sequence packet"));
    }

    // Second pass: annotations.
    let seq_len = doc.seq.len();
    let circular = doc.is_circular();
    for pkt in &packets {
        match pkt.kind {
            packet::FEATURES => {
                let root = parse_xml_packet(&pkt.data, "Features")?;
                let (features, root_rest) = read_features(&root, seq_len, circular)?;
                doc.features = features;
                raw.features_root = Some(root_rest);
            }
            packet::PRIMERS => {
                let root = parse_xml_packet(&pkt.data, "Primers")?;
                let (primers, root_rest) = read_primers(&root, seq_len, circular);
                doc.primers = primers;
                raw.primers_root = Some(root_rest);
            }
            packet::NOTES => {
                let root = parse_xml_packet(&pkt.data, "Notes")?;
                read_notes(&root, &mut doc.meta, &mut raw);
                raw.notes = Some(root);
            }
            packet::PROPERTIES => {
                let root = parse_xml_packet(&pkt.data, "AdditionalSequenceProperties")?;
                doc.ends = read_properties(&root);
                raw.properties = Some(root);
            }
            _ => {}
        }
    }

    raw.packets = packets;
    doc.raw.snapgene = Some(raw);
    doc.raw.snapgene.as_mut().unwrap().fingerprint = doc.sequence_fingerprint();
    Ok(doc)
}

fn parse_xml_packet(data: &[u8], what: &str) -> Result<Element> {
    let text = std::str::from_utf8(data).map_err(|_| parse_err(format!("{what} packet is not UTF-8")))?;
    xml::parse(text).map_err(|e| parse_err(format!("{what} packet: {e}")))
}

/// Parse a SnapGene range `a-b` (1-based inclusive) into a span. `a > b` runs through the origin.
fn parse_range_attr(text: &str, seq_len: usize, circular: bool) -> Option<Span> {
    let (a, b) = match text.split_once('-') {
        Some((a, b)) => (a.trim().parse::<usize>().ok()?, b.trim().parse::<usize>().ok()?),
        None => {
            let a = text.trim().parse::<usize>().ok()?;
            (a, a)
        }
    };
    if a == 0 || b == 0 || seq_len == 0 {
        return None;
    }
    let (a, b) = (a.min(seq_len), b.min(seq_len));
    if a <= b {
        Some(Span::new(a - 1, b))
    } else if circular {
        Some(Span::new(a - 1, b + seq_len))
    } else {
        Some(Span::new(b - 1, a))
    }
}

/// Format a span as a SnapGene range `a-b` (1-based inclusive).
fn format_range_attr(span: Span, seq_len: usize) -> String {
    let (first, last) = span.to_1based(seq_len);
    format!("{first}-{last}")
}

fn read_features(root: &Element, seq_len: usize, circular: bool) -> Result<(Vec<Feature>, Element)> {
    let mut features = Vec::new();
    let mut rest = Element { name: root.name.clone(), attrs: root.attrs.clone(), children: Vec::new() };
    for node in &root.children {
        match node {
            Node::Element(e) if e.name == "Feature" => {
                if let Some(feature) = read_feature(e, seq_len, circular) {
                    features.push(feature);
                }
            }
            Node::Element(e) => rest.push(e.clone()),
            _ => {}
        }
    }
    Ok((features, rest))
}

fn read_feature(e: &Element, seq_len: usize, circular: bool) -> Option<Feature> {
    let mut feature = Feature::new(String::new(), "misc_feature", Span::default(), Strand::None);
    feature.segments.clear();
    // Keep every attribute in order so the element can be written back as it was.
    feature.extra.attrs = e.attrs.clone();
    for (key, value) in &e.attrs {
        match key.as_str() {
            "name" => feature.name = value.clone(),
            "type" => feature.kind = value.clone(),
            "directionality" => {
                feature.strand = match value.trim() {
                    "1" => Strand::Forward,
                    "2" => Strand::Reverse,
                    "3" => Strand::Both,
                    _ => Strand::None,
                }
            }
            _ => {}
        }
    }
    for child in e.elements() {
        match child.name.as_str() {
            "Segment" => {
                let Some(span) = child.attr("range").and_then(|r| parse_range_attr(r, seq_len, circular)) else {
                    continue;
                };
                let mut segment = Segment::new(span);
                segment.extra = child.attrs.clone();
                for (key, value) in &child.attrs {
                    match key.as_str() {
                        // SnapGene also writes non-colours such as "noColor"; those stay in
                        // `extra` and are written back unchanged.
                        "color" => segment.color = Color::parse_hex(value),
                        "type" => segment.gap = value == "gap",
                        "name" => segment.name = Some(value.clone()),
                        "translated" => segment.translated = Some(value == "1"),
                        _ => {}
                    }
                }
                feature.segments.push(segment);
            }
            "Q" => {
                let Some(key) = child.attr("name") else { continue };
                let values = child.children_named("V").map(read_value).collect();
                feature.qualifiers.push(Qualifier { key: key.to_string(), values });
            }
            _ => feature.extra.children.push(child.clone()),
        }
    }
    if feature.segments.is_empty() {
        return None;
    }
    feature.normalize(seq_len, circular);
    Some(feature)
}

/// One `<V .../>` qualifier value. A value with a single text, int or predef attribute becomes
/// the matching typed value; anything else is kept as written.
fn read_value(v: &Element) -> QualValue {
    match v.attrs.as_slice() {
        [] => QualValue::Empty,
        [(k, value)] if k == "text" => QualValue::Text(value.clone()),
        [(k, value)] if k == "predef" => QualValue::Predef(value.clone()),
        [(k, value)] if k == "int" => match value.trim().parse::<i64>() {
            Ok(n) if n.to_string() == *value => QualValue::Int(n),
            _ => QualValue::Raw(v.attrs.clone()),
        },
        _ => QualValue::Raw(v.attrs.clone()),
    }
}

fn write_value(value: &QualValue) -> Element {
    let mut v = Element::new("V");
    match value {
        QualValue::Text(t) => v.set_attr("text", t.clone()),
        QualValue::Int(n) => v.set_attr("int", n.to_string()),
        QualValue::Predef(p) => v.set_attr("predef", p.clone()),
        QualValue::Empty => {}
        QualValue::Raw(attrs) => v.attrs = attrs.clone(),
    }
    v
}

/// Attributes for an element being written back.
///
/// `original` is the attribute list as read (empty for new elements). `modelled` lists the
/// attributes the model owns, with their current values (`None` = leave out). Original
/// attributes keep their order; modelled ones take their current value in place. Modelled
/// attributes the original lacked are slotted in at their position in SnapGene's `order`, but on
/// an element read from a file only when `add_missing` says so (so defaults SnapGene omitted
/// are not added).
fn merge_attrs(original: &[(String, String)], modelled: &[(&str, Option<String>, bool)], order: &[&str]) -> Vec<(String, String)> {
    let fresh = original.is_empty();
    let mut out = Vec::with_capacity(original.len() + modelled.len());
    for (key, value) in original {
        match modelled.iter().find(|(k, _, _)| k == key) {
            Some((_, Some(v), _)) => out.push((key.clone(), v.clone())),
            Some((_, None, _)) => {}
            None => out.push((key.clone(), value.clone())),
        }
    }
    let rank = |k: &str| order.iter().position(|o| *o == k);
    for (key, value, add_missing) in modelled {
        if original.iter().any(|(k, _)| k == key) {
            continue;
        }
        let Some(v) = value else { continue };
        if !(fresh || *add_missing) {
            continue;
        }
        // Insert before the first attribute SnapGene writes after this one.
        let index = match rank(key) {
            Some(r) => out.iter().position(|(k, _)| rank(k).is_some_and(|kr| kr > r)).unwrap_or(out.len()),
            None => out.len(),
        };
        out.insert(index, (key.to_string(), v.clone()));
    }
    out
}

fn read_primers(root: &Element, seq_len: usize, circular: bool) -> (Vec<Primer>, Element) {
    let mut primers = Vec::new();
    let mut rest = Element { name: root.name.clone(), attrs: root.attrs.clone(), children: Vec::new() };
    for node in &root.children {
        match node {
            Node::Element(e) if e.name == "Primer" => primers.push(read_primer(e, seq_len, circular)),
            Node::Element(e) => rest.push(e.clone()),
            _ => {}
        }
    }
    (primers, rest)
}

/// Primer binding-site locations are written `a-b` with `a` the 0-based first base and `b` the
/// 0-based last base (one less than feature ranges). Biopython's reader documents the same
/// offset, and the primers in SnapGene's own files match their sequences at it.
fn parse_primer_location(text: &str, seq_len: usize, circular: bool) -> Option<Span> {
    let (a, b) = text.split_once('-')?;
    let first = a.trim().parse::<usize>().ok()?;
    let last = b.trim().parse::<usize>().ok()?;
    if seq_len == 0 || first >= seq_len || last >= seq_len {
        return None;
    }
    if first <= last {
        Some(Span::new(first, last + 1))
    } else if circular {
        Some(Span::new(first, last + 1 + seq_len))
    } else {
        Some(Span::new(last, first + 1))
    }
}

fn format_primer_location(span: Span, seq_len: usize) -> String {
    let first = span.start;
    let last = (span.end - 1) % seq_len.max(1);
    format!("{first}-{last}")
}

fn read_primer(e: &Element, seq_len: usize, circular: bool) -> Primer {
    let mut primer = Primer::new(String::new(), String::new());
    primer.extra.attrs = e.attrs.clone();
    for (key, value) in &e.attrs {
        match key.as_str() {
            "name" => primer.name = value.clone(),
            "sequence" => primer.seq = value.clone(),
            "description" => primer.description = value.clone(),
            _ => {}
        }
    }
    for child in e.elements() {
        if child.name != "BindingSite" {
            primer.extra.children.push(child.clone());
            continue;
        }
        let Some(span) = child.attr("location").and_then(|l| parse_primer_location(l, seq_len, circular)) else {
            continue;
        };
        let mut site = BindingSite { span, strand: Strand::Forward, tm: None, extra: XmlExtra::default() };
        site.extra.attrs = child.attrs.clone();
        for (key, value) in &child.attrs {
            match key.as_str() {
                "boundStrand" => site.strand = if value.trim() == "1" { Strand::Reverse } else { Strand::Forward },
                "meltingTemperature" => site.tm = value.trim().parse().ok(),
                _ => {}
            }
        }
        site.extra.children = child.elements().cloned().collect();
        primer.sites.push(site);
    }
    primer
}

/// Convert a SnapGene date `2023.5.19` to ISO `2023-05-19`.
fn snapgene_date_to_iso(text: &str) -> Option<String> {
    let mut parts = text.trim().split('.');
    let year: u32 = parts.next()?.parse().ok()?;
    let month: u32 = parts.next()?.parse().ok()?;
    let day: u32 = parts.next()?.parse().ok()?;
    Some(format!("{year:04}-{month:02}-{day:02}"))
}

/// Convert an ISO date `2023-05-19` to SnapGene's `2023.5.19`.
fn iso_to_snapgene_date(text: &str) -> Option<String> {
    let mut parts = text.trim().split('-');
    let year: u32 = parts.next()?.parse().ok()?;
    let month: u32 = parts.next()?.parse().ok()?;
    let day: u32 = parts.next()?.parse().ok()?;
    Some(format!("{year}.{month}.{day}"))
}

/// Plain text from SnapGene's HTML snippets (`<html><body>...</body></html>`).
pub fn html_to_text(html: &str) -> String {
    let mut out = String::with_capacity(html.len());
    let mut rest = html;
    while let Some(lt) = rest.find('<') {
        out.push_str(&rest[..lt]);
        let Some(gt) = rest[lt..].find('>') else {
            out.push_str(&rest[lt..]);
            rest = "";
            break;
        };
        let tag = rest[lt + 1..lt + gt].trim().to_ascii_lowercase();
        // Line-breaking tags become newlines so paragraphs survive.
        if tag.starts_with("br") || tag == "/p" || tag == "/div" || tag == "/li" {
            out.push('\n');
        }
        rest = &rest[lt + gt + 1..];
    }
    out.push_str(rest);
    let text = xml::unescape(&out);
    text.lines().map(str::trim_end).collect::<Vec<_>>().join("\n").trim().to_string()
}

/// Wrap plain text as the HTML snippet SnapGene expects in descriptions.
pub fn text_to_html(text: &str) -> String {
    let body = xml::escape_text(text).replace('\n', "<br>");
    format!("<html><body>{body}</body></html>")
}

fn read_notes(root: &Element, meta: &mut Meta, raw: &mut SnapGeneRaw) {
    for child in root.elements() {
        let text = child.text();
        match child.name.as_str() {
            "Description" => {
                let plain = html_to_text(&text);
                raw.description_html = Some((text.clone(), plain.clone()));
                meta.description = plain;
            }
            "AccessionNumber" => meta.accession = Some(text.trim().to_string()).filter(|s| !s.is_empty()),
            "Organism" => meta.organism = Some(text.trim().to_string()).filter(|s| !s.is_empty()),
            "Type" => meta.synthetic = Some(text.trim().eq_ignore_ascii_case("synthetic")),
            "Created" => meta.created = snapgene_date_to_iso(&text),
            "LastModified" => meta.modified = snapgene_date_to_iso(&text),
            "CreatedBy" => meta.created_by = Some(text.trim().to_string()).filter(|s| !s.is_empty()),
            "Comments" => {
                let plain = html_to_text(&text);
                if !plain.is_empty() {
                    meta.comments.push(plain);
                }
            }
            "References" => {
                for r in child.children_named("Reference") {
                    meta.references.push(Reference {
                        authors: r.attr("authors").map(html_to_text).unwrap_or_default(),
                        title: r.attr("title").map(html_to_text).unwrap_or_default(),
                        journal: r.attr("journal").map(html_to_text).unwrap_or_default(),
                        pubmed: r.attr("pubMedID").map(str::to_string),
                    });
                }
            }
            _ => {}
        }
    }
}

/// Stickiness counts overhanging bases: positive for a 5' overhang, negative for a 3' one.
fn read_properties(root: &Element) -> Ends {
    let stick = |name: &str| root.child_text(name).and_then(|t| t.trim().parse::<i64>().ok()).unwrap_or(0);
    let to_overhang = |n: i64| match n {
        0 => Overhang::Blunt,
        n if n > 0 => Overhang::FivePrime(n as usize),
        n => Overhang::ThreePrime((-n) as usize),
    };
    Ends {
        left: to_overhang(stick("UpstreamStickiness")),
        right: to_overhang(stick("DownstreamStickiness")),
        left_modification: root.child_text("UpstreamModification"),
        right_modification: root.child_text("DownstreamModification"),
    }
}

// ---------------------------------------------------------------------------------------------
// Writing
//
// The rules below copy how SnapGene writes its own files.

/// Attribute orders SnapGene writes. New attributes are slotted into these positions.
const FEATURE_ATTR_ORDER: &[&str] = &[
    "recentID",
    "name",
    "directionality",
    "translateFirstCodonAsMet",
    "geneticCode",
    "translationMW",
    "type",
    "swappedSegmentNumbering",
    "subType",
    "allowSegmentOverlaps",
    "cleavageArrows",
    "readingFrame",
    "consecutiveTranslationNumbering",
    "consecutiveNumberingStartsFrom",
    "hitsStopCodon",
    "maxRunOn",
    "maxFusedRunOn",
    "visible",
    "prioritize",
    "originalName",
    "originalSequence",
    "isFavorite",
    "detectionMode",
];
const SEGMENT_ATTR_ORDER: &[&str] = &["name", "range", "color", "type", "translated", "translationNumberingStartsFrom"];
const PRIMER_ATTR_ORDER: &[&str] = &["recentID", "name", "sequence", "description", "visible", "minTailLength", "dateAdded"];
const SITE_ATTR_ORDER: &[&str] = &["simplified", "location", "boundStrand", "annealedBases", "meltingTemperature"];
/// Children of `<Notes>`, always in this order.
const NOTES_ORDER: &[&str] = &[
    "UUID",
    "Type",
    "ConfirmedExperimentally",
    "CustomMapLabel",
    "UseCustomMapLabel",
    "Product",
    "Description",
    "Created",
    "LastModified",
    "AccessionNumber",
    "CreatedBy",
    "Organism",
    "SequenceClass",
    "TransformedInto",
    "References",
    "Comments",
];
/// Children of `<AdditionalSequenceProperties>`, in the order SnapGene writes them.
const PROPERTIES_ORDER: &[&str] = &[
    "UpstreamEnzymeName",
    "DownstreamEnzymeName",
    "UpstreamStickiness",
    "DownstreamStickiness",
    "UpstreamModification",
    "DownstreamModification",
    "OriginalFirstBase",
];
/// Feature attributes SnapGene computes from the sequence. They are dropped after an edit.
const DERIVED_FEATURE_ATTRS: &[&str] = &["translationMW", "hitsStopCodon", "maxRunOn", "maxFusedRunOn"];

/// Write a document as a SnapGene file.
pub fn write(doc: &Document) -> Result<Vec<u8>> {
    let seq_type = match doc.molecule {
        Molecule::Dna => seq_type::DNA,
        Molecule::Rna => seq_type::RNA,
        Molecule::Protein => seq_type::PROTEIN,
    };
    let raw = doc.raw.snapgene.as_ref();
    let fresh = raw.is_none();
    let changed = doc.sequence_changed_since_load();
    let had = |kind: u8| raw.is_some_and(|r| r.packets.iter().any(|p| p.kind == kind));

    // After a sequence edit, data derived from the old sequence is stale: translations,
    // molecular weights, stop-codon positions and primer binding sites. Write from a copy with
    // those refreshed. New documents also get binding sites computed for their primers.
    let refreshed;
    let doc = if changed {
        refreshed = refresh_derived(doc, raw);
        &refreshed
    } else {
        doc
    };

    // When a regenerated XML packet means the same as the one read, write the original bytes,
    // so saving an unedited file changes nothing.
    let prefer_original = |new: RawPacket| -> RawPacket {
        let original = raw.and_then(|r| r.packets.iter().find(|p| p.kind == new.kind));
        match original {
            Some(orig) if same_xml(&orig.data, &new.data) => orig.clone(),
            _ => new,
        }
    };

    // The packets the model covers. Optional ones are only written when the original file had
    // them or there is something to put in them.
    let nucleic = doc.molecule != Molecule::Protein;
    let cookie = cookie_packet(seq_type, raw);
    let sequence = sequence_packet(doc);
    let properties =
        ((doc.molecule == Molecule::Dna && fresh) || had(packet::PROPERTIES)).then(|| prefer_original(properties_packet(doc, raw, changed && !fresh)));
    let features = (had(packet::FEATURES) || !doc.features.is_empty()).then(|| prefer_original(features_packet(doc, raw)));
    let primers = (nucleic && (fresh || had(packet::PRIMERS) || !doc.primers.is_empty())).then(|| prefer_original(primers_packet(doc, raw)));
    let notes = (fresh || had(packet::NOTES) || doc.meta != Meta::default()).then(|| prefer_original(notes_packet(doc, raw, changed)));
    // An empty restriction-site cache tells SnapGene to recompute it.
    let empty_cache = RawPacket { kind: packet::ENZYME_CACHE_3, data: Vec::new() };

    let mut out: Vec<RawPacket> = Vec::new();
    match raw {
        Some(raw) => {
            // Keep SnapGene's packet order: regenerated packets go where the originals stood,
            // everything else is copied, and stale packets are left out after an edit.
            let mut written = std::collections::HashSet::new();
            for pkt in &raw.packets {
                let slot = match pkt.kind {
                    packet::DNA | packet::RNA | packet::PROTEIN => packet::DNA,
                    other => other,
                };
                let regenerated = match pkt.kind {
                    packet::COOKIE => Some(Some(&cookie)),
                    packet::DNA | packet::RNA | packet::PROTEIN => Some(Some(&sequence)),
                    packet::PROPERTIES => Some(properties.as_ref()),
                    packet::FEATURES => Some(features.as_ref()),
                    packet::PRIMERS => Some(primers.as_ref()),
                    packet::NOTES => Some(notes.as_ref()),
                    _ => None,
                };
                if let Some(new_pkt) = regenerated {
                    if let Some(new_pkt) = new_pkt
                        && written.insert(slot)
                    {
                        out.push(new_pkt.clone());
                    }
                    continue;
                }
                if changed {
                    if STALE_AFTER_SEQUENCE_EDIT.contains(&pkt.kind) {
                        continue;
                    }
                    if pkt.kind == packet::ENZYME_CACHE_3 {
                        out.push(empty_cache.clone());
                        continue;
                    }
                    if pkt.kind == packet::CUSTOM_ENZYME_SETS {
                        out.push(without_highlighted_enzymes(pkt));
                        continue;
                    }
                }
                out.push(pkt.clone());
            }
            // Packets the original lacked (for example a first feature) go after the others,
            // except the cookie and sequence, which must lead.
            for (slot, pkt) in [
                (packet::COOKIE, Some(&cookie)),
                (packet::DNA, Some(&sequence)),
                (packet::PROPERTIES, properties.as_ref()),
                (packet::FEATURES, features.as_ref()),
                (packet::PRIMERS, primers.as_ref()),
                (packet::NOTES, notes.as_ref()),
            ] {
                if let Some(pkt) = pkt
                    && written.insert(slot)
                {
                    match slot {
                        packet::COOKIE => out.insert(0, pkt.clone()),
                        packet::DNA => out.insert(1.min(out.len()), pkt.clone()),
                        _ => out.push(pkt.clone()),
                    }
                }
            }
        }
        None => {
            // The smallest set SnapGene itself ships: cookie, sequence, empty site cache,
            // properties, features (when any), primers and notes.
            out.push(cookie);
            out.push(sequence);
            if doc.molecule == Molecule::Dna {
                out.push(empty_cache);
            }
            out.extend(properties);
            out.extend(features);
            out.extend(primers);
            out.extend(notes);
        }
    }

    let total: usize = out.iter().map(|p| p.data.len() + 5).sum();
    let mut bytes = Vec::with_capacity(total);
    for pkt in &out {
        let len = u32::try_from(pkt.data.len()).map_err(|_| Error::Invalid("packet larger than 4 GB".into()))?;
        bytes.push(pkt.kind);
        bytes.extend_from_slice(&len.to_be_bytes());
        bytes.extend_from_slice(&pkt.data);
    }
    Ok(bytes)
}

/// A copy of the document with sequence-derived data brought up to date.
fn refresh_derived(doc: &Document, raw: Option<&SnapGeneRaw>) -> Document {
    let mut doc = doc.clone();
    let snapshot = doc.clone();
    for feature in &mut doc.features {
        feature.extra.attrs.retain(|(k, _)| !DERIVED_FEATURE_ATTRS.contains(&k.as_str()));
        if feature.is_coding() && feature.qualifier("translation").is_some() {
            let protein = crate::translate::translate_feature(feature, &snapshot);
            feature.set_qualifier("translation", QualValue::Text(String::from_utf8_lossy(&protein).into_owned()));
        }
    }
    if !doc.primers.is_empty() {
        let rules = binding_rules(raw);
        crate::primers::refresh_binding_sites(&mut doc, &rules);
    }
    doc
}

/// Binding rules from the file's `HybridizationParams`, or SnapGene's defaults.
fn binding_rules(raw: Option<&SnapGeneRaw>) -> crate::primers::BindingRules {
    let mut rules = crate::primers::BindingRules::default();
    if let Some(params) = raw.and_then(|r| r.primers_root.as_ref()).and_then(|root| root.child("HybridizationParams")) {
        if let Some(n) = params.attr("minContinuousMatchLen").and_then(|v| v.parse().ok()) {
            rules.min_continuous_match = n;
        }
        if let Some(v) = params.attr("allowMismatch") {
            rules.allow_mismatch = v == "1";
        }
        if let Some(t) = params.attr("minMeltingTemperature").and_then(|v| v.parse().ok()) {
            rules.min_tm = t;
        }
    }
    rules
}

/// The custom-enzyme-sets packet without its `HighlightedEnzymeSet`, whose positions go stale
/// after an edit.
fn without_highlighted_enzymes(pkt: &RawPacket) -> RawPacket {
    let Ok(mut root) = parse_xml_packet(&pkt.data, "CustomEnzymeSets") else {
        return pkt.clone();
    };
    if root.child("HighlightedEnzymeSet").is_none() {
        return pkt.clone();
    }
    root.remove_children_named("HighlightedEnzymeSet");
    xml_packet(pkt.kind, &root, true)
}

/// Whether two XML packets mean the same thing (same canonical form).
fn same_xml(a: &[u8], b: &[u8]) -> bool {
    if a == b {
        return true;
    }
    let parse = |d: &[u8]| std::str::from_utf8(d).ok().and_then(|t| xml::parse(t).ok());
    match (parse(a), parse(b)) {
        (Some(x), Some(y)) => xml::canonical(&x) == xml::canonical(&y),
        _ => false,
    }
}

fn cookie_packet(seq_type: u16, raw: Option<&SnapGeneRaw>) -> RawPacket {
    let (export, import) = match raw {
        Some(r) if r.seq_type == seq_type => (r.export_version, r.import_version),
        _ => (DEFAULT_EXPORT_VERSION, DEFAULT_IMPORT_VERSION),
    };
    let mut data = b"SnapGene".to_vec();
    data.extend_from_slice(&seq_type.to_be_bytes());
    data.extend_from_slice(&export.to_be_bytes());
    data.extend_from_slice(&import.to_be_bytes());
    RawPacket { kind: packet::COOKIE, data }
}

fn sequence_packet(doc: &Document) -> RawPacket {
    let mut flags = 0u8;
    if doc.is_circular() {
        flags |= FLAG_CIRCULAR;
    }
    if doc.molecule == Molecule::Dna {
        if doc.strandedness == Strandedness::Double {
            flags |= FLAG_DOUBLE_STRANDED;
        }
        if doc.methylation.dam {
            flags |= FLAG_DAM;
        }
        if doc.methylation.dcm {
            flags |= FLAG_DCM;
        }
        if doc.methylation.ecoki {
            flags |= FLAG_ECOKI;
        }
    }
    let kind = match doc.molecule {
        Molecule::Dna => packet::DNA,
        Molecule::Rna => packet::RNA,
        Molecule::Protein => packet::PROTEIN,
    };
    let mut data = Vec::with_capacity(doc.seq.len() + 1);
    data.push(flags);
    data.extend_from_slice(&doc.seq);
    RawPacket { kind, data }
}

/// Serialise an XML packet. SnapGene puts a declaration and a trailing newline on features,
/// primers and enzyme packets, and neither on notes or properties.
fn xml_packet(kind: u8, root: &Element, declaration: bool) -> RawPacket {
    let mut text = String::new();
    if declaration {
        text.push_str("<?xml version=\"1.0\"?>");
    }
    root.write_to(&mut text);
    if declaration {
        text.push('\n');
    }
    RawPacket { kind, data: text.into_bytes() }
}

/// Set the text of child `name`, creating it at its canonical position (by `order`) when
/// missing. `newline` keeps a one-child-per-line layout by adding a line break after it.
fn set_child_ordered(root: &mut Element, name: &str, text: String, order: &[&str], newline: bool) {
    if let Some(child) = root.child_mut(name) {
        child.set_text(text);
        return;
    }
    let rank = |n: &str| order.iter().position(|o| *o == n);
    let my_rank = rank(name).unwrap_or(usize::MAX);
    let index =
        root.children.iter().position(|node| matches!(node, Node::Element(e) if rank(&e.name).is_some_and(|r| r > my_rank))).unwrap_or(root.children.len());
    let mut child = Element::new(name);
    child.set_text(text);
    if newline {
        root.children.insert(index, Node::Text("\n".into()));
    }
    root.children.insert(index, Node::Element(child));
}

/// Remove child `name` and the line break that follows it.
fn remove_child(root: &mut Element, name: &str) {
    while let Some(index) = root.children.iter().position(|n| matches!(n, Node::Element(e) if e.name == name)) {
        root.children.remove(index);
        if matches!(root.children.get(index), Some(Node::Text(t)) if t.trim().is_empty()) {
            root.children.remove(index);
        }
    }
}

fn properties_packet(doc: &Document, raw: Option<&SnapGeneRaw>, sequence_edited: bool) -> RawPacket {
    let mut root = raw.and_then(|r| r.properties.clone()).unwrap_or_else(|| Element::new("AdditionalSequenceProperties"));
    let stickiness = |o: Overhang| match o {
        Overhang::Blunt => 0i64,
        Overhang::FivePrime(n) => n as i64,
        Overhang::ThreePrime(n) => -(n as i64),
    };
    // Circular molecules have no ends, so SnapGene writes zero stickiness for them.
    let (left, right) = if doc.is_circular() { (0, 0) } else { (stickiness(doc.ends.left), stickiness(doc.ends.right)) };
    let modification = |m: &Option<String>| m.clone().unwrap_or_else(|| "Unmodified".to_string());
    set_child_ordered(&mut root, "UpstreamStickiness", left.to_string(), PROPERTIES_ORDER, false);
    set_child_ordered(&mut root, "DownstreamStickiness", right.to_string(), PROPERTIES_ORDER, false);
    set_child_ordered(&mut root, "UpstreamModification", modification(&doc.ends.left_modification), PROPERTIES_ORDER, false);
    set_child_ordered(&mut root, "DownstreamModification", modification(&doc.ends.right_modification), PROPERTIES_ORDER, false);
    // The original numbering anchor is a position, so it is meaningless after an edit.
    if sequence_edited {
        remove_child(&mut root, "OriginalFirstBase");
    }
    xml_packet(packet::PROPERTIES, &root, false)
}

/// The largest `recentID` in use plus one, so new elements get fresh ids.
fn next_recent_id<'a>(extras: impl Iterator<Item = &'a XmlExtra>, root: Option<&Element>) -> u64 {
    let mut next = root.and_then(|r| r.attr("nextValidID")).and_then(|v| v.parse::<u64>().ok()).unwrap_or(0);
    for extra in extras {
        if let Some(id) = extra.attr("recentID").and_then(|v| v.parse::<u64>().ok()) {
            next = next.max(id + 1);
        }
    }
    next
}

fn features_packet(doc: &Document, raw: Option<&SnapGeneRaw>) -> RawPacket {
    let seq_len = doc.seq.len();
    let original_root = raw.and_then(|r| r.features_root.as_ref());
    let mut next_id = next_recent_id(doc.features.iter().map(|f| &f.extra), original_root);
    let mut root = original_root.cloned().unwrap_or_else(|| Element::new("Features"));

    let mut feature_elements = Vec::with_capacity(doc.features.len());
    for feature in &doc.features {
        let mut e = Element::new("Feature");
        // Features we created get a fresh recentID; loaded ones keep theirs.
        let recent_id = match feature.extra.attr("recentID") {
            Some(id) => id.to_string(),
            None => {
                next_id += 1;
                (next_id - 1).to_string()
            }
        };
        let directionality = match feature.strand {
            Strand::None => None,
            Strand::Forward => Some("1".to_string()),
            Strand::Reverse => Some("2".to_string()),
            Strand::Both => Some("3".to_string()),
        };
        e.attrs = merge_attrs(
            &feature.extra.attrs,
            &[
                ("recentID", Some(recent_id), true),
                ("name", Some(feature.name.clone()), true),
                ("directionality", directionality, true),
                ("type", Some(feature.kind.clone()), true),
            ],
            FEATURE_ATTR_ORDER,
        );
        let default_color = crate::render::theme::feature_color(&feature.kind).hex();
        for segment in &feature.segments {
            // SnapGene always writes a colour: "noColor" for gaps, and the feature-type default
            // for segments that never had one. A placeholder read from the file is kept unless a
            // colour has been set since.
            let original_color = segment.extra.iter().find(|(k, _)| k == "color").map(|(_, v)| v.clone());
            let color = match (segment.color, original_color) {
                (Some(c), _) => Some(c.to_hex()),
                (None, Some(orig)) if Color::parse_hex(&orig).is_none() => Some(orig),
                (None, Some(_)) => None,
                (None, None) if segment.gap => Some("noColor".to_string()),
                (None, None) => Some(default_color.clone()),
            };
            let mut s = Element::new("Segment");
            s.attrs = merge_attrs(
                &segment.extra,
                &[
                    ("name", segment.name.clone(), true),
                    ("range", Some(format_range_attr(segment.span, seq_len)), true),
                    ("color", color, true),
                    ("type", Some(if segment.gap { "gap" } else { "standard" }.to_string()), true),
                    ("translated", segment.translated.map(|t| if t { "1" } else { "0" }.to_string()), true),
                ],
                SEGMENT_ATTR_ORDER,
            );
            e.push(s);
        }
        for qualifier in &feature.qualifiers {
            let mut q = Element::new("Q").with_attr("name", qualifier.key.clone());
            for value in &qualifier.values {
                q.push(write_value(value));
            }
            e.push(q);
        }
        for child in &feature.extra.children {
            e.push(child.clone());
        }
        feature_elements.push(e);
    }

    // Features come first under the root, then any other children the original root had.
    let others = std::mem::take(&mut root.children);
    root.children = feature_elements.into_iter().map(Node::Element).collect();
    root.children.extend(others);
    root.set_attr("nextValidID", next_id.to_string());
    xml_packet(packet::FEATURES, &root, true)
}

fn default_hybridization_params() -> Element {
    // SnapGene's defaults: a primer binds where at least 10 bases at its 3' end match
    // contiguously, mismatches are allowed elsewhere, and the annealed part must reach 40 C.
    Element::new("HybridizationParams")
        .with_attr("minContinuousMatchLen", "10")
        .with_attr("allowMismatch", "1")
        .with_attr("minMeltingTemperature", "40")
        .with_attr("showAdditionalFivePrimeMatches", "1")
        .with_attr("minimumFivePrimeAnnealing", "15")
}

/// A binding site computed by SlopGene, written the way SnapGene writes its own: the primer
/// split into components (an unpaired 5' tail, then the paired 3' part with its template
/// range), the paired bases, the Tm, and a `simplified` display twin.
fn fresh_site_elements(primer: &Primer, site: &BindingSite, seq_len: usize) -> [Element; 2] {
    let n = primer.seq.len();
    let annealed = site.span.len().min(n);
    let (tail, paired) = primer.seq.split_at(n - annealed);
    let location = format_primer_location(site.span, seq_len);
    let tm = site.tm.unwrap_or_else(|| crate::primers::melting_temperature(paired.as_bytes(), &crate::primers::TmConditions::default()));
    let mut full = Element::new("BindingSite")
        .with_attr("location", location.clone())
        .with_attr("boundStrand", if site.strand.is_reverse() { "1" } else { "0" })
        .with_attr("annealedBases", paired)
        .with_attr("meltingTemperature", format!("{}", tm.round() as i64));
    if !tail.is_empty() {
        full.push(Element::new("Component").with_attr("bases", tail));
    }
    full.push(Element::new("Component").with_attr("hybridizedRange", location).with_attr("bases", paired));
    let mut twin = full.clone();
    twin.attrs.insert(0, ("simplified".to_string(), "1".to_string()));
    [full, twin]
}

fn primers_packet(doc: &Document, raw: Option<&SnapGeneRaw>) -> RawPacket {
    let seq_len = doc.seq.len();
    let original_root = raw.and_then(|r| r.primers_root.as_ref());
    let mut next_id = next_recent_id(doc.primers.iter().map(|p| &p.extra), original_root);
    let mut root = original_root.cloned().unwrap_or_else(|| {
        let mut r = Element::new("Primers");
        r.push(default_hybridization_params());
        r
    });

    let mut primer_elements = Vec::with_capacity(doc.primers.len());
    for primer in &doc.primers {
        let mut e = Element::new("Primer");
        let recent_id = match primer.extra.attr("recentID") {
            Some(id) => id.to_string(),
            None => {
                next_id += 1;
                (next_id - 1).to_string()
            }
        };
        e.attrs = merge_attrs(
            &primer.extra.attrs,
            &[
                ("recentID", Some(recent_id), true),
                ("name", Some(primer.name.clone()), true),
                ("sequence", Some(primer.seq.clone()), true),
                ("description", Some(primer.description.clone()).filter(|d| !d.is_empty()), true),
            ],
            PRIMER_ATTR_ORDER,
        );
        for site in &primer.sites {
            if seq_len == 0 {
                break;
            }
            // Sites SlopGene computed carry no SnapGene attributes yet: write them in full.
            if site.extra.attrs.is_empty() {
                let [full, twin] = fresh_site_elements(primer, site, seq_len);
                e.push(full);
                e.push(twin);
                continue;
            }
            // Keep the stored temperature text when the value has not changed.
            let original_tm = site.extra.attr("meltingTemperature");
            let tm = site.tm.map(|t| match original_tm {
                Some(orig) if orig.trim().parse::<f64>().is_ok_and(|o| (o - t).abs() < 0.05) => orig.to_string(),
                _ => format!("{}", t.round() as i64),
            });
            let mut s = Element::new("BindingSite");
            s.attrs = merge_attrs(
                &site.extra.attrs,
                &[
                    ("location", Some(format_primer_location(site.span, seq_len)), true),
                    ("boundStrand", Some(if site.strand.is_reverse() { "1" } else { "0" }.to_string()), true),
                    ("meltingTemperature", tm, true),
                ],
                SITE_ATTR_ORDER,
            );
            for child in &site.extra.children {
                s.push(child.clone());
            }
            e.push(s);
        }
        for child in &primer.extra.children {
            e.push(child.clone());
        }
        primer_elements.push(e);
    }

    // HybridizationParams and other non-primer children stay ahead of the primers, as SnapGene
    // writes them.
    root.children.extend(primer_elements.into_iter().map(Node::Element));
    root.set_attr("nextValidID", next_id.to_string());
    xml_packet(packet::PRIMERS, &root, true)
}

/// `H:M:S` of the current UTC time, unpadded as SnapGene writes it.
fn utc_time_attr() -> String {
    let (_, _, _, h, m, s) = crate::clock::now_utc();
    format!("{h}:{m}:{s}")
}

fn notes_packet(doc: &Document, raw: Option<&SnapGeneRaw>, sequence_edited: bool) -> RawPacket {
    let meta = &doc.meta;
    let original = raw.and_then(|r| r.notes.as_ref());
    let today = crate::clock::today_iso();
    let mut root = original.cloned().unwrap_or_else(|| {
        // A new document: SnapGene's usual set of notes, one per line.
        let mut r = Element::new("Notes");
        r.children.push(Node::Text("\n".into()));
        let now = utc_time_attr();
        let created = meta.created.clone().unwrap_or_else(|| today.clone());
        let modified = meta.modified.clone().unwrap_or_else(|| today.clone());
        let mut add = |name: &str, text: String, utc: Option<&str>| {
            let mut child = Element::new(name);
            if let Some(t) = utc {
                child.set_attr("UTC", t);
            }
            child.set_text(text);
            r.children.push(Node::Element(child));
            r.children.push(Node::Text("\n".into()));
        };
        add("UUID", pseudo_uuid(doc), None);
        add("Type", if meta.synthetic == Some(false) { "Natural" } else { "Synthetic" }.into(), None);
        add("ConfirmedExperimentally", "0".into(), None);
        add("Created", iso_to_snapgene_date(&created).unwrap_or_default(), Some(&now));
        add("LastModified", iso_to_snapgene_date(&modified).unwrap_or_default(), Some(&now));
        add("SequenceClass", "UNA".into(), None);
        add("TransformedInto", "Unspecified".into(), None);
        r
    });

    // Description: keep SnapGene's original HTML when the text did not change.
    match raw.and_then(|r| r.description_html.as_ref()) {
        Some((html, plain)) if *plain == meta.description => set_child_ordered(&mut root, "Description", html.clone(), NOTES_ORDER, true),
        _ if meta.description.is_empty() => remove_child(&mut root, "Description"),
        _ => set_child_ordered(&mut root, "Description", text_to_html(&meta.description), NOTES_ORDER, true),
    }
    match &meta.accession {
        Some(acc) => set_child_ordered(&mut root, "AccessionNumber", acc.clone(), NOTES_ORDER, true),
        None => remove_child(&mut root, "AccessionNumber"),
    }
    if let Some(synthetic) = meta.synthetic {
        set_child_ordered(&mut root, "Type", if synthetic { "Synthetic" } else { "Natural" }.into(), NOTES_ORDER, true);
    }
    if let Some(by) = &meta.created_by {
        set_child_ordered(&mut root, "CreatedBy", by.clone(), NOTES_ORDER, true);
    }
    if let Some(organism) = &meta.organism
        && (original.is_none() || root.child("Organism").is_some())
    {
        set_child_ordered(&mut root, "Organism", organism.clone(), NOTES_ORDER, true);
    }

    // Last modified: stamped with the current time after a sequence edit, or when the caller
    // changed the date.
    let original_modified = original.and_then(|n| n.child_text("LastModified")).and_then(|t| snapgene_date_to_iso(&t));
    let stamp = if sequence_edited && original.is_some() {
        Some(today.clone())
    } else if meta.modified.is_some() && meta.modified != original_modified && original.is_some() {
        meta.modified.clone()
    } else {
        None
    };
    if let Some(date) = stamp.as_deref().and_then(iso_to_snapgene_date) {
        set_child_ordered(&mut root, "LastModified", date, NOTES_ORDER, true);
        if let Some(child) = root.child_mut("LastModified") {
            child.set_attr("UTC", utc_time_attr());
        }
    }
    let mut packet = xml_packet(packet::NOTES, &root, false);
    if original.is_none() {
        packet.data.push(b'\n');
    }
    packet
}

/// A UUID-shaped identifier for new documents, derived from the sequence, the time and a
/// process counter. SnapGene only needs it to be unique.
fn pseudo_uuid(doc: &Document) -> String {
    let nanos = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_nanos() as u64).unwrap_or(0);
    let a = doc.sequence_fingerprint() ^ nanos.rotate_left(17);
    let b = nanos ^ new_uid().wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ ((std::process::id() as u64) << 32);
    format!(
        "{:08x}-{:04x}-4{:03x}-{:04x}-{:012x}",
        (a >> 32) as u32,
        (a >> 16) as u16,
        (a as u16) & 0x0fff,
        ((b >> 48) as u16 & 0x3fff) | 0x8000,
        b & 0x0000_ffff_ffff_ffff
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_doc() -> Document {
        let mut doc = Document::new_dna("test", b"ATGAAACCCGGGTTTTAAGAATTC", Topology::Circular);
        let mut f = Feature::new("orf", "CDS", Span::new(0, 18), Strand::Forward);
        f.set_color(Some(Color::rgb(0x99, 0xcc, 0xff)));
        f.set_qualifier("note", QualValue::Text("a <b>note</b>".into()));
        doc.features.push(f);
        // A feature that runs through the origin.
        doc.features.push(Feature::new("wrap", "misc_feature", Span::new(20, 28), Strand::Reverse));
        let mut p = Primer::new("p1", "ATGAAACCCGGG");
        p.sites.push(BindingSite { span: Span::new(0, 12), strand: Strand::Forward, tm: Some(42.0), extra: XmlExtra::default() });
        doc.primers.push(p);
        doc.meta.description = "Line one\nLine & two".into();
        doc
    }

    #[test]
    fn write_then_read_round_trip() {
        let doc = sample_doc();
        let bytes = write(&doc).unwrap();
        assert!(is_snapgene(&bytes));
        let back = read(&bytes).unwrap();
        assert_eq!(back.seq, doc.seq);
        assert_eq!(back.topology, Topology::Circular);
        assert_eq!(back.features.len(), 2);
        assert_eq!(back.features[0].name, "orf");
        assert_eq!(back.features[0].segments[0].span, Span::new(0, 18));
        assert_eq!(back.features[0].color(), Some(Color::rgb(0x99, 0xcc, 0xff)));
        assert_eq!(back.features[0].qualifier_text("note").as_deref(), Some("a <b>note</b>"));
        assert_eq!(back.features[1].span(doc.len()), Span::new(20, 28));
        assert_eq!(back.features[1].strand, Strand::Reverse);
        assert_eq!(back.primers[0].sites[0].span, Span::new(0, 12));
        assert_eq!(back.meta.description, "Line one\nLine & two");
    }

    #[test]
    fn rewrite_unchanged_keeps_unknown_packets() {
        let doc = sample_doc();
        let mut bytes = write(&doc).unwrap();
        // Append an unknown packet and a history packet, as if SnapGene had written them.
        for (kind, payload) in [(0x42u8, b"opaque".as_slice()), (packet::HISTORY_TREE, b"history".as_slice())] {
            bytes.push(kind);
            bytes.extend_from_slice(&(payload.len() as u32).to_be_bytes());
            bytes.extend_from_slice(payload);
        }
        let loaded = read(&bytes).unwrap();
        let again = write(&loaded).unwrap();
        let kinds: Vec<u8> = split_packets(&again).unwrap().iter().map(|p| p.kind).collect();
        assert!(kinds.contains(&0x42));
        assert!(kinds.contains(&packet::HISTORY_TREE));

        // After a sequence edit the history goes, the unknown packet stays.
        let mut edited = loaded.clone();
        edited.seq.push(b'A');
        let kinds: Vec<u8> = split_packets(&write(&edited).unwrap()).unwrap().iter().map(|p| p.kind).collect();
        assert!(kinds.contains(&0x42));
        assert!(!kinds.contains(&packet::HISTORY_TREE));
    }

    #[test]
    fn html_conversion() {
        assert_eq!(html_to_text("<html><body>Standard <i>E. coli</i> vector</body></html>"), "Standard E. coli vector");
        assert_eq!(html_to_text("a<br>b &amp; c"), "a\nb & c");
        assert_eq!(text_to_html("a\nb<c"), "<html><body>a<br>b&lt;c</body></html>");
    }
}
