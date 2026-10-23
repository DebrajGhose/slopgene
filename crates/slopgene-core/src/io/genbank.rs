//! GenBank flat files (`.gb`, `.gbk`, `.genbank`), read and written.
//!
//! The reader is lenient: it accepts the long locus names, missing header fields and odd
//! spacing that plasmid tools (SnapGene, Benchling, ApE, pLannotate, Addgene) produce. The writer
//! follows the NCBI column layout so the files load everywhere.

use crate::error::{Result, parse_err};
use crate::model::*;
use crate::span::Span;

/// Qualifiers whose values are written without quotes.
const NUMERIC_QUALIFIERS: &[&str] = &["codon_start", "transl_table", "number", "estimated_length"];

/// Read every record in a GenBank file.
pub fn read(text: &str) -> Result<Vec<Document>> {
    let mut docs = Vec::new();
    let mut lines: Vec<&str> = Vec::new();
    for line in text.lines() {
        if line.starts_with("//") {
            if lines.iter().any(|l| !l.trim().is_empty()) {
                docs.push(read_record(&lines)?);
            }
            lines.clear();
        } else {
            lines.push(line);
        }
    }
    // Tolerate a missing final "//".
    if lines.iter().any(|l| l.starts_with("LOCUS")) {
        docs.push(read_record(&lines)?);
    }
    if docs.is_empty() {
        return Err(parse_err("no GenBank records found (expected a LOCUS line)"));
    }
    Ok(docs)
}

/// One raw feature as read from the feature table, before coordinates are resolved.
struct RawFeature {
    key: String,
    location: String,
    qualifiers: Vec<(String, Option<String>)>,
}

fn read_record(lines: &[&str]) -> Result<Document> {
    let mut doc = Document::default();
    let mut declared_len: Option<usize> = None;
    let mut raw_features: Vec<RawFeature> = Vec::new();
    let mut seq = Vec::new();

    // `section` is the header keyword we are inside, so continuation lines know where to go.
    let mut section = String::new();
    let mut in_features = false;
    let mut in_origin = false;
    let mut definition = String::new();
    let mut comment = String::new();
    let mut reference: Option<Reference> = None;

    // State for the qualifier currently being read inside the feature table.
    let mut open_qualifier: Option<(String, String)> = None;

    let flush_qualifier = |raw_features: &mut Vec<RawFeature>, open: &mut Option<(String, String)>| {
        if let Some((key, value)) = open.take()
            && let Some(f) = raw_features.last_mut()
        {
            f.qualifiers.push((key, Some(value)));
        }
    };

    for &line in lines {
        if in_origin {
            // Sequence lines are a position number followed by blocks of ten bases.
            seq.extend(line.bytes().filter(|b| b.is_ascii_alphabetic()));
            continue;
        }

        if in_features {
            // The table ends at the next keyword in column 1.
            if !line.starts_with(' ') && !line.trim().is_empty() {
                flush_qualifier(&mut raw_features, &mut open_qualifier);
                in_features = false;
            } else {
                read_feature_line(line, &mut raw_features, &mut open_qualifier);
                continue;
            }
        }

        // Header lines: a keyword in columns 1-12, or a continuation of the previous one.
        let keyword = line.get(..12).unwrap_or(line).trim();
        let value = line.get(12..).unwrap_or("").trim_end();
        let is_continuation = line.starts_with("            ") || (keyword.is_empty() && !line.trim().is_empty());
        if is_continuation {
            match section.as_str() {
                "DEFINITION" => {
                    definition.push(' ');
                    definition.push_str(value.trim());
                }
                "COMMENT" => {
                    comment.push('\n');
                    comment.push_str(value.trim());
                }
                "AUTHORS" | "TITLE" | "JOURNAL" => {
                    if let Some(r) = reference.as_mut() {
                        let target = match section.as_str() {
                            "AUTHORS" => &mut r.authors,
                            "TITLE" => &mut r.title,
                            _ => &mut r.journal,
                        };
                        target.push(' ');
                        target.push_str(value.trim());
                    }
                }
                "KEYWORDS" => {
                    if let Some(k) = doc.meta.keywords.as_mut() {
                        k.push(' ');
                        k.push_str(value.trim());
                    }
                }
                _ => {}
            }
            continue;
        }
        match keyword {
            "LOCUS" => {
                section = "LOCUS".into();
                parse_locus(line, &mut doc, &mut declared_len);
            }
            "DEFINITION" => {
                section = "DEFINITION".into();
                definition = value.trim().to_string();
            }
            "ACCESSION" => {
                section = "ACCESSION".into();
                let acc = value.split_whitespace().next().unwrap_or("").to_string();
                doc.meta.accession = Some(acc).filter(|a| !a.is_empty() && a != ".");
            }
            "VERSION" => {
                section = "VERSION".into();
                let v = value.split_whitespace().next().unwrap_or("").to_string();
                doc.meta.version = Some(v).filter(|v| !v.is_empty() && v != ".");
            }
            "KEYWORDS" => {
                section = "KEYWORDS".into();
                let k = value.trim().trim_end_matches('.').to_string();
                doc.meta.keywords = Some(k).filter(|k| !k.is_empty());
            }
            "SOURCE" => section = "SOURCE".into(),
            "ORGANISM" => {
                section = "ORGANISM".into();
                let o = value.trim().to_string();
                doc.meta.organism = Some(o).filter(|o| !o.is_empty() && o != ".");
            }
            "REFERENCE" => {
                section = "REFERENCE".into();
                if let Some(r) = reference.take() {
                    doc.meta.references.push(r);
                }
                reference = Some(Reference::default());
            }
            "AUTHORS" | "TITLE" | "JOURNAL" => {
                section = keyword.to_string();
                if let Some(r) = reference.as_mut() {
                    let target = match keyword {
                        "AUTHORS" => &mut r.authors,
                        "TITLE" => &mut r.title,
                        _ => &mut r.journal,
                    };
                    *target = value.trim().to_string();
                }
            }
            "PUBMED" => {
                if let Some(r) = reference.as_mut() {
                    r.pubmed = Some(value.trim().to_string());
                }
            }
            "COMMENT" => {
                section = "COMMENT".into();
                if !comment.is_empty() {
                    comment.push('\n');
                }
                comment.push_str(value.trim());
            }
            "FEATURES" => {
                section = "FEATURES".into();
                in_features = true;
            }
            "ORIGIN" => {
                in_origin = true;
            }
            _ => section = keyword.to_string(),
        }
    }
    flush_qualifier(&mut raw_features, &mut open_qualifier);
    if let Some(r) = reference.take() {
        doc.meta.references.push(r);
    }

    doc.seq = seq;
    if let Some(declared) = declared_len
        && declared != doc.seq.len()
        && !doc.seq.is_empty()
    {
        return Err(parse_err(format!("LOCUS says {declared} bp but the sequence has {} bases", doc.seq.len())));
    }
    let definition = definition.trim().trim_end_matches('.').trim().to_string();
    doc.meta.description = definition;
    if !comment.trim().is_empty() {
        doc.meta.comments.push(comment.trim().to_string());
    }

    // Resolve features now that the length and topology are known.
    let seq_len = doc.seq.len();
    let circular = doc.is_circular();
    for raw in raw_features {
        let qualifiers = to_qualifiers(raw.qualifiers);
        if raw.key == "source" {
            // The source feature describes the whole record; fold it into the metadata.
            if let Some(org) = qualifiers.iter().find(|q| q.key == "organism").and_then(|q| q.values.first()) {
                let org = org.as_text();
                if doc.meta.organism.is_none() && !org.is_empty() && org != "." {
                    doc.meta.organism = Some(org);
                }
            }
            doc.raw.genbank_source = Some(qualifiers);
            continue;
        }
        let Ok((segments, strand)) = parse_location(&raw.location, seq_len, circular) else {
            continue;
        };
        if segments.is_empty() {
            continue;
        }
        let mut feature = Feature::new(String::new(), raw.key.clone(), Span::default(), strand);
        feature.segments = segments.into_iter().map(Segment::new).collect();
        feature.qualifiers = qualifiers;
        feature.name = feature_label(&feature);
        apply_color_qualifiers(&mut feature);
        feature.normalize(seq_len, circular);
        doc.features.push(feature);
    }
    Ok(doc)
}

/// Handle one line of the feature table.
///
/// A line is one of: the start of a feature (key indented by about 5 columns, then the
/// location), the start of a qualifier (`/key=value`, indented by 21), or a continuation of the
/// location, of an unquoted value or of a quoted value that has not closed yet.
fn read_feature_line(line: &str, features: &mut Vec<RawFeature>, open: &mut Option<(String, String)>) {
    let content = line.trim();
    if content.is_empty() {
        return;
    }

    // Inside a quoted value everything is text until the closing quote. Translations are
    // joined without spaces; other text is joined with one space.
    if let Some((key, value)) = open.as_mut() {
        if key != "translation" {
            value.push(' ');
        }
        value.push_str(content);
        if quote_closed(value) {
            let finished = finish_quoted(value);
            let key = std::mem::take(key);
            if let Some(f) = features.last_mut() {
                f.qualifiers.push((key, Some(finished)));
            }
            *open = None;
        }
        return;
    }

    // A feature key sits near the left margin; qualifiers and continuations are indented ~21.
    let indent = line.len() - line.trim_start().len();
    if indent <= 10 && !content.starts_with('/') {
        let (key, location) = content.split_once(char::is_whitespace).unwrap_or((content, ""));
        features.push(RawFeature { key: key.to_string(), location: location.trim().to_string(), qualifiers: Vec::new() });
        return;
    }

    if let Some(q) = content.strip_prefix('/') {
        let Some(f) = features.last_mut() else { return };
        match q.split_once('=') {
            None => f.qualifiers.push((q.trim().to_string(), None)),
            Some((key, value)) => {
                let key = key.trim().to_string();
                let value = value.trim();
                if value.starts_with('"') {
                    if quote_closed(value) {
                        f.qualifiers.push((key, Some(finish_quoted(value))));
                    } else {
                        *open = Some((key, value.to_string()));
                    }
                } else {
                    f.qualifiers.push((key, Some(value.to_string())));
                }
            }
        }
        return;
    }

    // Continuation of the location (before any qualifier) or of an unquoted value.
    if let Some(f) = features.last_mut() {
        if f.qualifiers.is_empty() {
            f.location.push_str(content);
        } else if let Some((_, Some(last))) = f.qualifiers.last_mut() {
            last.push(' ');
            last.push_str(content);
        }
    }
}

/// Whether a value that started with `"` has reached its closing quote. Inside GenBank values a
/// literal quote is written as two quotes, so count them.
fn quote_closed(value: &str) -> bool {
    let quotes = value.bytes().filter(|&b| b == b'"').count();
    quotes >= 2 && quotes % 2 == 0
}

/// Strip the outer quotes and turn doubled quotes back into single ones.
fn finish_quoted(value: &str) -> String {
    let inner = value.trim();
    let inner = inner.strip_prefix('"').unwrap_or(inner);
    let inner = inner.strip_suffix('"').unwrap_or(inner);
    inner.replace("\"\"", "\"")
}

fn to_qualifiers(raw: Vec<(String, Option<String>)>) -> Vec<Qualifier> {
    let mut out: Vec<Qualifier> = Vec::new();
    for (key, value) in raw {
        let value = match value {
            None => QualValue::Empty,
            Some(v) => match (NUMERIC_QUALIFIERS.contains(&key.as_str()), v.trim().parse::<i64>()) {
                (true, Ok(n)) => QualValue::Int(n),
                _ => QualValue::Text(v),
            },
        };
        match out.iter_mut().find(|q| q.key == key) {
            Some(q) => q.values.push(value),
            None => out.push(Qualifier { key, values: vec![value] }),
        }
    }
    out
}

/// The name to show for a feature: the first of /label, /gene, /product, /standard_name,
/// /locus_tag, /note (shortened) or the feature key.
fn feature_label(feature: &Feature) -> String {
    for key in ["label", "gene", "product", "standard_name", "locus_tag", "allele", "bound_moiety"] {
        if let Some(v) = feature.qualifier_text(key).filter(|v| !v.trim().is_empty()) {
            return v.trim().to_string();
        }
    }
    if let Some(note) = feature.qualifier_text("note").filter(|v| !v.trim().is_empty()) {
        let short: String = note.split_whitespace().take(4).collect::<Vec<_>>().join(" ");
        return short;
    }
    feature.kind.clone()
}

/// Pick up colours that other tools store as qualifiers (ApE, Benchling, SnapGene exports).
fn apply_color_qualifiers(feature: &mut Feature) {
    let forward = feature.qualifier_text("ApEinfo_fwdcolor");
    let reverse = feature.qualifier_text("ApEinfo_revcolor");
    let plain = feature.qualifier_text("color").or_else(|| feature.qualifier_text("colour"));
    let chosen = if feature.strand.is_reverse() { reverse.or(forward) } else { forward.or(reverse) }.or(plain);
    if let Some(color) = chosen.and_then(|c| Color::parse_hex(&c)) {
        feature.set_color(Some(color));
    }
}

fn parse_locus(line: &str, doc: &mut Document, declared_len: &mut Option<usize>) {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    // Find the length: the number in front of "bp" or "aa".
    let unit_index = tokens.iter().position(|t| *t == "bp" || *t == "aa");
    match unit_index {
        Some(ui) if ui >= 2 => {
            *declared_len = tokens[ui - 1].parse().ok();
            doc.name = tokens[1..ui - 1].join("_");
            if tokens[ui] == "aa" {
                doc.molecule = Molecule::Protein;
            }
        }
        _ => {
            doc.name = tokens.get(1).copied().unwrap_or("").to_string();
        }
    }
    for token in tokens.iter().skip(unit_index.map(|u| u + 1).unwrap_or(2)) {
        let lower = token.to_ascii_lowercase();
        if lower == "circular" {
            doc.topology = Topology::Circular;
        } else if lower == "linear" {
            doc.topology = Topology::Linear;
        } else if lower.contains("rna") && doc.molecule != Molecule::Protein {
            doc.molecule = Molecule::Rna;
            doc.strandedness = Strandedness::Single;
        } else if lower.starts_with("ss-") {
            doc.strandedness = Strandedness::Single;
        } else if lower.len() == 3 && token.chars().all(|c| c.is_ascii_uppercase()) && lower != "dna" && lower != "rna" {
            doc.meta.division = Some(token.to_string());
        } else if let Some(date) = genbank_date_to_iso(token) {
            doc.meta.modified = Some(date);
        }
    }
}

const MONTHS: [&str; 12] = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"];

/// `30-MAY-2013` -> `2013-05-30`.
fn genbank_date_to_iso(text: &str) -> Option<String> {
    let mut parts = text.split('-');
    let day: u32 = parts.next()?.parse().ok()?;
    let month = parts.next()?.to_ascii_uppercase();
    let year: u32 = parts.next()?.parse().ok()?;
    let month_index = MONTHS.iter().position(|m| *m == month)? + 1;
    Some(format!("{year:04}-{month_index:02}-{day:02}"))
}

/// `2013-05-30` -> `30-MAY-2013`.
fn iso_to_genbank_date(text: &str) -> Option<String> {
    let mut parts = text.split('-');
    let year: u32 = parts.next()?.parse().ok()?;
    let month: usize = parts.next()?.parse().ok()?;
    let day: u32 = parts.next()?.parse().ok()?;
    Some(format!("{day:02}-{}-{year:04}", MONTHS.get(month.checked_sub(1)?)?))
}

// ---------------------------------------------------------------------------------------------
// Locations

/// One piece of a location in biological order, with its own strand.
#[derive(Debug, Clone, Copy)]
struct Piece {
    first: usize,
    last: usize,
    reverse: bool,
}

/// Parse a GenBank location into top-strand-ordered spans and an overall strand.
///
/// Supports `a..b`, `a`, `a^b`, `<a..>b`, `complement(...)`, `join(...)`, `order(...)`, and
/// joins that run through the origin of a circular sequence. Pieces on other records
/// (`X12345.1:1..20`) are skipped.
pub fn parse_location(text: &str, seq_len: usize, circular: bool) -> Result<(Vec<Span>, Strand)> {
    let compact: String = text.chars().filter(|c| !c.is_whitespace()).collect();
    let mut parser = LocationParser { src: compact.as_bytes(), pos: 0 };
    let pieces = parser.parse()?;
    if parser.pos != parser.src.len() {
        return Err(parse_err(format!("trailing text in location '{text}'")));
    }
    if pieces.is_empty() {
        return Ok((Vec::new(), Strand::Forward));
    }

    let reverse = pieces.iter().filter(|p| p.reverse).count() * 2 > pieces.len();
    let strand = if reverse { Strand::Reverse } else { Strand::Forward };

    // Biological order -> top-strand order: reverse-strand features list their pieces 3'->5'
    // along the top strand, so flip them.
    let mut ordered = pieces;
    if reverse {
        ordered.reverse();
    }
    let mut spans = Vec::with_capacity(ordered.len());
    for piece in ordered {
        let first = piece.first.clamp(1, seq_len.max(1));
        let last = piece.last.clamp(1, seq_len.max(1));
        let span = if first <= last {
            Span::new(first - 1, last)
        } else if circular {
            Span::new(first - 1, last + seq_len)
        } else {
            Span::new(last - 1, first)
        };
        spans.push(span);
    }
    Ok((spans, strand))
}

struct LocationParser<'a> {
    src: &'a [u8],
    pos: usize,
}

impl LocationParser<'_> {
    fn peek_word(&self, word: &str) -> bool {
        self.src[self.pos..].starts_with(word.as_bytes())
    }

    fn expect(&mut self, byte: u8) -> Result<()> {
        if self.src.get(self.pos) == Some(&byte) {
            self.pos += 1;
            Ok(())
        } else {
            Err(parse_err(format!("expected '{}' at position {} of location", byte as char, self.pos + 1)))
        }
    }

    fn number(&mut self) -> Result<usize> {
        // Partial markers '<' and '>' are accepted and ignored.
        while matches!(self.src.get(self.pos), Some(b'<') | Some(b'>')) {
            self.pos += 1;
        }
        let start = self.pos;
        while self.src.get(self.pos).is_some_and(|b| b.is_ascii_digit()) {
            self.pos += 1;
        }
        std::str::from_utf8(&self.src[start..self.pos])
            .ok()
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| parse_err(format!("expected a number at position {} of location", start + 1)))
    }

    /// location := complement(location) | join(list) | order(list) | range
    fn parse(&mut self) -> Result<Vec<Piece>> {
        for word in ["complement(", "join(", "order(", "bond("] {
            if self.peek_word(word) {
                self.pos += word.len();
                let mut pieces = Vec::new();
                loop {
                    pieces.extend(self.parse()?);
                    if self.src.get(self.pos) == Some(&b',') {
                        self.pos += 1;
                    } else {
                        break;
                    }
                }
                self.expect(b')')?;
                if word == "complement(" {
                    // complement() flips the strand and the biological order of what it wraps.
                    pieces.reverse();
                    for p in &mut pieces {
                        p.reverse = !p.reverse;
                    }
                }
                return Ok(pieces);
            }
        }
        // A piece on another record ("X12345.1:10..20") cannot be placed here; skip it.
        let rest = &self.src[self.pos..];
        let end = rest.iter().position(|&b| b == b',' || b == b')').unwrap_or(rest.len());
        if rest[..end].contains(&b':') {
            self.pos += end;
            return Ok(Vec::new());
        }
        let first = self.number()?;
        let last = if self.peek_word("..") {
            self.pos += 2;
            self.number()?
        } else if self.src.get(self.pos) == Some(&b'^') {
            // Between two bases: record it as the single base after the gap.
            self.pos += 1;
            let after = self.number()?;
            return Ok(vec![Piece { first: after, last: after, reverse: false }]);
        } else if self.src.get(self.pos) == Some(&b'.') {
            // Old "a.b" (one base somewhere in a..b): take the whole range.
            self.pos += 1;
            self.number()?
        } else {
            first
        };
        Ok(vec![Piece { first, last, reverse: false }])
    }
}

/// Format spans as a GenBank location, using `join` for several pieces or a span through the
/// origin, and `complement` for reverse-strand features.
pub fn format_location(spans: &[Span], strand: Strand, seq_len: usize) -> String {
    let mut parts = Vec::new();
    for span in spans {
        for (start, end) in span.pieces(seq_len) {
            if end == start + 1 {
                parts.push(format!("{end}"));
            } else {
                parts.push(format!("{}..{}", start + 1, end));
            }
        }
    }
    let inner = if parts.len() == 1 { parts.remove(0) } else { format!("join({})", parts.join(",")) };
    if strand.is_reverse() { format!("complement({inner})") } else { inner }
}

// ---------------------------------------------------------------------------------------------
// Writing

/// Write one document as a GenBank record.
pub fn write(doc: &Document) -> String {
    let mut out = String::new();
    let seq_len = doc.seq.len();
    let name: String = doc.name.chars().map(|c| if c.is_whitespace() { '_' } else { c }).collect::<String>();
    let name = if name.is_empty() { "unnamed".to_string() } else { name };
    // Strand prefix (columns 45-47) and molecule type (columns 48-53).
    let (unit, strand_prefix, mol) = match doc.molecule {
        Molecule::Protein => ("aa", "", ""),
        Molecule::Rna => ("bp", "", "RNA"),
        Molecule::Dna if doc.strandedness == Strandedness::Single => ("bp", "ss-", "DNA"),
        Molecule::Dna => ("bp", "", "DNA"),
    };
    let division = doc.meta.division.clone().unwrap_or_else(|| "SYN".to_string());
    let date = doc.meta.modified.as_deref().or(doc.meta.created.as_deref()).and_then(iso_to_genbank_date).unwrap_or_else(|| "01-JAN-1980".to_string());

    // LOCUS line in the NCBI layout: name in columns 13-28, length right-aligned in 30-40,
    // unit in 42-43, molecule in 45-53, topology in 56-63, division in 65-67, date in 69-79.
    out.push_str(&format!(
        "LOCUS       {:<16} {:>11} {} {:<3}{:<6}  {:<8} {} {}\n",
        name,
        seq_len,
        unit,
        strand_prefix,
        mol,
        doc.topology.name(),
        division,
        date
    ));

    let definition = if doc.meta.description.trim().is_empty() { ".".to_string() } else { one_line(&doc.meta.description) };
    push_wrapped(&mut out, "DEFINITION", &definition);
    push_wrapped(&mut out, "ACCESSION", doc.meta.accession.as_deref().unwrap_or("."));
    push_wrapped(&mut out, "VERSION", doc.meta.version.as_deref().unwrap_or("."));
    push_wrapped(&mut out, "KEYWORDS", doc.meta.keywords.as_deref().unwrap_or("."));
    let organism = doc.meta.organism.clone().unwrap_or_else(|| ".".to_string());
    push_wrapped(&mut out, "SOURCE", &organism);
    push_wrapped(&mut out, "  ORGANISM", &organism);
    for (ii, r) in doc.meta.references.iter().enumerate() {
        push_wrapped(&mut out, "REFERENCE", &format!("{}  (bases 1 to {seq_len})", ii + 1));
        push_wrapped(&mut out, "  AUTHORS", if r.authors.is_empty() { "." } else { &r.authors });
        push_wrapped(&mut out, "  TITLE", if r.title.is_empty() { "Direct Submission" } else { &r.title });
        push_wrapped(&mut out, "  JOURNAL", if r.journal.is_empty() { "Unpublished" } else { &r.journal });
        if let Some(pm) = &r.pubmed {
            push_wrapped(&mut out, "   PUBMED", pm);
        }
    }
    for comment in &doc.meta.comments {
        let mut first = true;
        for line in comment.lines() {
            push_wrapped(&mut out, if first { "COMMENT" } else { "" }, line);
            first = false;
        }
    }

    out.push_str("FEATURES             Location/Qualifiers\n");
    // A source feature first, as GenBank expects.
    let source_quals = doc.raw.genbank_source.clone().unwrap_or_else(|| {
        let mol_type = match doc.molecule {
            Molecule::Rna => "other RNA",
            _ => "other DNA",
        };
        vec![
            Qualifier { key: "organism".into(), values: vec![QualValue::Text(organism.clone())] },
            Qualifier { key: "mol_type".into(), values: vec![QualValue::Text(mol_type.into())] },
        ]
    });
    if seq_len > 0 {
        push_feature(&mut out, "source", &format!("1..{seq_len}"), &source_quals);
    }
    for feature in &doc.features {
        let spans: Vec<Span> = feature.segments.iter().filter(|s| !s.gap).map(|s| s.span).collect();
        if spans.is_empty() {
            continue;
        }
        let location = format_location(&spans, feature.strand, seq_len);
        let mut quals = Vec::new();
        // The feature name travels as /label unless a /label already says the same thing.
        let has_label = feature.qualifier_text("label").is_some_and(|l| l == feature.name);
        if !feature.name.is_empty() && !has_label {
            quals.push(Qualifier { key: "label".into(), values: vec![QualValue::Text(feature.name.clone())] });
        }
        for q in &feature.qualifiers {
            if q.key == "label" && !has_label {
                continue;
            }
            quals.push(q.clone());
        }
        if let Some(color) = feature.color()
            && feature.qualifier("ApEinfo_fwdcolor").is_none()
        {
            quals.push(Qualifier { key: "ApEinfo_fwdcolor".into(), values: vec![QualValue::Text(color.to_hex())] });
            quals.push(Qualifier { key: "ApEinfo_revcolor".into(), values: vec![QualValue::Text(color.to_hex())] });
        }
        push_feature(&mut out, &feature.kind, &location, &quals);
    }
    for primer in &doc.primers {
        for site in &primer.sites {
            if site.extra.attr("simplified") == Some("1") {
                continue;
            }
            let location = format_location(&[site.span], site.strand, seq_len);
            let quals = vec![
                Qualifier { key: "label".into(), values: vec![QualValue::Text(primer.name.clone())] },
                Qualifier { key: "note".into(), values: vec![QualValue::Text(format!("sequence: {}", primer.seq))] },
            ];
            push_feature(&mut out, "primer_bind", &location, &quals);
        }
    }

    out.push_str("ORIGIN\n");
    let lower = doc.seq.to_ascii_lowercase();
    for (ii, chunk) in lower.chunks(60).enumerate() {
        out.push_str(&format!("{:>9}", ii * 60 + 1));
        for block in chunk.chunks(10) {
            out.push(' ');
            out.push_str(std::str::from_utf8(block).unwrap_or(""));
        }
        out.push('\n');
    }
    out.push_str("//\n");
    out
}

fn one_line(text: &str) -> String {
    text.split_whitespace().collect::<Vec<_>>().join(" ")
}

/// Write a header keyword with its value wrapped at 79 columns (continuations indented 12).
fn push_wrapped(out: &mut String, keyword: &str, value: &str) {
    let lines = wrap_words(value, 79 - 12);
    for (ii, line) in lines.iter().enumerate() {
        let key = if ii == 0 { keyword } else { "" };
        out.push_str(&format!("{key:<12}{line}\n"));
    }
    if lines.is_empty() {
        out.push_str(&format!("{keyword:<12}\n"));
    }
}

/// Greedy word wrap; words longer than the width are hard-broken.
fn wrap_words(text: &str, width: usize) -> Vec<String> {
    let mut lines = Vec::new();
    let mut current = String::new();
    for word in text.split_whitespace() {
        let mut word = word.to_string();
        while word.chars().count() > width {
            if !current.is_empty() {
                lines.push(std::mem::take(&mut current));
            }
            let split: String = word.chars().take(width).collect();
            word = word.chars().skip(width).collect();
            lines.push(split);
        }
        if current.is_empty() {
            current = word;
        } else if current.chars().count() + 1 + word.chars().count() <= width {
            current.push(' ');
            current.push_str(&word);
        } else {
            lines.push(std::mem::replace(&mut current, word));
        }
    }
    if !current.is_empty() {
        lines.push(current);
    }
    lines
}

fn push_feature(out: &mut String, key: &str, location: &str, quals: &[Qualifier]) {
    // Location: break after commas so lines stay within 79 columns.
    let mut line = format!("     {key:<16}");
    let mut first = true;
    for part in split_keep_commas(location) {
        if !first && line.len() + part.len() > 79 {
            out.push_str(&line);
            out.push('\n');
            line = " ".repeat(21);
        }
        line.push_str(&part);
        first = false;
    }
    out.push_str(&line);
    out.push('\n');

    for q in quals {
        for value in &q.values {
            let text = match value {
                QualValue::Empty => {
                    out.push_str(&format!("{}/{}\n", " ".repeat(21), q.key));
                    continue;
                }
                QualValue::Int(n) => format!("/{}={}", q.key, n),
                other => {
                    let t = other.as_text();
                    if NUMERIC_QUALIFIERS.contains(&q.key.as_str()) && t.trim().parse::<i64>().is_ok() {
                        format!("/{}={}", q.key, t.trim())
                    } else {
                        format!("/{}=\"{}\"", q.key, one_line(&t).replace('"', "\"\""))
                    }
                }
            };
            let width = 79 - 21;
            let wrapped = if q.key == "translation" {
                // Protein text has no spaces; break it into fixed-width lines.
                let chars: Vec<char> = text.chars().collect();
                chars.chunks(width).map(|c| c.iter().collect::<String>()).collect()
            } else {
                wrap_words(&text, width)
            };
            for piece in wrapped {
                out.push_str(&" ".repeat(21));
                out.push_str(&piece);
                out.push('\n');
            }
        }
    }
}

fn split_keep_commas(text: &str) -> Vec<String> {
    let mut parts = Vec::new();
    let mut current = String::new();
    for ch in text.chars() {
        current.push(ch);
        if ch == ',' {
            parts.push(std::mem::take(&mut current));
        }
    }
    if !current.is_empty() {
        parts.push(current);
    }
    parts
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = r#"LOCUS       test_plasmid1           5767 bp    DNA     circular SYN 18-AUG-2026
DEFINITION  A test plasmid
            with a two-line definition.
ACCESSION   test_plasmid1
KEYWORDS    .
SOURCE
  ORGANISM  .
            .
COMMENT     Annotated with pLannotate v1.2.2
FEATURES             Location/Qualifiers
     source          1..5767
                     /organism="synthetic DNA construct"
     promoter        complement(join(5588..5767,1..24))
                     /note="pLannotate"
                     /label="CMV promoter"
     CDS             602..903
                     /label="TcR (fragment)"
                     /codon_start=1
                     /translation="MKSNNALIVILGTVTLDAVGIGLVMPVLPGLLRDIVHSDSIASHYGVLLALYALMQFLCA
                     PVLGALSDRF"
                     /note="a note with ""quotes"" in it"
     misc_feature    join(10..20,30..40)
                     /pseudo
ORIGIN
"#;

    fn sample_with_sequence() -> String {
        let mut text = SAMPLE.to_string();
        let seq = "acgt".repeat(5767 / 4 + 1);
        let seq = &seq[..5767];
        for (ii, chunk) in seq.as_bytes().chunks(60).enumerate() {
            text.push_str(&format!("{:>9}", ii * 60 + 1));
            for block in chunk.chunks(10) {
                text.push(' ');
                text.push_str(std::str::from_utf8(block).unwrap());
            }
            text.push('\n');
        }
        text.push_str("//\n");
        text
    }

    #[test]
    fn reads_plannotate_style_record() {
        let docs = read(&sample_with_sequence()).unwrap();
        let doc = &docs[0];
        assert_eq!(doc.name, "test_plasmid1");
        assert_eq!(doc.len(), 5767);
        assert!(doc.is_circular());
        assert_eq!(doc.meta.description, "A test plasmid with a two-line definition");
        assert_eq!(doc.meta.organism.as_deref(), Some("synthetic DNA construct"));
        assert_eq!(doc.features.len(), 3);

        // The promoter runs through the origin on the reverse strand.
        let promoter = &doc.features[0];
        assert_eq!(promoter.name, "CMV promoter");
        assert_eq!(promoter.strand, Strand::Reverse);
        assert_eq!(promoter.span(doc.len()), Span::new(5587, 5767 + 24));

        let cds = &doc.features[1];
        assert_eq!(cds.qualifier("codon_start").unwrap().values[0], QualValue::Int(1));
        assert!(cds.qualifier_text("translation").unwrap().ends_with("CAPVLGALSDRF"));
        assert_eq!(cds.qualifier_text("note").unwrap(), "a note with \"quotes\" in it");

        let joined = &doc.features[2];
        assert_eq!(joined.segments.len(), 2);
        assert_eq!(joined.qualifier("pseudo").unwrap().values[0], QualValue::Empty);
    }

    #[test]
    fn write_read_round_trip() {
        let doc = read(&sample_with_sequence()).unwrap().remove(0);
        let text = write(&doc);
        let back = read(&text).unwrap().remove(0);
        assert_eq!(back.seq, doc.seq);
        assert_eq!(back.features.len(), doc.features.len());
        for (a, b) in doc.features.iter().zip(&back.features) {
            assert_eq!(a.name, b.name);
            assert_eq!(a.strand, b.strand);
            assert_eq!(a.span(doc.len()), b.span(doc.len()));
        }
        assert!(text.lines().all(|l| l.len() <= 80), "a line is longer than 80 columns");
    }

    #[test]
    fn location_forms() {
        let (spans, strand) = parse_location("complement(join(5588..5767,1..24))", 5767, true).unwrap();
        assert_eq!(strand, Strand::Reverse);
        assert_eq!(spans, vec![Span::new(5587, 5767), Span::new(0, 24)]);

        let (spans, strand) = parse_location("join(complement(30..40),complement(10..20))", 100, false).unwrap();
        assert_eq!(strand, Strand::Reverse);
        assert_eq!(spans, vec![Span::new(9, 20), Span::new(29, 40)]);

        let (spans, _) = parse_location("<1..>50", 100, false).unwrap();
        assert_eq!(spans, vec![Span::new(0, 50)]);

        let (spans, _) = parse_location("42", 100, false).unwrap();
        assert_eq!(spans, vec![Span::new(41, 42)]);

        assert_eq!(format_location(&[Span::new(5587, 5767 + 24)], Strand::Reverse, 5767), "complement(join(5588..5767,1..24))");
    }
}
