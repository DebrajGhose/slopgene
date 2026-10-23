//! SlopGene's JSON document schema, used for `--json` output and as a file format.
//!
//! The schema is meant to be read and written by scripts and agents, so it uses the same
//! 1-based inclusive coordinates as the rest of the user-facing surface. A span that runs
//! through the origin of a circular sequence has `start > end`.
//!
//! ```json
//! { "schema": "slopgene.document/1", "name": "pUC19", "molecule": "dna",
//!   "topology": "circular", "length": 2686, "sequence": "TCGC...",
//!   "features": [ { "name": "AmpR", "type": "CDS", "strand": "reverse",
//!                   "start": 1626, "end": 2486, "length": 861, "color": "#ccff33",
//!                   "segments": [ { "start": 1626, "end": 2486 } ],
//!                   "qualifiers": { "gene": ["bla"] } } ],
//!   "primers": [ { "name": "M13 fwd", "sequence": "GTAAAACGACGGCCAGT",
//!                  "sites": [ { "start": 1, "end": 17, "strand": "forward" } ] } ] }
//! ```

use crate::error::{Error, Result};
use crate::model::*;
use crate::span::Span;
use serde_json::{Map, Value, json};

pub const SCHEMA: &str = "slopgene.document/1";

fn span_json(span: Span, seq_len: usize) -> (usize, usize) {
    span.to_1based(seq_len)
}

fn overhang_json(o: Overhang) -> Value {
    match o {
        Overhang::Blunt => json!({ "type": "blunt", "length": 0 }),
        Overhang::FivePrime(n) => json!({ "type": "5'", "length": n }),
        Overhang::ThreePrime(n) => json!({ "type": "3'", "length": n }),
    }
}

/// JSON for one feature.
pub fn feature_to_value(feature: &Feature, seq_len: usize) -> Value {
    let (start, end) = span_json(feature.span(seq_len), seq_len);
    let mut obj = Map::new();
    obj.insert("name".into(), json!(feature.name));
    obj.insert("type".into(), json!(feature.kind));
    obj.insert("strand".into(), json!(feature.strand.name()));
    obj.insert("start".into(), json!(start));
    obj.insert("end".into(), json!(end));
    obj.insert("length".into(), json!(feature.len()));
    if let Some(color) = feature.color() {
        obj.insert("color".into(), json!(color.to_hex()));
    }
    let segments: Vec<Value> = feature
        .segments
        .iter()
        .map(|s| {
            let (a, b) = span_json(s.span, seq_len);
            let mut seg = Map::new();
            seg.insert("start".into(), json!(a));
            seg.insert("end".into(), json!(b));
            if s.gap {
                seg.insert("gap".into(), json!(true));
            }
            if let Some(name) = &s.name {
                seg.insert("name".into(), json!(name));
            }
            if let Some(color) = s.color {
                seg.insert("color".into(), json!(color.to_hex()));
            }
            Value::Object(seg)
        })
        .collect();
    obj.insert("segments".into(), Value::Array(segments));
    let mut quals = Map::new();
    for q in &feature.qualifiers {
        let values: Vec<Value> = q
            .values
            .iter()
            .map(|v| match v {
                QualValue::Int(n) => json!(n),
                QualValue::Empty => Value::Null,
                other => json!(other.as_text()),
            })
            .collect();
        quals.insert(q.key.clone(), Value::Array(values));
    }
    obj.insert("qualifiers".into(), Value::Object(quals));
    Value::Object(obj)
}

/// JSON for one primer.
pub fn primer_to_value(primer: &Primer, seq_len: usize) -> Value {
    let sites: Vec<Value> = primer
        .sites
        .iter()
        .filter(|s| s.extra.attr("simplified") != Some("1"))
        .map(|s| {
            let (a, b) = span_json(s.span, seq_len);
            json!({ "start": a, "end": b, "strand": s.strand.name(), "tm": s.tm })
        })
        .collect();
    json!({
        "name": primer.name,
        "sequence": primer.seq,
        "length": primer.seq.len(),
        "description": primer.description,
        "sites": sites,
    })
}

/// JSON for a whole document. `include_sequence` false leaves out the bases (for summaries).
pub fn document_to_value(doc: &Document, include_sequence: bool) -> Value {
    let seq_len = doc.len();
    let mut obj = Map::new();
    obj.insert("schema".into(), json!(SCHEMA));
    obj.insert("name".into(), json!(doc.name));
    obj.insert("molecule".into(), json!(doc.molecule.name()));
    obj.insert("topology".into(), json!(doc.topology.name()));
    obj.insert("strandedness".into(), json!(if doc.strandedness == Strandedness::Double { "double" } else { "single" }));
    obj.insert("length".into(), json!(seq_len));
    if include_sequence {
        obj.insert("sequence".into(), json!(doc.seq_str()));
    }
    obj.insert("description".into(), json!(doc.meta.description));
    if doc.molecule == Molecule::Dna {
        obj.insert("methylation".into(), json!({ "dam": doc.methylation.dam, "dcm": doc.methylation.dcm, "ecoki": doc.methylation.ecoki }));
        if !doc.is_circular() {
            obj.insert("ends".into(), json!({ "left": overhang_json(doc.ends.left), "right": overhang_json(doc.ends.right) }));
        }
    }
    let meta = &doc.meta;
    let mut m = Map::new();
    let mut put = |key: &str, value: &Option<String>| {
        if let Some(v) = value {
            m.insert(key.into(), json!(v));
        }
    };
    put("accession", &meta.accession);
    put("version", &meta.version);
    put("organism", &meta.organism);
    put("keywords", &meta.keywords);
    put("created", &meta.created);
    put("modified", &meta.modified);
    put("created_by", &meta.created_by);
    put("division", &meta.division);
    if let Some(s) = meta.synthetic {
        m.insert("synthetic".into(), json!(s));
    }
    if !meta.comments.is_empty() {
        m.insert("comments".into(), json!(meta.comments));
    }
    if !meta.references.is_empty() {
        let refs: Vec<Value> =
            meta.references.iter().map(|r| json!({ "authors": r.authors, "title": r.title, "journal": r.journal, "pubmed": r.pubmed })).collect();
        m.insert("references".into(), Value::Array(refs));
    }
    obj.insert("meta".into(), Value::Object(m));
    obj.insert("features".into(), Value::Array(doc.features.iter().map(|f| feature_to_value(f, seq_len)).collect()));
    obj.insert("primers".into(), Value::Array(doc.primers.iter().map(|p| primer_to_value(p, seq_len)).collect()));
    Value::Object(obj)
}

fn bad(msg: impl Into<String>) -> Error {
    Error::Parse(format!("JSON document: {}", msg.into()))
}

fn get_str<'a>(obj: &'a Value, key: &str) -> Option<&'a str> {
    obj.get(key).and_then(Value::as_str)
}

fn get_usize(obj: &Value, key: &str) -> Option<usize> {
    obj.get(key).and_then(Value::as_u64).map(|n| n as usize)
}

fn read_span(obj: &Value, seq_len: usize, circular: bool) -> Result<Span> {
    let start = get_usize(obj, "start").ok_or_else(|| bad("missing 'start'"))?;
    let end = get_usize(obj, "end").ok_or_else(|| bad("missing 'end'"))?;
    Span::from_1based(start, end, seq_len, circular)
}

/// Read documents from JSON: a single document object, an array of them, or an object with a
/// `documents` array.
pub fn documents_from_value(value: &Value) -> Result<Vec<Document>> {
    match value {
        Value::Array(items) => items.iter().map(document_from_value).collect(),
        Value::Object(obj) if obj.contains_key("documents") => documents_from_value(obj.get("documents").unwrap()),
        Value::Object(_) => Ok(vec![document_from_value(value)?]),
        _ => Err(bad("expected an object or an array")),
    }
}

/// Read one document from JSON in the schema written by [`document_to_value`]. Only
/// `sequence` is required.
pub fn document_from_value(value: &Value) -> Result<Document> {
    let seq_text = get_str(value, "sequence").ok_or_else(|| bad("missing 'sequence'"))?;
    let molecule = match get_str(value, "molecule").unwrap_or("dna").to_ascii_lowercase().as_str() {
        "rna" => Molecule::Rna,
        "protein" | "aa" => Molecule::Protein,
        _ => Molecule::Dna,
    };
    let seq = if molecule == Molecule::Protein {
        crate::alphabet::clean_protein(seq_text).map_err(bad)?
    } else {
        crate::alphabet::clean_nucleotides(seq_text).map_err(bad)?
    };
    let topology = get_str(value, "topology").and_then(Topology::parse).unwrap_or(Topology::Linear);
    let mut doc = Document::new_dna(get_str(value, "name").unwrap_or(""), seq, topology);
    doc.molecule = molecule;
    if molecule != Molecule::Dna || get_str(value, "strandedness") == Some("single") {
        doc.strandedness = Strandedness::Single;
    }
    doc.meta.description = get_str(value, "description").unwrap_or("").to_string();
    if let Some(m) = value.get("meta") {
        let text = |key: &str| get_str(m, key).map(str::to_string);
        doc.meta.accession = text("accession");
        doc.meta.version = text("version");
        doc.meta.organism = text("organism");
        doc.meta.keywords = text("keywords");
        doc.meta.created = text("created");
        doc.meta.modified = text("modified");
        doc.meta.created_by = text("created_by");
        doc.meta.division = text("division");
        doc.meta.synthetic = m.get("synthetic").and_then(Value::as_bool);
        if let Some(comments) = m.get("comments").and_then(Value::as_array) {
            doc.meta.comments = comments.iter().filter_map(|c| c.as_str().map(str::to_string)).collect();
        }
    }
    if let Some(meth) = value.get("methylation") {
        let flag = |k: &str| meth.get(k).and_then(Value::as_bool).unwrap_or(false);
        doc.methylation = Methylation { dam: flag("dam"), dcm: flag("dcm"), ecoki: flag("ecoki") };
    }
    if let Some(ends) = value.get("ends") {
        let read_end = |v: Option<&Value>| -> Overhang {
            let Some(v) = v else { return Overhang::Blunt };
            let n = get_usize(v, "length").unwrap_or(0);
            match get_str(v, "type").unwrap_or("blunt") {
                "5'" | "5" | "5prime" | "five_prime" if n > 0 => Overhang::FivePrime(n),
                "3'" | "3" | "3prime" | "three_prime" if n > 0 => Overhang::ThreePrime(n),
                _ => Overhang::Blunt,
            }
        };
        doc.ends.left = read_end(ends.get("left"));
        doc.ends.right = read_end(ends.get("right"));
    }

    let seq_len = doc.len();
    let circular = doc.is_circular();
    if let Some(features) = value.get("features").and_then(Value::as_array) {
        for (ii, f) in features.iter().enumerate() {
            let strand = get_str(f, "strand").and_then(Strand::parse).unwrap_or(Strand::Forward);
            let kind = get_str(f, "type").unwrap_or("misc_feature");
            let name = get_str(f, "name").unwrap_or(kind);
            let mut feature = Feature::new(name, kind, Span::default(), strand);
            feature.segments.clear();
            match f.get("segments").and_then(Value::as_array).filter(|s| !s.is_empty()) {
                Some(segments) => {
                    for s in segments {
                        let mut segment = Segment::new(read_span(s, seq_len, circular).map_err(|e| bad(format!("feature {}: {e}", ii + 1)))?);
                        segment.gap = s.get("gap").and_then(Value::as_bool).unwrap_or(false);
                        segment.name = get_str(s, "name").map(str::to_string);
                        segment.color = get_str(s, "color").and_then(Color::parse_hex);
                        feature.segments.push(segment);
                    }
                }
                None => {
                    let span = read_span(f, seq_len, circular).map_err(|e| bad(format!("feature {}: {e}", ii + 1)))?;
                    feature.segments.push(Segment::new(span));
                }
            }
            if let Some(color) = get_str(f, "color").and_then(Color::parse_hex) {
                for s in feature.segments.iter_mut().filter(|s| s.color.is_none()) {
                    s.color = Some(color);
                }
            }
            if let Some(quals) = f.get("qualifiers").and_then(Value::as_object) {
                for (key, values) in quals {
                    let list = match values {
                        Value::Array(items) => items.clone(),
                        other => vec![other.clone()],
                    };
                    let values = list
                        .iter()
                        .map(|v| match v {
                            Value::Null => QualValue::Empty,
                            Value::Number(n) if n.is_i64() => QualValue::Int(n.as_i64().unwrap()),
                            Value::String(s) => QualValue::Text(s.clone()),
                            other => QualValue::Text(other.to_string()),
                        })
                        .collect();
                    feature.qualifiers.push(Qualifier { key: key.clone(), values });
                }
            }
            feature.normalize(seq_len, circular);
            doc.features.push(feature);
        }
    }
    if let Some(primers) = value.get("primers").and_then(Value::as_array) {
        for p in primers {
            let mut primer = Primer::new(get_str(p, "name").unwrap_or(""), get_str(p, "sequence").unwrap_or(""));
            primer.description = get_str(p, "description").unwrap_or("").to_string();
            if let Some(sites) = p.get("sites").and_then(Value::as_array) {
                for s in sites {
                    primer.sites.push(BindingSite {
                        span: read_span(s, seq_len, circular)?,
                        strand: get_str(s, "strand").and_then(Strand::parse).unwrap_or(Strand::Forward),
                        tm: s.get("tm").and_then(Value::as_f64),
                        extra: XmlExtra::default(),
                    });
                }
            }
            doc.primers.push(primer);
        }
    }
    Ok(doc)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn json_round_trip() {
        let mut doc = Document::new_dna("demo", b"ATGCATGCATGCATGC", Topology::Circular);
        let mut f = Feature::new("wrap", "misc_feature", Span::new(12, 20), Strand::Reverse);
        f.set_color(Some(Color::rgb(10, 20, 30)));
        f.set_qualifier("note", QualValue::Text("hello".into()));
        f.set_qualifier("codon_start", QualValue::Int(1));
        doc.features.push(f);
        let value = document_to_value(&doc, true);
        assert_eq!(value["features"][0]["start"], 13);
        assert_eq!(value["features"][0]["end"], 4);
        let back = document_from_value(&value).unwrap();
        assert_eq!(back.seq, doc.seq);
        assert_eq!(back.features[0].span(16), Span::new(12, 20));
        assert_eq!(back.features[0].strand, Strand::Reverse);
        assert_eq!(back.features[0].color(), Some(Color::rgb(10, 20, 30)));
        assert_eq!(back.features[0].qualifier("codon_start").unwrap().values[0], QualValue::Int(1));
    }
}
