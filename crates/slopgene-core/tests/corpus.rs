//! Round-trip tests over real files.
//!
//! The corpus is not part of the repository (SnapGene's sample files are not ours to
//! redistribute, and lab files are private). Point `SLOPGENE_CORPUS` at one or more folders,
//! separated by `;` on Windows or `:` elsewhere. When it is unset, the SnapGene installation
//! folder is used if present. With no corpus the test passes without checking anything.
//!
//! Run with `cargo test -p slopgene-core --test corpus -- --nocapture` to see the report.

use slopgene_core::io::{self, Format, snapgene};
use slopgene_core::xml::{self, canonical};
use slopgene_core::*;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

fn corpus_dirs() -> Vec<PathBuf> {
    if let Ok(list) = std::env::var("SLOPGENE_CORPUS") {
        let sep = if cfg!(windows) { ';' } else { ':' };
        return list.split(sep).filter(|s| !s.is_empty()).map(PathBuf::from).collect();
    }
    let defaults = ["C:\\Program Files\\SnapGene", "/Applications/SnapGene.app/Contents/Resources", "/opt/gslbiotech/snapgene"];
    defaults.iter().map(PathBuf::from).filter(|p| p.exists()).collect()
}

fn collect(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(entries) = std::fs::read_dir(dir) else { return };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            collect(&path, out);
        } else if let Some(ext) = path.extension().and_then(|e| e.to_str())
            && matches!(ext.to_ascii_lowercase().as_str(), "dna" | "prot" | "rna" | "gb" | "gbk" | "genbank")
        {
            out.push(path);
        }
    }
}

/// Compare the model of two documents and describe the first difference.
fn compare_models(a: &Document, b: &Document) -> Option<String> {
    if a.seq != b.seq {
        return Some("sequence differs".into());
    }
    if a.topology != b.topology || a.strandedness != b.strandedness || a.molecule != b.molecule {
        return Some("topology, strandedness or molecule differs".into());
    }
    if a.features.len() != b.features.len() {
        return Some(format!("{} features vs {}", a.features.len(), b.features.len()));
    }
    for (fa, fb) in a.features.iter().zip(&b.features) {
        let same = fa.name == fb.name
            && fa.kind == fb.kind
            && fa.strand == fb.strand
            && fa.segments.len() == fb.segments.len()
            && fa.segments.iter().zip(&fb.segments).all(|(x, y)| x.span == y.span && x.color == y.color && x.gap == y.gap)
            && fa.qualifiers == fb.qualifiers;
        if !same {
            return Some(format!("feature '{}' differs:\n    before: {fa:?}\n    after:  {fb:?}", fa.name));
        }
    }
    if a.primers.len() != b.primers.len() {
        return Some(format!("{} primers vs {}", a.primers.len(), b.primers.len()));
    }
    for (pa, pb) in a.primers.iter().zip(&b.primers) {
        if pa.name != pb.name || pa.seq != pb.seq || pa.sites.len() != pb.sites.len() {
            return Some(format!("primer '{}' differs", pa.name));
        }
        for (sa, sb) in pa.sites.iter().zip(&pb.sites) {
            if sa.span != sb.span || sa.strand != sb.strand {
                return Some(format!("primer '{}' binding site differs", pa.name));
            }
        }
    }
    if a.meta.description != b.meta.description {
        return Some("description differs".into());
    }
    if a.ends != b.ends {
        return Some("ends differ".into());
    }
    None
}

/// For a primer site, how many of the annealed bases match the primer (as read on its strand).
fn site_match(doc: &Document, primer: &Primer, site: &BindingSite) -> (usize, usize) {
    let template = site.span.extract(&doc.seq).to_ascii_uppercase();
    let primer_seq = primer.seq.to_ascii_uppercase().into_bytes();
    // The annealed part is the 3' end of the primer.
    let probe: Vec<u8> = if site.strand == Strand::Reverse {
        let rc = alphabet::revcomp(&primer_seq);
        rc[..template.len().min(rc.len())].to_vec()
    } else {
        primer_seq[primer_seq.len().saturating_sub(template.len())..].to_vec()
    };
    let matches = probe.iter().zip(&template).filter(|(p, t)| alphabet::iupac_matches(**p, **t) || alphabet::iupac_matches(**t, **p)).count();
    (matches, template.len())
}

/// Insert bases into every SnapGene file, write it, read it back and check the result is
/// consistent: features away from the edit still cover the same bases, stale history and
/// derived packets are gone, and every recomputed primer site matches the new template.
#[test]
fn corpus_edit_round_trip() {
    let dirs = corpus_dirs();
    let mut files = Vec::new();
    for dir in &dirs {
        collect(dir, &mut files);
    }
    let mut checked = 0;
    let mut failures = Vec::new();
    for path in files.iter().filter(|p| p.extension().is_some_and(|e| e.eq_ignore_ascii_case("dna"))) {
        let Ok(doc) = io::read_one(path) else { continue };
        if doc.len() < 100 {
            continue;
        }
        // Insert 7 bases a third of the way in.
        let pos = doc.len() / 3;
        let mut edited = doc.clone();
        edit::insert(&mut edited, pos, b"ACGTACG").unwrap();
        let bytes = io::to_bytes(&edited, Format::SnapGene).unwrap();
        let back = match io::read_bytes(&bytes, Some(Format::SnapGene)) {
            Ok(mut d) => d.remove(0),
            Err(e) => {
                failures.push(format!("{}: re-read failed: {e}", path.display()));
                continue;
            }
        };
        checked += 1;
        if back.seq != edited.seq {
            failures.push(format!("{}: sequence differs after edit", path.display()));
        }
        let seq_len = doc.len();
        for (before, after) in doc.features.iter().zip(&back.features) {
            let span = before.span(seq_len);
            if span.contains(pos % seq_len, seq_len) || span.contains((pos + seq_len - 1) % seq_len, seq_len) {
                continue;
            }
            if before.extract(&doc.seq) != after.extract(&back.seq) {
                failures.push(format!("{}: feature '{}' moved off its bases", path.display(), before.name));
            }
        }
        let kinds: Vec<u8> = snapgene::split_packets(&bytes).unwrap().iter().map(|p| p.kind).collect();
        for stale in [0x02u8, 0x07, 0x0B, 0x1B, 0x22] {
            if kinds.contains(&stale) {
                failures.push(format!("{}: stale packet 0x{stale:02x} kept after an edit", path.display()));
            }
        }
        for primer in &back.primers {
            for site in &primer.sites {
                let (matches, len) = site_match(&back, primer, site);
                if matches != len {
                    failures.push(format!("{}: primer {} site does not match after edit", path.display(), primer.name));
                }
            }
        }
    }
    eprintln!("edited and re-read {checked} SnapGene files");
    for failure in failures.iter().take(20) {
        eprintln!("FAIL {failure}");
    }
    assert!(failures.is_empty(), "{} failures", failures.len());
}

#[test]
fn corpus_round_trip() {
    let dirs = corpus_dirs();
    if dirs.is_empty() {
        eprintln!("no corpus found; set SLOPGENE_CORPUS to run this test");
        return;
    }
    let mut files = Vec::new();
    for dir in &dirs {
        collect(dir, &mut files);
    }
    files.sort();

    let mut read_ok = 0;
    let mut failures: Vec<String> = Vec::new();
    let mut xml_diffs: BTreeMap<String, usize> = BTreeMap::new();
    let mut primer_sites = 0;
    let mut primer_sites_exact = 0;
    let mut primer_mismatch_examples = Vec::new();
    let mut byte_identical = 0;
    let mut snapgene_files = 0;

    for path in &files {
        let format = io::Format::from_path(path).unwrap();
        let docs = match io::read_path(path) {
            Ok(d) => d,
            Err(Error::Unsupported(_)) => continue,
            Err(e) => {
                failures.push(format!("read {}: {e}", path.display()));
                continue;
            }
        };
        read_ok += 1;
        for doc in &docs {
            // Primer coordinates: the stored site should hold the primer's annealed 3' end.
            for primer in &doc.primers {
                for site in &primer.sites {
                    primer_sites += 1;
                    let (matches, len) = site_match(doc, primer, site);
                    if matches == len {
                        primer_sites_exact += 1;
                    } else if primer_mismatch_examples.len() < 12 {
                        // Would the site match better one base to either side?
                        let shifted: Vec<String> = [-1i64, 1]
                            .iter()
                            .map(|&d| {
                                let start = (site.span.start as i64 + d).rem_euclid(doc.len() as i64) as usize;
                                let moved = BindingSite { span: Span::new(start, start + site.span.len()), ..site.clone() };
                                let (m, l) = site_match(doc, primer, &moved);
                                format!("shift {d:+}: {m}/{l}")
                            })
                            .collect();
                        primer_mismatch_examples.push(format!(
                            "{}: primer {} ({} nt) site {} ({}): {matches}/{len} bases match; {}; attrs {:?}",
                            path.file_name().unwrap().to_string_lossy(),
                            primer.name,
                            primer.seq.len(),
                            site.span.format(doc.len()),
                            site.strand.name(),
                            shifted.join(", "),
                            site.extra.attrs
                        ));
                    }
                }
            }

            // Round trip through the same format and check the model survives.
            let bytes = match io::to_bytes(doc, format) {
                Ok(b) => b,
                Err(e) => {
                    failures.push(format!("write {}: {e}", path.display()));
                    continue;
                }
            };
            let back = match io::read_bytes(&bytes, Some(format)) {
                Ok(mut d) => d.remove(0),
                Err(e) => {
                    failures.push(format!("re-read {}: {e}", path.display()));
                    continue;
                }
            };
            if let Some(diff) = compare_models(doc, &back) {
                failures.push(format!("round trip {}: {diff}", path.display()));
            }

            // For SnapGene files, also compare the regenerated XML packets with the originals.
            if format == Format::SnapGene {
                let original_bytes = std::fs::read(path).unwrap();
                snapgene_files += 1;
                if original_bytes == bytes {
                    byte_identical += 1;
                }
                let original = snapgene::split_packets(&original_bytes).unwrap();
                let rewritten = snapgene::split_packets(&bytes).unwrap();
                let kinds_before: Vec<u8> = original.iter().map(|p| p.kind).collect();
                let kinds_after: Vec<u8> = rewritten.iter().map(|p| p.kind).collect();
                if kinds_before != kinds_after {
                    *xml_diffs.entry("packet order or set".into()).or_default() += 1;
                    eprintln!("  {}: packets {:02x?} -> {:02x?}", path.file_name().unwrap().to_string_lossy(), kinds_before, kinds_after);
                }
                for kind in [0x05u8, 0x06, 0x08, 0x0A] {
                    let find = |list: &[RawPacket]| list.iter().find(|p| p.kind == kind).map(|p| p.data.clone());
                    let (Some(a), Some(b)) = (find(&original), find(&rewritten)) else { continue };
                    if a == b {
                        continue;
                    }
                    let parse = |d: &[u8]| xml::parse(std::str::from_utf8(d).unwrap_or(""));
                    match (parse(&a), parse(&b)) {
                        (Ok(x), Ok(y)) if canonical(&x) == canonical(&y) => {
                            *xml_diffs.entry(format!("0x{kind:02x} same DOM, different bytes")).or_default() += 1;
                        }
                        (Ok(x), Ok(y)) => {
                            let key = format!("0x{kind:02x} DOM differs");
                            let count = xml_diffs.entry(key.clone()).or_default();
                            *count += 1;
                            if *count <= 2 {
                                let (cx, cy) = (canonical(&x), canonical(&y));
                                let at = cx.bytes().zip(cy.bytes()).position(|(p, q)| p != q).unwrap_or(cx.len().min(cy.len()));
                                let lo = at.saturating_sub(80);
                                eprintln!(
                                    "  {} 0x{kind:02x} first difference:\n    original:  ...{}\n    rewritten: ...{}",
                                    path.file_name().unwrap().to_string_lossy(),
                                    &cx[lo..(at + 120).min(cx.len())],
                                    &cy[lo..(at + 120).min(cy.len())]
                                );
                            }
                        }
                        _ => *xml_diffs.entry(format!("0x{kind:02x} unparsable")).or_default() += 1,
                    }
                }
            }
        }
    }

    eprintln!("corpus: {} files found, {read_ok} read", files.len());
    eprintln!("SnapGene files rewritten byte for byte: {byte_identical}/{snapgene_files}");
    eprintln!("primer binding sites: {primer_sites_exact}/{primer_sites} match the template exactly at the stored coordinates");
    for example in &primer_mismatch_examples {
        eprintln!("  {example}");
    }
    for (what, count) in &xml_diffs {
        eprintln!("  xml: {what}: {count}");
    }
    for failure in &failures {
        eprintln!("FAIL {failure}");
    }
    assert!(failures.is_empty(), "{} corpus failures", failures.len());
}
