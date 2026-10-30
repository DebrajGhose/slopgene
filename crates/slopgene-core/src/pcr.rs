//! PCR simulation.
//!
//! The product runs from the forward primer's 5' end to the reverse primer's 5' end. Primer
//! bases are copied from the primers, not the template, so 5' tails (restriction sites,
//! assembly overlaps) and deliberate mismatches end up in the product, as they do at the bench.

use crate::alphabet;
use crate::edit;
use crate::error::{Error, Result};
use crate::model::*;
use crate::primers::{Binding, BindingRules, find_binding_sites};
use crate::span::Span;

/// A simulated PCR product.
#[derive(Clone, Debug)]
pub struct PcrProduct {
    pub doc: Document,
    /// Template region between the two annealed 3' ends, inclusive of the annealed bases.
    pub template_span: Span,
    pub forward_site: Binding,
    pub reverse_site: Binding,
    /// Other primer-pair combinations that would also give a product (possible side products).
    pub alternatives: usize,
}

/// Run a PCR on `template` with a forward and a reverse primer (both written 5'->3').
pub fn pcr(template: &Document, forward: &Primer, reverse: &Primer, rules: &BindingRules) -> Result<PcrProduct> {
    if !template.molecule.is_nucleic() {
        return Err(Error::Invalid("PCR needs a DNA template".into()));
    }
    let seq_len = template.len();
    let circular = template.is_circular();

    // A primer can be listed either way round, so collect sites of both primers on both
    // strands and pair every top-strand site with every bottom-strand site.
    let mut tops: Vec<(Binding, &Primer)> = Vec::new();
    let mut bottoms: Vec<(Binding, &Primer)> = Vec::new();
    for primer in [forward, reverse] {
        for site in find_binding_sites(&template.seq, circular, primer.seq.as_bytes(), rules) {
            match site.strand {
                Strand::Forward => tops.push((site, primer)),
                _ => bottoms.push((site, primer)),
            }
        }
    }
    if tops.is_empty() || bottoms.is_empty() {
        return Err(Error::Invalid(format!(
            "no product: {} binds {} top-strand and {} bottom-strand sites in total",
            if tops.is_empty() && bottoms.is_empty() { "neither primer" } else { "the pair" },
            tops.len(),
            bottoms.len()
        )));
    }

    // Each pairing gives a product spanning from the top site's start to the bottom site's end,
    // read along the top strand (through the origin on circular templates).
    // A primer with one place it binds.
    type Bound<'a> = (Binding, &'a Primer);
    let mut products: Vec<(usize, Span, &Bound, &Bound)> = Vec::new();
    for top in &tops {
        for bottom in &bottoms {
            let start = top.0.span.start;
            let mut end = bottom.0.span.end;
            if circular {
                if end <= start {
                    end += seq_len;
                }
            } else if end <= start {
                continue;
            }
            let span = Span::new(start, end);
            if span.len() > if circular { seq_len } else { seq_len + 1 } || span.len() < top.0.span.len().max(bottom.0.span.len()) {
                continue;
            }
            products.push((span.len(), span, top, bottom));
        }
    }
    products.sort_by_key(|p| p.0);
    let Some(&(_, span, top, bottom)) = products.first() else {
        return Err(Error::Invalid("the primers bind but do not face each other, so there is no product".into()));
    };
    let alternatives = products.len() - 1;

    // Assemble: whole forward primer + template between the annealed parts + reverse complement
    // of the whole reverse primer.
    let (top_site, top_primer) = (&top.0, top.1);
    let (bottom_site, bottom_primer) = (&bottom.0, bottom.1);
    let top_full = alphabet::upper(top_primer.seq.as_bytes());
    let bottom_full_rc = alphabet::revcomp(&alphabet::upper(bottom_primer.seq.as_bytes()));
    let top_tail = top_full.len() - top_site.annealed;
    let inner_start = top_site.span.end;
    let inner_end = if circular && bottom_site.span.start < inner_start { bottom_site.span.start + seq_len } else { bottom_site.span.start };

    let mut product_seq = top_full.clone();
    if inner_end >= inner_start {
        let inner = Span::new(inner_start, inner_end);
        if !inner.is_empty() {
            product_seq.extend(inner.extract(&template.seq));
        }
        product_seq.extend_from_slice(&bottom_full_rc);
    } else {
        // The annealed parts overlap (a very short product): share the overlapping bases.
        let overlap = inner_start - inner_end;
        if overlap > bottom_full_rc.len() {
            return Err(Error::Invalid("the primers overlap more than their own length".into()));
        }
        product_seq.extend_from_slice(&bottom_full_rc[overlap..]);
    }

    // Features from the template region, shifted past the forward primer's 5' tail.
    let mut product =
        Document::new_dna(format!("{} PCR {}-{}", template.name, top_primer.name, bottom_primer.name).trim().to_string(), product_seq, Topology::Linear);
    // A polymerase makes unmethylated DNA, so a site that methylation blocks on the template is
    // cut on the product.
    product.methylation = Methylation::NONE;
    if let Ok(region) = edit::extract(template, span) {
        for mut feature in region.features {
            for s in &mut feature.segments {
                s.span = Span::new(s.span.start + top_tail, s.span.end + top_tail);
            }
            let fits = feature.segments.iter().all(|s| s.span.end <= product.len());
            if fits {
                product.features.push(feature);
            }
        }
    }
    let product_len = product.len();
    let mut add_primer = |primer: &Primer, span: Span, strand: Strand| {
        let mut copy = Primer::new(primer.name.clone(), primer.seq.clone());
        copy.description = primer.description.clone();
        copy.sites.push(BindingSite { span, strand, tm: None, extra: XmlExtra::default() });
        product.primers.push(copy);
    };
    add_primer(top_primer, Span::new(0, top_full.len()), Strand::Forward);
    add_primer(bottom_primer, Span::new(product_len - bottom_full_rc.len(), product_len), Strand::Reverse);
    product.meta.description = format!("PCR product of {} with {} and {} ({} bp)", template.name, top_primer.name, bottom_primer.name, product_len);
    product.meta.synthetic = Some(true);
    Ok(PcrProduct { doc: product, template_span: span, forward_site: top_site.clone(), reverse_site: bottom_site.clone(), alternatives })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn product_carries_primer_tails() {
        // Template: 20 T, a 60 bp region, 20 T.
        let region = b"ATGGCTAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGAC";
        let mut seq = vec![b'T'; 20];
        seq.extend_from_slice(region);
        seq.extend(vec![b'T'; 20]);
        let mut template = Document::new_dna("tmpl", seq, Topology::Linear);
        template.features.push(Feature::new("gfp-part", "CDS", Span::new(20, 80), Strand::Forward));

        let fwd = Primer::new("F", "GGATCCATGGCTAGCAAGGGCGAGGAG"); // BamHI tail
        let rev_rc = &region[region.len() - 21..];
        let rev = Primer::new("R", format!("GAATTC{}", String::from_utf8(alphabet::revcomp(rev_rc)).unwrap()));
        let product = pcr(&template, &fwd, &rev, &BindingRules::default()).unwrap();
        // The reverse tail appears as its reverse complement at the 3' end of the top strand
        // (GAATTC is its own reverse complement, so it reads the same).
        let mut expected = b"GGATCC".to_vec();
        expected.extend_from_slice(region);
        expected.extend(alphabet::revcomp(b"GAATTC"));
        assert_eq!(product.doc.seq, expected);
        assert_eq!(product.doc.features.len(), 1);
        assert_eq!(product.doc.features[0].segments[0].span, Span::new(6, 66));
        assert_eq!(product.doc.primers.len(), 2);
    }

    #[test]
    fn product_across_origin() {
        let region = b"GCTAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTC";
        // Put the region across the origin of a circular template.
        let mut seq = region[24..].to_vec();
        seq.extend(vec![b'A'; 40]);
        seq.extend_from_slice(&region[..24]);
        let template = Document::new_dna("circ", seq, Topology::Circular);
        let fwd = Primer::new("F", String::from_utf8(region[..20].to_vec()).unwrap());
        let rev = Primer::new("R", String::from_utf8(alphabet::revcomp(&region[region.len() - 20..])).unwrap());
        let product = pcr(&template, &fwd, &rev, &BindingRules::default()).unwrap();
        assert_eq!(product.doc.seq, region.to_vec());
    }
}
