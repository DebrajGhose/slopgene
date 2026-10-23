//! Stretches of sequence that may run through the origin of a circular molecule.

use crate::error::{Result, invalid};

/// A stretch of sequence, 0-based and half-open: it covers bases `start..end`.
///
/// On a circular sequence of length `L` a span may run through the origin. Such a span keeps
/// counting past the end, so `end` lies in `L+1..=start+L`. For example, on a 100 bp plasmid the
/// span 90..110 covers bases 90..100 and then 0..10. A span never covers more than `L` bases, and
/// `start` is always below `L` once the span is normalised.
///
/// User-facing coordinates are 1-based and inclusive (like SnapGene and GenBank); convert with
/// [`Span::from_1based`] and [`Span::to_1based`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub struct Span {
    pub start: usize,
    pub end: usize,
}

impl Span {
    pub const fn new(start: usize, end: usize) -> Span {
        Span { start, end }
    }

    /// Number of bases covered.
    pub fn len(&self) -> usize {
        self.end.saturating_sub(self.start)
    }

    pub fn is_empty(&self) -> bool {
        self.end <= self.start
    }

    /// Whether the span runs through the origin of a sequence of length `seq_len`.
    pub fn wraps(&self, seq_len: usize) -> bool {
        self.end > seq_len
    }

    /// Bring `start` into `0..seq_len`, keeping the length.
    pub fn normalize(self, seq_len: usize) -> Span {
        if seq_len == 0 {
            return self;
        }
        let start = self.start % seq_len;
        Span::new(start, start + self.len().min(seq_len))
    }

    /// The one or two linear pieces covered, each as `start..end` inside `0..seq_len`.
    pub fn pieces(&self, seq_len: usize) -> Vec<(usize, usize)> {
        if self.end <= seq_len { vec![(self.start, self.end)] } else { vec![(self.start, seq_len), (0, self.end - seq_len)] }
    }

    /// Whether base `pos` (in `0..seq_len`) lies inside the span.
    pub fn contains(&self, pos: usize, seq_len: usize) -> bool {
        if pos >= self.start && pos < self.end {
            return true;
        }
        // A wrapping span also covers 0..(end - L) at the start of the sequence.
        self.end > seq_len && pos + seq_len < self.end
    }

    /// Whether every base of `other` lies inside this span.
    pub fn contains_span(&self, other: &Span, seq_len: usize) -> bool {
        if other.is_empty() {
            return self.contains(other.start % seq_len.max(1), seq_len) || other.start == self.end;
        }
        if other.len() > self.len() {
            return false;
        }
        // Shift `other` so it starts inside this span's unwrapped coordinates, then compare ends.
        let mut start = other.start;
        if start < self.start && seq_len > 0 {
            start += seq_len;
        }
        start >= self.start && start + other.len() <= self.end
    }

    /// Whether the two spans share at least one base.
    pub fn overlaps(&self, other: &Span, seq_len: usize) -> bool {
        for (a0, a1) in self.pieces(seq_len) {
            for (b0, b1) in other.pieces(seq_len) {
                if a0 < b1 && b0 < a1 {
                    return true;
                }
            }
        }
        false
    }

    /// Copy the covered bases out of `seq`, following the span through the origin if needed.
    pub fn extract(&self, seq: &[u8]) -> Vec<u8> {
        let seq_len = seq.len();
        if self.end <= seq_len {
            seq[self.start..self.end].to_vec()
        } else {
            let mut out = seq[self.start..].to_vec();
            out.extend_from_slice(&seq[..self.end - seq_len]);
            out
        }
    }

    /// 1-based inclusive `(first, last)`. For a span through the origin `first > last`.
    pub fn to_1based(&self, seq_len: usize) -> (usize, usize) {
        let first = self.start + 1;
        let last = if self.end > seq_len { self.end - seq_len } else { self.end };
        (first, last)
    }

    /// Build a span from 1-based inclusive coordinates.
    ///
    /// On a circular sequence `first > last` means the span runs through the origin. On a
    /// linear sequence that is an error. `first == last + 1` on a circular sequence describes the
    /// whole molecule starting at `first`.
    pub fn from_1based(first: usize, last: usize, seq_len: usize, circular: bool) -> Result<Span> {
        if first == 0 || last == 0 {
            return Err(invalid("coordinates are 1-based; 0 is not a position"));
        }
        if first > seq_len || last > seq_len {
            return Err(invalid(format!("{first}..{last} lies outside the sequence (length {seq_len})")));
        }
        if first <= last {
            Ok(Span::new(first - 1, last))
        } else if circular {
            Ok(Span::new(first - 1, last + seq_len))
        } else {
            Err(invalid(format!("start {first} is after end {last} on a linear sequence")))
        }
    }

    /// Format as `first..last` (1-based inclusive), the convention used across SlopGene.
    pub fn format(&self, seq_len: usize) -> String {
        let (first, last) = self.to_1based(seq_len);
        format!("{first}..{last}")
    }
}

/// Parse a user-facing range such as `120..450`, `120-450` or `120:450` (1-based inclusive).
/// A single number is a one-base range.
pub fn parse_range(text: &str, seq_len: usize, circular: bool) -> Result<Span> {
    let text = text.trim();
    let parts: Vec<&str> = if let Some((a, b)) = text.split_once("..") {
        vec![a, b]
    } else if let Some((a, b)) = text.split_once(':') {
        vec![a, b]
    } else if let Some((a, b)) = text.split_once('-') {
        vec![a, b]
    } else {
        vec![text]
    };
    let parse = |s: &str| s.trim().replace([',', '_'], "").parse::<usize>().map_err(|_| invalid(format!("'{s}' is not a position")));
    let first = parse(parts[0])?;
    let last = if parts.len() == 2 { parse(parts[1])? } else { first };
    Span::from_1based(first, last, seq_len, circular)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wrapping_span_basics() {
        let span = Span::from_1based(91, 10, 100, true).unwrap();
        assert_eq!(span, Span::new(90, 110));
        assert_eq!(span.len(), 20);
        assert!(span.wraps(100));
        assert_eq!(span.pieces(100), vec![(90, 100), (0, 10)]);
        assert!(span.contains(95, 100));
        assert!(span.contains(5, 100));
        assert!(!span.contains(50, 100));
        assert_eq!(span.to_1based(100), (91, 10));
        assert_eq!(span.format(100), "91..10");
    }

    #[test]
    fn extract_through_origin() {
        let seq = b"ABCDEFGHIJ";
        assert_eq!(Span::new(8, 12).extract(seq), b"IJAB".to_vec());
        assert_eq!(Span::new(2, 5).extract(seq), b"CDE".to_vec());
    }

    #[test]
    fn linear_rejects_backwards_range() {
        assert!(Span::from_1based(50, 10, 100, false).is_err());
        assert!(Span::from_1based(0, 10, 100, false).is_err());
        assert!(Span::from_1based(5, 101, 100, false).is_err());
    }

    #[test]
    fn overlaps_and_containment() {
        let a = Span::new(90, 110);
        assert!(a.overlaps(&Span::new(5, 15), 100));
        assert!(!a.overlaps(&Span::new(20, 30), 100));
        assert!(a.contains_span(&Span::new(2, 8), 100));
        assert!(a.contains_span(&Span::new(95, 105), 100));
        assert!(!a.contains_span(&Span::new(5, 15), 100));
    }

    #[test]
    fn range_parsing() {
        assert_eq!(parse_range("10..20", 100, false).unwrap(), Span::new(9, 20));
        assert_eq!(parse_range("10-20", 100, false).unwrap(), Span::new(9, 20));
        assert_eq!(parse_range("1,000..1,010", 2000, false).unwrap(), Span::new(999, 1010));
        assert_eq!(parse_range("7", 100, false).unwrap(), Span::new(6, 7));
        assert_eq!(parse_range("95..5", 100, true).unwrap(), Span::new(94, 105));
    }
}
