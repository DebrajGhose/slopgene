//! A small XML DOM, just enough for the XML that SnapGene embeds in its files.
//!
//! SnapGene stores features, primers, notes and settings as short XML documents. Round-tripping
//! a file without losing anything means keeping every element and attribute we do not model, in
//! order. A tiny DOM that keeps attribute order and text nodes does that.

use crate::error::{Result, parse_err};

/// One XML element: its name, attributes in document order, and children.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct Element {
    pub name: String,
    pub attrs: Vec<(String, String)>,
    pub children: Vec<Node>,
}

/// A child of an element.
#[derive(Clone, Debug, PartialEq)]
pub enum Node {
    Element(Element),
    Text(String),
    CData(String),
    Comment(String),
}

impl Element {
    pub fn new(name: impl Into<String>) -> Element {
        Element { name: name.into(), attrs: Vec::new(), children: Vec::new() }
    }

    /// Value of attribute `name`, if present.
    pub fn attr(&self, name: &str) -> Option<&str> {
        self.attrs.iter().find(|(k, _)| k == name).map(|(_, v)| v.as_str())
    }

    /// Set attribute `name`, replacing an existing value in place or appending a new one.
    pub fn set_attr(&mut self, name: &str, value: impl Into<String>) {
        let value = value.into();
        match self.attrs.iter_mut().find(|(k, _)| k == name) {
            Some(slot) => slot.1 = value,
            None => self.attrs.push((name.to_string(), value)),
        }
    }

    /// Builder form of [`Element::set_attr`].
    pub fn with_attr(mut self, name: &str, value: impl Into<String>) -> Element {
        self.set_attr(name, value);
        self
    }

    pub fn remove_attr(&mut self, name: &str) -> Option<String> {
        let index = self.attrs.iter().position(|(k, _)| k == name)?;
        Some(self.attrs.remove(index).1)
    }

    /// Child elements, skipping text and comments.
    pub fn elements(&self) -> impl Iterator<Item = &Element> {
        self.children.iter().filter_map(|node| match node {
            Node::Element(e) => Some(e),
            _ => None,
        })
    }

    /// Child elements with the given name.
    pub fn children_named<'a>(&'a self, name: &'a str) -> impl Iterator<Item = &'a Element> + 'a {
        self.elements().filter(move |e| e.name == name)
    }

    /// First child element with the given name.
    pub fn child(&self, name: &str) -> Option<&Element> {
        self.elements().find(|e| e.name == name)
    }

    pub fn child_mut(&mut self, name: &str) -> Option<&mut Element> {
        self.children.iter_mut().find_map(|node| match node {
            Node::Element(e) if e.name == name => Some(e),
            _ => None,
        })
    }

    /// Text of the element: its direct text and CDATA children joined together.
    pub fn text(&self) -> String {
        let mut out = String::new();
        for node in &self.children {
            match node {
                Node::Text(t) | Node::CData(t) => out.push_str(t),
                _ => {}
            }
        }
        out
    }

    /// Text of the first child element named `name`.
    pub fn child_text(&self, name: &str) -> Option<String> {
        self.child(name).map(|e| e.text())
    }

    /// Replace all children with a single text node.
    pub fn set_text(&mut self, text: impl Into<String>) {
        self.children = vec![Node::Text(text.into())];
    }

    pub fn push(&mut self, child: Element) {
        self.children.push(Node::Element(child));
    }

    /// Remove every child element named `name`.
    pub fn remove_children_named(&mut self, name: &str) {
        self.children.retain(|node| !matches!(node, Node::Element(e) if e.name == name));
    }

    /// Set the text of child `name`, creating the child (at the end) when it is missing.
    pub fn set_child_text(&mut self, name: &str, text: impl Into<String>) {
        let text = text.into();
        match self.child_mut(name) {
            Some(child) => child.set_text(text),
            None => {
                let mut child = Element::new(name);
                child.set_text(text);
                self.push(child);
            }
        }
    }

    /// Serialise the element (no XML declaration).
    pub fn to_xml(&self) -> String {
        let mut out = String::new();
        self.write_to(&mut out);
        out
    }

    pub fn write_to(&self, out: &mut String) {
        out.push('<');
        out.push_str(&self.name);
        for (key, value) in &self.attrs {
            out.push(' ');
            out.push_str(key);
            out.push_str("=\"");
            out.push_str(&escape_attr(value));
            out.push('"');
        }
        if self.children.is_empty() {
            out.push_str("/>");
            return;
        }
        out.push('>');
        for node in &self.children {
            match node {
                Node::Element(e) => e.write_to(out),
                Node::Text(t) => out.push_str(&escape_text(t)),
                Node::CData(t) => {
                    out.push_str("<![CDATA[");
                    out.push_str(t);
                    out.push_str("]]>");
                }
                Node::Comment(t) => {
                    out.push_str("<!--");
                    out.push_str(t);
                    out.push_str("-->");
                }
            }
        }
        out.push_str("</");
        out.push_str(&self.name);
        out.push('>');
    }
}

/// Escape text for use inside a double-quoted attribute.
///
/// Follows SnapGene's style: `&`, `<` and `"` are escaped, `>` is left alone (it is legal
/// there), so files that were not edited come back byte for byte.
pub fn escape_attr(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for ch in text.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '"' => out.push_str("&quot;"),
            '\n' => out.push_str("&#10;"),
            '\r' => out.push_str("&#13;"),
            '\t' => out.push_str("&#9;"),
            other => out.push(other),
        }
    }
    out
}

/// Escape text for use as element content (`&` and `<`; `>` only where it would close a
/// CDATA-like `]]>`).
pub fn escape_text(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for ch in text.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' if out.ends_with("]]") => out.push_str("&gt;"),
            other => out.push(other),
        }
    }
    out
}

/// A canonical text form of an element for comparing meaning rather than bytes: attributes
/// sorted, comments and whitespace-only text dropped, text trimmed.
pub fn canonical(e: &Element) -> String {
    let mut out = String::new();
    canonical_into(e, &mut out);
    out
}

fn canonical_into(e: &Element, out: &mut String) {
    let mut attrs = e.attrs.clone();
    attrs.sort();
    out.push('<');
    out.push_str(&e.name);
    for (k, v) in attrs {
        out.push(' ');
        out.push_str(&k);
        out.push_str("=\"");
        out.push_str(&escape_attr(&v));
        out.push('"');
    }
    out.push('>');
    for child in &e.children {
        match child {
            Node::Element(c) => canonical_into(c, out),
            Node::Text(t) | Node::CData(t) if !t.trim().is_empty() => out.push_str(&escape_text(t.trim())),
            _ => {}
        }
    }
    out.push_str("</");
    out.push_str(&e.name);
    out.push('>');
}

/// Replace XML character and entity references with the characters they stand for.
/// Unknown entities are kept as written so nothing is lost.
pub fn unescape(text: &str) -> String {
    if !text.contains('&') {
        return text.to_string();
    }
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(amp) = rest.find('&') {
        out.push_str(&rest[..amp]);
        let after = &rest[amp + 1..];
        // Entities are short; give up on anything without a ';' in the next 12 bytes.
        let Some(semi) = after.find(';').filter(|&s| s <= 12) else {
            out.push('&');
            rest = after;
            continue;
        };
        let entity = &after[..semi];
        let replacement = match entity {
            "lt" => Some('<'),
            "gt" => Some('>'),
            "amp" => Some('&'),
            "quot" => Some('"'),
            "apos" => Some('\''),
            _ if entity.starts_with("#x") || entity.starts_with("#X") => u32::from_str_radix(&entity[2..], 16).ok().and_then(char::from_u32),
            _ if entity.starts_with('#') => entity[1..].parse::<u32>().ok().and_then(char::from_u32),
            _ => None,
        };
        match replacement {
            Some(ch) => out.push(ch),
            None => {
                out.push('&');
                out.push_str(entity);
                out.push(';');
            }
        }
        rest = &after[semi + 1..];
    }
    out.push_str(rest);
    out
}

/// Parse an XML document and return its root element.
///
/// The prolog (XML declaration, comments, doctype) is skipped. Text after the root is ignored.
pub fn parse(text: &str) -> Result<Element> {
    let mut parser = Parser { src: text.as_bytes(), pos: 0, text };
    parser.skip_bom();
    loop {
        parser.skip_whitespace();
        if parser.starts_with("<?") {
            parser.skip_past("?>")?;
        } else if parser.starts_with("<!--") {
            parser.skip_past("-->")?;
        } else if parser.starts_with("<!") {
            parser.skip_past(">")?;
        } else if parser.starts_with("<") {
            return parser.parse_element();
        } else if parser.pos >= parser.src.len() {
            return Err(parse_err("XML document has no root element"));
        } else {
            return Err(parse_err(format!("unexpected text before XML root at byte {}", parser.pos)));
        }
    }
}

struct Parser<'a> {
    src: &'a [u8],
    pos: usize,
    text: &'a str,
}

impl<'a> Parser<'a> {
    fn skip_bom(&mut self) {
        if self.src.starts_with(&[0xEF, 0xBB, 0xBF]) {
            self.pos = 3;
        }
    }

    fn starts_with(&self, prefix: &str) -> bool {
        self.src[self.pos..].starts_with(prefix.as_bytes())
    }

    fn skip_whitespace(&mut self) {
        while self.pos < self.src.len() && self.src[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
    }

    /// Move past the next occurrence of `marker`.
    fn skip_past(&mut self, marker: &str) -> Result<()> {
        match self.text[self.pos..].find(marker) {
            Some(offset) => {
                self.pos += offset + marker.len();
                Ok(())
            }
            None => Err(parse_err(format!("unterminated XML construct, expected '{marker}'"))),
        }
    }

    /// Read an element or attribute name.
    fn parse_name(&mut self) -> Result<String> {
        let start = self.pos;
        while self.pos < self.src.len() {
            let b = self.src[self.pos];
            if b.is_ascii_whitespace() || matches!(b, b'/' | b'>' | b'=' | b'<') {
                break;
            }
            self.pos += 1;
        }
        if self.pos == start {
            return Err(parse_err(format!("expected an XML name at byte {start}")));
        }
        Ok(self.text[start..self.pos].to_string())
    }

    fn parse_element(&mut self) -> Result<Element> {
        // We are at '<'.
        self.pos += 1;
        let name = self.parse_name()?;
        let mut element = Element::new(name);

        // Attributes, up to '>' or '/>'.
        loop {
            self.skip_whitespace();
            if self.pos >= self.src.len() {
                return Err(parse_err(format!("unterminated start tag <{}>", element.name)));
            }
            if self.starts_with("/>") {
                self.pos += 2;
                return Ok(element);
            }
            if self.src[self.pos] == b'>' {
                self.pos += 1;
                break;
            }
            let key = self.parse_name()?;
            self.skip_whitespace();
            if self.pos >= self.src.len() || self.src[self.pos] != b'=' {
                return Err(parse_err(format!("attribute '{key}' has no value")));
            }
            self.pos += 1;
            self.skip_whitespace();
            let quote = *self.src.get(self.pos).ok_or_else(|| parse_err("unterminated attribute"))?;
            if quote != b'"' && quote != b'\'' {
                return Err(parse_err(format!("attribute '{key}' value is not quoted")));
            }
            self.pos += 1;
            let start = self.pos;
            while self.pos < self.src.len() && self.src[self.pos] != quote {
                self.pos += 1;
            }
            if self.pos >= self.src.len() {
                return Err(parse_err(format!("unterminated value for attribute '{key}'")));
            }
            let value = unescape(&self.text[start..self.pos]);
            self.pos += 1;
            element.attrs.push((key, value));
        }

        // Content, up to the matching end tag.
        loop {
            if self.pos >= self.src.len() {
                return Err(parse_err(format!("element <{}> is never closed", element.name)));
            }
            if self.starts_with("</") {
                self.pos += 2;
                let end_name = self.parse_name()?;
                if end_name != element.name {
                    return Err(parse_err(format!("end tag </{end_name}> does not match <{}>", element.name)));
                }
                self.skip_whitespace();
                if self.src.get(self.pos) != Some(&b'>') {
                    return Err(parse_err(format!("malformed end tag </{end_name}>")));
                }
                self.pos += 1;
                return Ok(element);
            } else if self.starts_with("<![CDATA[") {
                let start = self.pos + 9;
                self.skip_past("]]>")?;
                element.children.push(Node::CData(self.text[start..self.pos - 3].to_string()));
            } else if self.starts_with("<!--") {
                let start = self.pos + 4;
                self.skip_past("-->")?;
                element.children.push(Node::Comment(self.text[start..self.pos - 3].to_string()));
            } else if self.starts_with("<?") {
                self.skip_past("?>")?;
            } else if self.src[self.pos] == b'<' {
                let child = self.parse_element()?;
                element.children.push(Node::Element(child));
            } else {
                let start = self.pos;
                while self.pos < self.src.len() && self.src[self.pos] != b'<' {
                    self.pos += 1;
                }
                element.children.push(Node::Text(unescape(&self.text[start..self.pos])));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_snapgene_style_features() {
        let xml = r##"<?xml version="1.0"?><Features nextValidID="2"><Feature name="lac &lt;i&gt;op" type="protein_bind"><Segment range="489-505" color="#31849b"/><Q name="note"><V text="&lt;html>&lt;body>hi&lt;/body>&lt;/html>"/></Q></Feature></Features>"##;
        let root = parse(xml).unwrap();
        assert_eq!(root.name, "Features");
        assert_eq!(root.attr("nextValidID"), Some("2"));
        let feature = root.child("Feature").unwrap();
        assert_eq!(feature.attr("name"), Some("lac <i>op"));
        let v = feature.child("Q").unwrap().child("V").unwrap();
        assert_eq!(v.attr("text"), Some("<html><body>hi</body></html>"));
    }

    #[test]
    fn round_trip_keeps_order_and_text() {
        let xml = "<Notes>\n<UUID>abc</UUID>\n<Created UTC=\"0:8:26\">2023.5.19</Created>\n</Notes>";
        let root = parse(xml).unwrap();
        assert_eq!(root.child_text("Created").as_deref(), Some("2023.5.19"));
        assert_eq!(root.child("Created").unwrap().attr("UTC"), Some("0:8:26"));
        let again = parse(&root.to_xml()).unwrap();
        assert_eq!(again, root);
    }

    #[test]
    fn numeric_entities_and_unknown_entities() {
        assert_eq!(unescape("a&#65;&#x42;c"), "aABc");
        assert_eq!(unescape("x &nbsp; y"), "x &nbsp; y");
        assert_eq!(unescape("AT&T"), "AT&T");
    }

    #[test]
    fn mismatched_tags_are_errors() {
        assert!(parse("<a><b></a>").is_err());
        assert!(parse("<a>").is_err());
    }
}
