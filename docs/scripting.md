# Scripting SlopGene

SlopGene embeds the [Rhai](https://rhai.rs) scripting language for work that one command cannot do alone.
A script opens files, inspects and edits documents, runs analyses and saves the results.
The same engine runs from the command line and from the console of the desktop app.
`slopgene run --api` lists every function described here.

## Running a script

The `run` command takes a script file, or a line of code after `-e`.

```
slopgene run examples/scripts/library_summary.rhai --arg dir=plasmids --arg recursive=true
slopgene run -e 'open("pUC19.dna").length'
slopgene run -e 'open("pUC19.dna").cut_counts("EcoRI,BamHI")' --json
```

Each `--arg key=value` adds an entry to the `ARGS` map.
Its values are always strings.
A script tests for an argument with `"key" in ARGS`.

In text mode, `print` writes each line as the script runs.
The value of the script's last expression is printed at the end, unless it is empty.
With `--json`, the command prints one JSON object instead.
That object has the fields `ok`, `output` (everything the script printed) and `result` (the last value).
A document in `result` appears as SlopGene document JSON without its sequence.

A failing script exits with status 1.
Its message names the problem and the line.
With `--json` the message arrives as `{"ok": false, "error": "..."}` on standard output.

### In the desktop app

Ctrl+J opens the script console.
There, the variable `DOC` holds the document in the active tab.
Changes a script makes to `DOC` are applied to the open document when the script ends.
Ctrl+Z undoes them in one step.

## The language in brief

Rhai reads like a small mix of JavaScript and Rust.
This section covers what SlopGene scripts use most.
The [Rhai book](https://rhai.rs/book/) describes the whole language.

```rhai
let doc = open("pUC19.dna");             // variables
const TAB = "\t";                        // constants
let names = [];                          // arrays
let info = #{ name: doc.name, bp: doc.length };  // maps (object literals)

if doc.is_circular && doc.length > 2000 {
    print(`${doc.name} is a plasmid of ${doc.length} bp`);
}

for f in doc.features {                  // loop over an array
    if f.type == "CDS" { names.push(f.name); }
}
print(join(names, ", "));

fn inside(pos, start, end) {             // functions
    if start <= end { pos >= start && pos <= end } else { pos >= start || pos <= end }
}
```

A few details matter in practice.

- Backtick strings interpolate `${...}`, but they keep backslashes as written. Write a tab or a newline in a double-quoted string, as in `const TAB = "\t";`, and interpolate that constant.
- `doc.subseq(1, 10)` and `subseq(doc, 1, 10)` are the same call. Both styles change the variable in place for editing functions.
- Assigning a document to another variable copies it. Editing the copy leaves the original unchanged.
- The variable of a `for` loop is a copy of each element. To edit documents held in an array, index the array, as in `docs[0].insert(1, "GG")`.
- Integers are 64-bit (`int`) and decimals are `float`. `x.to_int()` and `x.to_float()` convert between them.
- `try { ... } catch (e) { ... }` catches an error, such as a file that cannot be opened.

## Conventions

SlopGene's scripts follow the same conventions as its command line.

- **Positions are 1-based and inclusive**, as in SnapGene and GenBank. `doc.subseq(100, 250)` returns 151 bases.
- **Ranges can run through the origin** of a circular molecule. On a 3,000 bp plasmid, `doc.subseq(2950, 50)` returns 101 bases.
- **Strands** are the strings `"forward"`, `"reverse"`, `"both"` and `"none"`.
- **Enzymes** are given as a comma-separated list (`"EcoRI,BamHI"`) or as a set name. The sets are `common`, `6+`, `8+`, `neb`, `thermo`, `golden-gate`, `type-iis`, `blunt` and `all`.
- **Features** are found by exact name in any case, or by `"#N"` for the Nth feature. A name that several features share is an error that suggests `"#N"`.
- **Cut positions** are top-strand positions, as SnapGene labels them. A cut at position 396 falls after base 396.
- **Primer binding sites** are recomputed after every edit that changes the sequence.
- **Paths** are relative to the current folder. Forward slashes work on every platform.

## Functions

### Files and folders

| Function | What it does |
|---|---|
| `glob(pattern) -> [string]` | Paths that match a pattern such as `"maps/**/*.dna"`, sorted. `*` matches within one folder and `**` crosses folders. |
| `files(dir) -> [string]` | Sequence files in a folder, sorted. |
| `files(dir, recursive) -> [string]` | Sequence files in a folder and, when recursive is true, in its subfolders. |
| `exists(path) -> bool` | Whether a file or folder exists. |
| `file_name(path) -> string` | The last part of a path, such as `"pUC19.dna"`. |
| `file_stem(path) -> string` | The file name without its extension, such as `"pUC19"`. |
| `read_text(path) -> string` | The contents of a text file. |
| `write_text(path, text)` | Write text to a file, replacing the file in one step. |
| `mkdir(path)` | Create a folder and any missing parent folders. |

`files` recognises the extensions `.dna`, `.rna`, `.prot`, `.gb`, `.gbk`, `.genbank`, `.gbff`, `.fa`, `.fasta`, `.fna`, `.faa`, `.fas`, `.json` and `.ape`.

### Opening and creating documents

| Function | What it does |
|---|---|
| `open(path) -> Doc` | Open a SnapGene, GenBank, FASTA or JSON file (the first record when it holds several). |
| `open_all(path) -> [Doc]` | Every record in a file. |
| `new_dna(name, seq) -> Doc` | A new linear DNA document. |
| `new_dna(name, seq, topology) -> Doc` | A new DNA document that is `"linear"` or `"circular"`. |
| `from_json(text) -> Doc` | A document from SlopGene document JSON, the format `to_json` writes. |
| `doc.copy() -> Doc` | A copy of the document. Assigning a document to a new variable also copies it. |

The format of a file is detected from its contents, so the extension does not need to match.

### Document properties

| Property | What it holds |
|---|---|
| `doc.name` | The name, which can be assigned. |
| `doc.description` | The description, which can be assigned. |
| `doc.length` | The number of bases or residues. |
| `doc.len() -> int` | The number of bases or residues, the same as `doc.length`. |
| `doc.seq` | The sequence as a string. |
| `doc.topology` | `"linear"` or `"circular"`. Assigning either value changes the topology. |
| `doc.molecule` | `"dna"`, `"rna"` or `"protein"`. |
| `doc.is_circular` | Whether the molecule is circular. |
| `doc.gc` | GC content in percent. |
| `doc.path` | The file the document was opened from or last saved to, or `""` for a new document. |
| `doc.features` | Every feature as a map (see below). |
| `doc.primers` | Every primer as a map with name, sequence, length, description and sites. |
| `doc.to_map() -> map` | The whole document as a map with the content of `to_json`. |
| `doc.to_json() -> string` | The whole document as SlopGene document JSON. |
| `doc.to_string() -> string` | The name, length and topology on one line, which is what `print(doc)` shows. |

A feature map has the keys `name`, `type`, `strand`, `start`, `end`, `length`, `segments`, `qualifiers` and `index`.
It has `color` as `"#rrggbb"` when the feature has a colour.
`index` is the feature's number for `"#N"` lookups.
`qualifiers` maps each key to an array of values, as in `f.qualifiers.gene[0]`.
Each primer site has `start`, `end`, `strand` and `tm`, the melting temperature of the bases that anneal there.

### Reading parts of a document

| Function | What it does |
|---|---|
| `doc.subseq(start, end) -> string` | The bases from start to end, through the origin when start > end on a circular molecule. |
| `doc.feature(name) -> map` | One feature, by name in any case or by `"#N"` for the Nth feature. |
| `doc.has_feature(name) -> bool` | Whether any feature has this name. |
| `doc.features_of_type(type) -> [map]` | The features of one type, such as `"CDS"`. |
| `doc.feature_seq(name) -> string` | A feature's bases, read 5' to 3' along its strand. |
| `doc.translate_feature(name) -> string` | The protein a CDS feature encodes. |
| `doc.extract(start, end) -> Doc` | A region as a new linear document, with the features that overlap it clipped to the region. |

### Editing in place

| Function | What it does |
|---|---|
| `doc.insert(pos, bases)` | Insert bases so that the first new base is at pos, from 1 to length + 1. |
| `doc.delete(start, end)` | Delete the bases from start to end. |
| `doc.replace(start, end, bases)` | Replace the bases from start to end with new bases. |
| `doc.invert(start, end)` | Reverse-complement a region in place, flipping the features wholly inside it. |
| `doc.reverse_complement()` | Reverse-complement the whole molecule, flipping every feature and primer with it. |
| `doc.set_origin(pos)` | Make pos the first base of a circular molecule. |
| `doc.add_feature(name, type, start, end, strand) -> map` | Add a feature whose strand is `"forward"`, `"reverse"`, `"both"` or `"none"`. |
| `doc.add_feature(spec) -> map` | Add a feature from a map with name, type, start, end, strand, color and qualifiers. |
| `doc.update_feature(name, changes)` | Change name, type, strand, color or range (`"100..250"`), and set any other key as a qualifier. |
| `doc.remove_feature(name) -> int` | Remove every feature with this name and return how many were removed. |
| `doc.add_primer(name, seq) -> map` | Add a primer and find where it binds. |
| `doc.remove_primer(name) -> int` | Remove every primer with this name and return how many were removed. |
| `doc.refresh_primers()` | Recompute where every primer binds. |

Features move with every edit, so they keep covering the same bases.
A replacement of the same length is a point mutation, and nothing moves.
Edits change only the document in memory until it is saved.

```rhai
let doc = open("construct.dna");
doc.add_feature(#{ name: "His tag", type: "CDS", start: 801, end: 818, strand: "forward",
                   color: "#515faa", qualifiers: #{ note: "6xHis" } });
doc.update_feature("His tag", #{ range: "801..821" });
doc.save();
```

### Saving and drawing

| Function | What it does |
|---|---|
| `doc.save() -> string` | Write the document back to the file it came from and return the path. |
| `doc.save(path) -> string` | Write to a file in the format its extension names (.dna, .gbk, .fasta or .json). |
| `doc.save(path, format) -> string` | Write in a named format: `"dna"`, `"genbank"`, `"fasta"` or `"json"`. |
| `doc.render(path) -> string` | Draw the map as SVG, or as PNG when the path ends in .png. |
| `doc.render(path, options) -> string` | Draw the map with options for width, height, enzymes, dark and colors. |

Writes are atomic, so another program never sees a half-written file.
An unedited SnapGene file is written back byte for byte.

The render options are all optional.

| Option | Meaning | Default |
|---|---|---|
| `width`, `height` | Size of the drawing in pixels | 900, 700 |
| `enzymes` | Enzymes whose unique sites are marked: a set, a list or `"none"` | `"common"` |
| `dark` | Draw for a dark background | `false` |
| `colors` | `"palette"` for SlopGene's colours, or `"file"` for the colours stored in the file | `"palette"` |

```rhai
doc.render("pUC19.png", #{ width: 1200, height: 900, enzymes: "EcoRI,BamHI,HindIII", dark: true });
```

### Restriction enzymes

| Function | What it does |
|---|---|
| `doc.sites(enzymes) -> [map]` | Every cutting site as a map with enzyme, start, end, strand and cut, in position order. |
| `doc.cut_counts(enzymes) -> map` | How many times each enzyme cuts, keyed by enzyme name. |
| `doc.unique_cutters(enzymes) -> [map]` | The enzymes that cut exactly once, as maps with enzyme and cut. |
| `doc.digest(enzymes) -> [map]` | The fragments, each with length, start, end, the enzyme and end type on each side, and the fragment itself as doc. |
| `enzyme(name) -> map` | One enzyme with its name, site, site_with_cuts and suppliers. |
| `enzymes(spec) -> [string]` | The enzyme names of a set or of a comma-separated list. |

A digest fragment has the keys `length`, `start`, `end`, `left_enzyme`, `right_enzyme`, `left_end`, `right_end` and `doc`.
The end types read `"blunt"`, `"5' overhang (4 nt)"` or `"3' overhang (4 nt)"`.
`left_enzyme` is empty (`()`) where a fragment ends at the end of a linear molecule.
`doc` is the fragment as a document with exact sticky ends, ready for `ligate`.

`site_with_cuts` shows the cut positions in the site, as in `G^AATTC` for EcoRI.
A site cut outside itself (a type IIS enzyme) shows the cut offsets, as in `GGTCTC(1/5)` for BsaI.

### Analysis

| Function | What it does |
|---|---|
| `doc.find(pattern) -> [map]` | Matches of a DNA motif on both strands, as maps with start, end, strand and mismatches. |
| `doc.find(pattern, mismatches) -> [map]` | Matches of a DNA motif with up to this many mismatched bases. |
| `doc.find_protein(motif) -> [map]` | Matches of a protein motif in all six reading frames. |
| `doc.orfs() -> [map]` | Open reading frames of at least 75 codons, as maps with start, end, strand, frame and codons. |
| `doc.orfs(min_codons) -> [map]` | Open reading frames of at least this many codons. |
| `align(reference, read) -> map` | The alignment of a read to a reference: verdict, identity, coverage, structural changes and each difference. |
| `doc.annotate(library) -> [map]` | Add the library features found at 95% identity or more, and return the added features. |
| `doc.annotate(library, min_identity) -> [map]` | Add the library features found at this identity in percent or more. |

DNA motifs may use IUPAC codes such as `N`, `R` and `Y`.

`align` finds the read in either orientation and through the origin of a circular reference.
It chains shared 15-mers across insertions and deletions of any size, as SnapGene and minimap2 do.
A read from a circular record, or nearly as long as a circular reference, is treated as a whole plasmid.
It returns the same map as `slopgene align --json`.

| Key | Meaning |
|---|---|
| `verdict` | `"exact match"`, `"structural change"`, `"not this plasmid"` and the other verdicts of the [agent guide](agent-guide.md) |
| `identity` | Matching columns as a percentage of aligned columns, outside structural changes |
| `orientation` | `"forward"`, or `"reverse"` when the read matches the bottom strand |
| `reference_start` | The reference position where the alignment begins |
| `reference_covered`, `coverage` | Reference bases the alignment accounts for, as a count and as a percentage |
| `read_length`, `read_aligned` | Read bases in all, and inside the alignment |
| `whole_read` | Whether the read was treated as a whole circular plasmid |
| `events` | One map per structural change of 50 bp or more (see below) |
| `differences` | One map per smaller difference (see below) |
| `method` | The alignment method and its settings: `name`, `seed_length`, `chaining`, `alignment`, `scores` and `structural_min_bp` |

Each event has `type` (`"deletion"`, `"insertion"`, `"duplication"`, `"inversion"` or `"replacement"`), `length`, `read_bases`, `features` and `description`.
It has `reference_start` and `reference_end`, or `after` for a change that adds bases at one point.

Each difference has `type` (`"substitution"`, `"insertion"` or `"deletion"`), `position`, `reference`, `read`, `features`, `effects` and `homopolymer`.
`position` is the first reference base involved.
For an insertion it is the reference base the extra bases follow.
`features` names the features the difference falls in.
`effects` describes the change to each protein it touches, such as `"AmpR: L191V"`, `"AmpR: L191 silent"` or `"lacZα: frameshift at codon 12"`.
`homopolymer` is true for an insertion or deletion of one repeated base inside a run of four or more of that base.
Such indels are the usual nanopore sequencing error.

`annotate` takes an array of library sources.
Each source is a document or the path of a file.
Every feature of 30 bp or more in the library becomes something to look for.
`files(dir, true)` passes a whole folder of trusted maps.

```rhai
let read = open("results/clone7.gbk");
let added = read.annotate(files("maps", true));
print(`${added.len()} features added`);
read.save("results/clone7_annotated.dna");
```

### Cloning

| Function | What it does |
|---|---|
| `doc.nearest_cuts(start, end, enzymes) -> map` | The enzymes that cut nearest the region on each side, nearest first, and those that cut inside it. |
| `doc.nearest_cuts(start, end, enzymes, options) -> map` | The same, with the options vector (a document) and count (enzymes on each side). |
| `doc.excision(start, end, left, right) -> map` | What cutting with `left` before the region and `right` after it gives: the piece that holds the region, its ends and its direction. |
| `doc.excision(start, end, left, right, options) -> map` | The same, with the option vector, on which the pair is checked too. |
| `doc.pcr(fwd, rev) -> Doc` | The PCR product of two primers, given as sequences or as names of primers in the document. |
| `ligate(left, right) -> Doc` | Join two linear fragments whose facing ends are compatible. |
| `circularize(fragment) -> Doc` | Close a linear fragment whose two ends are compatible. |
| `restriction_clone(vector, insert, enzymes) -> [Doc]` | Cut both molecules with the same enzymes and ligate, giving one product per orientation that fits. |
| `assemble(fragments) -> Doc` | Join fragments in order through end overlaps of at least 15 bp and close the circle. |
| `assemble(fragments, min_overlap) -> Doc` | Join fragments in order through end overlaps of at least this length. |
| `golden_gate(parts, enzyme) -> Doc` | Golden Gate assembly with a type IIS enzyme, the vector first and the parts in any order. |

PCR products carry the primers' own bases, so 5' tails appear in the product.
`assemble` models Gibson, NEBuilder HiFi and In-Fusion assembly.
Products keep the features of their parts.

```rhai
let vector = open("pUC57.dna");
let insert = open("mTagBFP2.dna");
let products = restriction_clone(vector, insert, "SalI,KpnI");
products[0].save("pUC57-mTagBFP2.dna");
```

`nearest_cuts` answers the question a map asks first: which enzymes cut close to this part?
It lists the enzymes that cut before the region and those that cut after it, nearest first.
Each enzyme appears once on each side, at its cut nearest the region.
The lists hold 8 enzymes each unless `count` says otherwise.
The result is a map with `before`, `after` and `inside`.
Each entry of `before` and `after` is a map with `enzyme`, `cut` (the bases before the cut), `distance` (the bases between the cut and the region), `end` (the end the cut leaves, such as `"5' AATT"`) and `cuts` (how often the enzyme cuts the map).
With a `vector`, each entry also has `vector_cuts`, the number of times the enzyme cuts the vector.
`inside` names the enzymes that cut inside the region, which would cut it too.
Nothing is paired or ranked, so the script chooses one enzyme on each side.

`excision` then says what a chosen pair cuts out.
The result is a map with `enzymes`, `left` and `right` (each with `enzyme`, `cut`, `extra` and `end`), `start`, `end`, `length`, `directional` and `other_cuts`.
`directional` is true when the two ends differ, so the insert goes in one way round.
`other_cuts` counts the pair's cuts elsewhere in the map.
With a `vector`, the result also carries a `vector` map with `cuts`, `dropped` (the bases the vector loses between the cuts) and `problem` (why the pair cannot open the vector, or `()`).
A pair that cannot cut the region out is an error that says why.

```rhai
let donor = open("donor.dna");
let gfp = donor.feature("GFP");
let near = donor.nearest_cuts(gfp.start, gfp.end, "common", #{ vector: open("pUC19.dna") });
for cut in near.before {
    print(`${cut.enzyme}\t${cut.distance} bp before\tcuts the vector ${cut.vector_cuts} time(s)`);
}
let piece = donor.excision(gfp.start, gfp.end, near.before[0].enzyme, near.after[0].enzyme);
print(`${piece.length} bp, ${piece.left.end} / ${piece.right.end}`);
```

### Cloning protocols

| Function | What it does |
|---|---|
| `protocol_run(path) -> map` | Run a protocol file and report every step, as `slopgene protocol run --json` does. |
| `protocol_new(path, vector, insert) -> map` | Write a protocol that cuts a piece out of the insert's map and puts it into the vector, then run it. |
| `protocol_new(path, vector, insert, options) -> map` | The same, with the options enzymes, vector_enzymes, cut_out, replace, name, keep_phosphates, pcr (copy the insert out by PCR), tm (the melting temperature its primer design aims for) and force. |
| `protocol_add(path, map) -> map` | Add a map to cut a piece from. Its cut waits for a part and enzymes. |
| `protocol_add(path, map, options) -> map` | The same, with the options cut_out (the part), enzymes and pcr (copy the part out by PCR). |
| `protocol_cut(path, step, options) -> map` | Change a cut's part, its enzymes or both, with the options part and enzymes (`""` clears either). For a PCR step, the part it copies and the enzymes on its primers' tails. |
| `protocol_keep(path, pick, piece) -> map` | Make a pick keep another piece of its cut, counted from 1 along the cut map. |
| `protocol_put(path, piece, into) -> map` | Put a kept piece into a map, which is cut with the piece's enzymes, or into a join. |
| `protocol_take_out(path, piece, join) -> map` | Take a piece out of a join. |
| `protocol_remove_map(path, step) -> map` | Take a map out of the protocol with every step that works on it. |
| `protocol_leave_out(path, step) -> map` | Leave out a step that has one input, and pass that input on. |
| `protocol_dephosphorylate(path, piece, join) -> map` | Take the 5' phosphates off a piece on its way into a join. |
| `protocol_turn(path, piece, join, way) -> map` | Say which way a piece goes into a join: "forward", "reverse" or "either". |
| `protocol_copy_out(path, step) -> map` | Copy a map's part out by PCR instead of cutting it out: a PCR step goes between the map and its cut, with the cut's enzymes on the primers' tails. |
| `protocol_cut_out(path, step) -> map` | Cut a map's part out instead of copying it by PCR, the reverse of `protocol_copy_out`. |
| `protocol_primers(path, step, options) -> map` | The primers of a PCR step: `pair` keeps pair N of those the run designs, `forward` and `reverse` give your own, `tm` aims the design at a melting temperature, and `clear: true` designs them anew. |

These functions make the changes the app's protocol bench makes with the mouse.
Each one changes the file, writes it back and runs it.
A change that cannot be made is an error that says why, and the file stays as it was.
Steps are named by their ids, such as `cut_donor`, `insert` and `product`.
The result is the run as `slopgene protocol run --json` reports it, with `path`.
After a change the result also has `changed`, a few words on what changed.
Maps are paths from the working folder, and the protocol stores them relative to its own folder.

`protocol_put` onto a map cuts that map with the piece's own enzymes, so the ends match.
It keeps the map's largest piece as the backbone and takes its 5' phosphates off.
`protocol_leave_out(path, "cip")` puts the phosphates back.

A PCR lane copies the part out of its map instead of cutting it out.
`protocol_new` and `protocol_add` make one with `pcr: true`, and `protocol_copy_out` turns a cut lane into one.
Its steps are the map, the PCR (`pcr_donor`), the cut that trims the primers' tails (`cut_pcr`) and the pick (`insert`).
The PCR's enzymes go on the primers' tails, the forward primer's first, and the cut uses the same enzymes.
The run reports the PCR's `primers`, each with its sequence and where it anneals, and the `primer_pairs` it designed, best first.
The primers are designed anew on each run until `protocol_primers` keeps a pair, after which every run checks them against the part.
`tm` aims the design at a melting temperature for the part of each primer that anneals, 60 °C unless set, and the run reports each primer's own `tm`.

A join's `byproducts` list the other molecules its tube makes, such as the insert the other way round or the backbone closed on itself.
Each one has a `description` in plain words and the `pieces` it holds in order, each with `turned_round`.
`can_choose` is true when it holds every piece in order, so `protocol_turn` can make it the product.

```rhai
let file = "egfp_into_vector.json";
protocol_new(file, "vector.dna", "donor.dna", #{ cut_out: "EGFP", enzymes: "SpeI,SalI" });
// A second donor's EGFP goes into the same new plasmid.
let added = protocol_add(file, "donor2.dna", #{ cut_out: "EGFP", enzymes: "SpeI,SalI" });
let run = protocol_put(file, "insert_2", "product");
print(run.changed);
for step in run.steps {
    print(`${step.id}\t${step.status}`);
}
```

### Sequences and oligos

| Function | What it does |
|---|---|
| `revcomp(seq) -> string` | The reverse complement. |
| `complement(seq) -> string` | The complement, in the same direction. |
| `gc(seq) -> float` | GC content in percent. |
| `translate(seq) -> string` | The protein from the first base, with the standard genetic code. |
| `translate(seq, table) -> string` | The protein from the first base, with the NCBI genetic code of this table number. |
| `tm(seq) -> float` | The nearest-neighbour melting temperature at SnapGene's default conditions. |
| `primer_report(seq) -> map` | Length, tm, gc_percent, self_dimer, three_prime_dimer, hairpin_stem and warnings for an oligo. |

`tm` uses the parameters of SantaLucia (1998) with 50 mM Na+ and 250 nM oligo.
These are the conditions SnapGene assumes by default.

### Utilities

| Function | What it does |
|---|---|
| `today() -> string` | Today's date as YYYY-MM-DD. |
| `join(items, separator) -> string` | The items of an array joined into one string. |
| `to_json(value) -> string` | Any value as JSON text. |
| `parse_json(text) -> any` | JSON text as a map, an array or a simple value. |

Rhai's own functions are also available.
They include `print`, `len`, `push`, `map`, `filter`, `sort`, `contains`, `split`, `to_upper` and `sub_string`.

## Examples

The folder [examples/scripts](../examples/scripts) holds complete scripts.

| Script | What it does |
|---|---|
| `library_summary.rhai` | One tab-separated line per file: length, topology, GC, CDS count and unique common cutters |
| `cloning_sites.rhai` | Enzymes that cut a plasmid once, inside a chosen feature such as the MCS |
| `golden_gate_check.rhai` | Parts that carry internal sites of a Golden Gate enzyme |

### Mark a motif in every file

This script adds a feature wherever a motif occurs, in every map of a folder.
It skips files that already carry the feature.

```rhai
let motif = "GGTCTC";
for path in files("maps", true) {
    let doc = open(path);
    if doc.has_feature("BsaI site") { continue; }
    let hits = doc.find(motif);
    for hit in hits {
        doc.add_feature("BsaI site", "misc_feature", hit.start, hit.end, hit.strand);
    }
    if hits.len() > 0 {
        doc.save();
        print(`${file_name(path)}: ${hits.len()} site(s)`);
    }
}
```

### Return structured results to an agent

A script can build maps and arrays and return them as its last value.
With `--json`, that value arrives as JSON in the `result` field.

```rhai
let rows = [];
for path in glob("maps/*.dna") {
    let doc = open(path);
    rows.push(#{ file: file_name(path), length: doc.length, cutters: doc.cut_counts("EcoRI,NotI") });
}
rows
```

```
slopgene run report.rhai --json | jq '.result[] | select(.cutters.NotI == 0) | .file'
```
