# SlopGene guide for scripts and AI agents

SlopGene reads, edits and analyses DNA, RNA and protein files from the command line.
This guide covers what a script or an agent needs to use it reliably.
Run `slopgene <command> --help` for every option of a command.

## Conventions

- **Coordinates are 1-based and inclusive**, as in SnapGene and GenBank. `100..250` covers 151 bases.
- **Ranges can run through the origin** of a circular sequence. On a 3000 bp plasmid, `2950..50` covers 101 bases.
- **Strands** are `forward`, `reverse`, `both` or `none`.
- **Formats** are detected from the file contents. They are SnapGene (`.dna`, `.rna`, `.prot`), GenBank (`.gb`, `.gbk`), FASTA and SlopGene JSON.
- **File arguments** accept files, folders and glob patterns such as `"plasmids/**/*.dna"`. Quote patterns so the shell passes them through. SlopGene expands them itself on every platform.

## Output

Every command prints text for people by default.
Add `--json` to get one JSON document, or `--jsonl` to get one JSON object per line.
Use `--jsonl` when a command covers many files.

A failed command exits with status 1.
In JSON modes it prints `{"ok": false, "error": "..."}` to standard output.
In text mode the error goes to standard error.

## Reading

| Task | Command |
|---|---|
| Summarise files | `slopgene info "lab/**/*.dna" --jsonl` |
| Everything in one file | `slopgene show pUC19.dna --json` |
| Sequence of a region | `slopgene seq pUC19.dna --range 396..454` |
| Sequence of a feature, along its strand | `slopgene seq pUC19.dna --feature AmpR` |
| Protein of a CDS | `slopgene translate pUC19.dna --feature AmpR --json` |
| Features of a type | `slopgene features "lab/*.dna" --type CDS --json` |
| Primers with Tm and warnings | `slopgene primers construct.dna --json` |
| Open reading frames | `slopgene orfs construct.dna --min-codons 100 --json` |
| Motif search, with mismatches | `slopgene find GGTCTC lab/ -r --mismatches 1 --json` |
| Protein motif in six frames | `slopgene find DYKDDDDK construct.dna --protein` |

`slopgene show --json` returns the full document in the schema `slopgene.document/1`.
Its main keys are `name`, `molecule`, `topology`, `length`, `sequence`, `features` and `primers`.
Each feature has `name`, `type`, `strand`, `start`, `end`, `length`, `color`, `segments` and `qualifiers`.
Qualifiers map a key to a list of values.

## Restriction enzymes

| Task | Command |
|---|---|
| Unique cutters among common enzymes | `slopgene sites pUC19.dna --unique` |
| Chosen enzymes, with positions | `slopgene sites pUC19.dna -e EcoRI,BamHI,HindIII --json` |
| Which enzymes cut nowhere | `slopgene sites insert.dna -e 6+ --noncutters` |
| Digest, with fragment ends | `slopgene digest pUC19.dna -e EcoRI,PvuII --json` |
| Digest and virtual gel | `slopgene digest pUC19.dna -e EcoRI --gel gel.png --uncut` |
| Enzyme details | `slopgene enzymes info BsaI` |
| Enzyme sets | `slopgene enzymes sets` |

Enzyme arguments take names (`EcoRI,BamHI`) or a set name.
The sets are `common`, `6+`, `8+`, `neb`, `thermo`, `golden-gate`, `type-iis`, `blunt` and `all`.
Site positions are top-strand cut positions, as SnapGene labels them.
The enzyme data comes from REBASE.

Methylation counts.
A plasmid grown in the usual E. coli strains is methylated by Dam (GATC), Dcm (CCWGG) and EcoKI, and an enzyme whose site holds a methylated base may not cut there: XbaI, ClaI and BclI are blocked by Dam where a GATC overlaps their site, StuI by Dcm, and DpnI cuts only Dam-methylated GATC.
SlopGene reads each file's methylation state (`.dna` files store it, and imported GenBank and FASTA DNA is taken to be Dam+ Dcm+ EcoKI+, as a new file is), leaves the blocked sites out of `sites`, `digest`, `excise`, cloning and protocols, and lists them under `blocked` with the reason.
A PCR product is unmethylated, so a site that methylation blocks on the template cuts on the product.
`slopgene edit plasmid.dna methylation "dam- dcm-"` changes the state (`none` clears it and `ecoli` sets all three), `info` and `show` print it, and `enzymes info` gives each enzyme's response to Dam, Dcm, CpG, EcoBI and EcoKI as REBASE records it.
To use a newer REBASE release, run `slopgene enzymes import-rebase` and point the environment variable `SLOPGENE_ENZYMES` at the result.

## Primers and PCR

| Task | Command |
|---|---|
| Tm and quality checks | `slopgene tm GTAAAACGACGGCCAGT --json` |
| Simulate a PCR | `slopgene pcr template.dna --fwd ATGGTG... --rev TTACTT... -o product.dna` |
| Use primers stored in the file | `slopgene pcr template.dna --fwd "M13 fwd" --rev "M13 rev"` |
| Design primers for a region | `slopgene design-primers template.dna --range 120..980 --fwd-tail GGTCTCA` |

Melting temperatures and primer checks follow Primer3.
Tm is SantaLucia's (1998) nearest-neighbour model with Primer3's salt correction, at Primer3's default conditions: 50 mM monovalent salt, 1.5 mM Mg2+, 0.6 mM dNTPs and 50 nM oligo.
Change them with `--na`, `--mg`, `--dntp` and `--oligo`.
`tm` reports Primer3's `self_any` and `self_end` pairing scores (Primer3 accepts up to 8 and 3) and the `end_stability` of the last five bases in kcal/mol.
`design-primers` and the PCR steps of protocols design primers as Primer3 does: every length from 18 to 27 bases is judged by Primer3's limits (Tm within 3 °C of the target, GC 20 to 80%, no run over 5, the two pairing scores, no other binding site on the template) and ranked by Primer3's penalty, one point per degree from the target and one per base from 20.
When no primer passes, the least faulty is returned with its faults as `warnings`.
Primer3's thermodynamic hairpin and dimer model is not used; the pairing scores are Primer3's alignment scores.
PCR products carry the primers' own bases, so 5' tails appear in the product.

## Cloning

| Task | Command |
|---|---|
| The enzymes that cut nearest a region on each side | `slopgene excise donor.dna --feature GFP --json` |
| The same, with how often each cuts the vector | `slopgene excise donor.dna --range 1201..1950 --vector pUC19.dna` |
| What a chosen pair cuts out | `slopgene excise donor.dna --feature GFP --left NcoI --right EagI --vector pUC19.dna` |
| Restriction cloning | `slopgene clone pUC57.dna insert.dna -e SalI,KpnI -o product.dna` |
| Gibson, HiFi or In-Fusion assembly | `slopgene assemble backbone.dna part1.dna part2.dna -o product.dna` |
| Golden Gate assembly | `slopgene golden-gate pGGA.dna part1.dna part2.dna -e BsaI -o product.dna` |

`excise` starts from the region to cut out, given as `--range` or `--feature`.
It lists the enzymes that cut nearest the region on each side, nearest first.
Each enzyme appears once on each side, at its cut nearest the region.
`--count` sets how many enzymes each side lists (8 unless given).
`-e` chooses the enzymes to consider, a set or a list (`common` unless given).
Nothing is paired or ranked, so the caller chooses one enzyme before the region and one after it.
The JSON has `before` and `after`, each entry with `enzyme`, `cut`, `distance`, `end` and `cuts` (how often the enzyme cuts the map).
With `--vector`, each entry also has `vector_cuts`.
`inside` names the enzymes that cut inside the region, which would cut it too.

`--left A --right B` shows what that pair cuts out instead.
The JSON's `excision` has `enzymes`, `left` and `right` (each with `enzyme`, `cut`, `extra` and `end`), `start`, `end`, `length`, `directional` and `other_cuts`.
With `--vector`, it also has `vector` (`cuts`, `dropped`, `problem`).
A pair that cannot cut the region out is an error that says why.
Pass the chosen pair to `protocol new -e` to build the protocol.

`clone` cuts both molecules with the same enzymes.
It keeps the vector's largest fragment as the backbone and every insert fragment whose ends fit it.
Both orientations are tried, so an insert with identical ends gives two products, saved as `product_1.dna` and `product_2.dna`.
`assemble` takes the fragments in order and joins each to the next through an end overlap of at least 15 bp (`--min-overlap`).
It closes the circle unless `--linear` is given.
`golden-gate` takes the destination vector first and the parts in any order.
Products keep the features of their parts.

Each product's JSON has `name`, `length`, `topology`, `features` and `notes`.
The notes describe what happened, such as `insert 721 bp, reversed`.
Without `-o`, the JSON also carries the product's `sequence`.

## Protocols

A protocol builds a molecule step by step and keeps every choice in a file.
Each run starts again from the maps, so a protocol gives the new product after a map is corrected.

| Task | Command |
|---|---|
| Write a restriction-cloning protocol and run it | `slopgene protocol new clone.json --vector pUC19.dna --insert donor.dna --cut-out GFP -e EcoRI,BamHI` |
| The same, with the insert copied out of its map by PCR | `slopgene protocol new clone.json --vector pUC19.dna --insert donor.dna --cut-out GFP -e EcoRI,BamHI --pcr` |
| Run a protocol | `slopgene protocol run clone.json --json` |
| Settle a step that waits for a choice | `slopgene protocol run clone.json --choose insert=2 --save` |
| Save the molecules the steps make | `slopgene protocol run clone.json --out-dir steps/` |

`protocol new` writes eight steps.
It opens the vector, cuts it, keeps its largest fragment and dephosphorylates it.
It then opens the donor, cuts it, keeps the insert and ligates the two.
`--cut-out` names the part of the donor to cut out, as a feature's name or a range such as `389..1108`.
The insert is then the fragment that holds that part, whatever enzymes cut it.
`--replace` names the part of the vector that the insert replaces, in the same way.
`-e` may be left out, and the cuts then wait for their enzymes.
`--vector-enzymes` cuts the vector with other enzymes than the insert, such as BglII for a BamHI insert.
`--keep-phosphates` leaves the dephosphorylation out.
`--pcr` copies the insert out of its map by PCR instead of cutting it out.
The donor is then followed by a PCR step (`pcr_donor`), a cut that trims the primers' tails (`cut_pcr`) and the pick.
The enzymes go on the primers' 5' tails, the forward primer's first, and the cut uses the same enzymes.
The file names each map relative to its own folder when the map lies inside it.

Twelve commands change a protocol file the way the app's bench does with the mouse.
Each one writes the file back and runs it.
A change that cannot be made is an error that says why, and the file stays as it was.
The JSON is the run's JSON with `changed`, a few words on what changed, and `saved`.

| Task | Command |
|---|---|
| Add a map to cut a piece from | `slopgene protocol add clone.json donor2.dna --cut-out RFP -e EcoRI,BamHI` |
| Change a cut's part, its enzymes or both | `slopgene protocol cut clone.json cut_donor --part GFP -e SpeI,SalI` |
| Keep another piece of a cut, counted from 1 along the cut map | `slopgene protocol keep clone.json insert --piece 2` |
| Put a kept piece into a map or a join | `slopgene protocol put clone.json insert_2 --into product` |
| Take a piece out of a join | `slopgene protocol take-out clone.json insert_2 --join product` |
| Take a map out with every step that works on it | `slopgene protocol remove-map clone.json donor_2` |
| Leave out a step that has one input | `slopgene protocol leave-out clone.json cip` |
| Take a piece's 5' phosphates off on its way into a join | `slopgene protocol dephosphorylate clone.json backbone --join product` |
| Say which way a piece goes into a join | `slopgene protocol turn clone.json insert --join product --way reverse` |
| Copy a map's part out by PCR instead of cutting it out | `slopgene protocol copy-out clone.json donor` |
| Cut it out again, without the PCR | `slopgene protocol cut-out clone.json donor` |
| Keep a pair of the primers a PCR designed | `slopgene protocol primers clone.json pcr_donor --pair 1` |
| Use primers of your own, or design them anew on each run | `slopgene protocol primers clone.json pcr_donor --fwd TAAGCAGAATTC... --rev ...`, `--clear` |
| Aim the primer design at a melting temperature | `slopgene protocol primers clone.json pcr_donor --tm 58` |
| Draw the bench as a picture | `slopgene protocol render clone.json -o clone.png` |

`add` gives its steps new ids, such as `donor_2`, `cut_donor_2` and `insert_2`.
Without `-e`, the new cut waits for its enzymes.
`add --pcr` copies the part out by PCR, with the steps `pcr_donor_2`, `cut_pcr_2` and `insert_2`.
`cut --no-part` forgets the part, and `cut -e ""` takes the enzymes off.
`cut` on a PCR step changes the part it copies and the enzymes on its primers' tails, and the cut that trims the tails changes with it.
New enzymes on the insert's cut carry over to a vector cut that used the same enzymes.
The cut of a second insert keeps its own enzymes.
New enzymes on a vector's cut carry over to the PCRs whose tails had the same enzymes, since a tail can carry any site, while a donor's cut keeps its own.
`copy-out` and `cut-out` take any step of the map's lane, and switch the lane between a cut and a PCR with the same part and enzymes.
`primers --pair N` keeps pair N of the pairs `protocol run` lists, 1 being the best.
`--fwd` and `--rev` give two primers of your own, written 5' to 3' with their tails, and `--clear` designs them anew on each run.
`--tm` sets the melting temperature the design aims for, for the part of each primer that anneals, 60 °C unless set, and designs the primers anew at it, or with `--pair` keeps a pair designed at it.
`protocol new --pcr --tm 58` sets it from the start.
`render` draws the protocol as the bench the app shows: each map as a small ring with its features, its cut under it with the piece kept in green and the rest in grey, a PCR as its product, and the join with the new plasmid, with the enzymes and treatments written on the arrows.
`-o` names an .svg or .png file (by default the protocol's name with .svg), and `--dark` and `--scale` work as for `slopgene render`.
`put --into` a map cuts that map with the piece's own enzymes, so the ends match.
It keeps the map's largest piece as the backbone, takes its 5' phosphates off, and adds a join of the two.
`--keep-phosphates` leaves the phosphates on.
`put --into` a join adds the piece after the join's other pieces.
`turn --way` takes `forward` (as given), `reverse` (turned round) or `either`.
With `either`, the join waits for a choice when both ways fit.
The first piece of a join sets the direction, so it cannot be turned.

A protocol file is JSON with the schema `slopgene.protocol/1`:

```json
{
  "schema": "slopgene.protocol/1",
  "name": "geneX into pUC19",
  "steps": [
    {"id": "vector", "type": "open", "params": {"path": "pUC19.dna"}},
    {"id": "cut_vector", "in": ["vector"], "type": "digest", "params": {"enzymes": ["EcoRI", "BamHI"]}},
    {"id": "backbone", "in": ["cut_vector"], "type": "pick", "params": {"largest": true}},
    {"id": "cip", "in": ["backbone"], "type": "dephosphorylate"},
    {"id": "donor", "type": "open", "params": {"path": "donor.dna"}},
    {"id": "cut_donor", "in": ["donor"], "type": "digest", "params": {"enzymes": ["EcoRI", "BamHI"], "around": "geneX"}},
    {"id": "insert", "in": ["cut_donor"], "type": "pick", "params": {"contains": "geneX"}},
    {"id": "product", "in": ["cip", "insert"], "type": "ligate"}
  ]
}
```

Each step takes the outputs of earlier steps through `in`.
A step id names all of that step's outputs, and `id#2` names its second output.

| Type | Parameters | What it does |
|---|---|---|
| `open` | `path`, or `document` (a `slopgene.document/1` object) | starts from a map |
| `digest` | `enzymes`, `around` (a feature's name or a range) | cuts one molecule and lists its fragments in order along it |
| `amplify` | `region` (a feature's name or a range), `forward_enzyme`, `reverse_enzyme`, `forward_primer`, `reverse_primer` | copies the region by PCR, with the enzymes' sites on the primers' 5' tails |
| `pick` | `index`, `size` with `tolerance`, `ends`, `contains`, `holds` (a range of the cut molecule), `largest`, `smallest` | keeps the one fragment that meets every criterion given |
| `dephosphorylate` | none | removes both 5' phosphates of a linear fragment |
| `flip` | none | turns a molecule end to end |
| `ligate` | `circular` (true unless set), `orientations` | joins its inputs in order and closes the product into a circle |

A digest's `around` names the part of the molecule the cut is about.
The run says which fragment holds that part, and warns when an enzyme cuts inside it.
A pick's `contains` matches the fragment that carries the most of the named feature, so a few of its bases in another fragment's overhang do not count.
A digest without enzymes waits for them.
A PCR without a region waits for it.
A PCR designs primers for its region anew on each run until `forward_primer` and `reverse_primer` are set, which `protocol primers` does, and then checks them against the region on every run.
Its `target_tm` is the melting temperature the design aims for, 60 unless set.
Its JSON lists the `primers` used, each with its sequence, its melting temperature and where it anneals, and the `primer_pairs` designed, best first.
Each tail is six spare bases and the enzyme's site, so the enzyme can cut near the end of the product.
A PCR product carries no 5' phosphates, as primers are made without them, and the cut that trims its tails leaves phosphates on the piece kept.
A PCR lane in a file reads like this:

```json
    {"id": "pcr_donor", "in": ["donor"], "type": "amplify", "params": {"region": "geneX", "forward_enzyme": "EcoRI", "reverse_enzyme": "BamHI"}},
    {"id": "cut_pcr", "in": ["pcr_donor"], "type": "digest", "params": {"enzymes": ["EcoRI", "BamHI"], "around": "geneX"}},
    {"id": "insert", "in": ["cut_pcr"], "type": "pick", "params": {"contains": "geneX"}}
```

A step chooses by itself only when exactly one option fits.
A pick that several fragments match waits for a choice.
A ligation whose pieces join more than one way round waits too.
The JSON of such a step lists its `choices`, each with the molecules it gives and their junctions.
`--choose STEP=N` settles it, and `--save` writes the choice into the file as `index` or `orientations`.
Every run checks a pinned `index` against the pick's criteria, so a changed map cannot slip through.

Ligation checks every junction as ligase would.
The two ends must be the same kind and length, with overhangs that pair.
At least one of the two must carry a 5' phosphate.
A dephosphorylated vector therefore cannot close on itself.

A ligation also lists the other molecules its tube makes as `byproducts`.
Ligase joins every pair of ends that fit, whatever direction the step chose.
An insert that fits both ways therefore also goes in the other way round.
A backbone that keeps its phosphates can close on itself.
A circle of several inserts can also form with fewer of them or with them in another order.
Every by-product holds the join's first piece, which is the backbone.
Only circles that hold the backbone give colonies.
Each by-product has a `description`, the `pieces` it holds in order (`from` and `turned_round`), its `molecule`, `junctions` and `parts`.
`can_choose` is true when a by-product holds every piece in order.
`protocol turn` then makes it the product.
Several copies of one insert can also join in a row.
The run does not list those.
For a circle with more than three inserts, the run lists only the other directions.

Each step in the JSON has a `status` of `done`, `needs_choice`, `failed` or `skipped`.
Each step also has a `headline` and `details`, which say in plain words what the step does and what came of it.
The run's `summary` tells the whole protocol in a few sentences.
A failed step has an `error` with a stable `code`, a `message` and `fixes`.
A skipped step names the step it waits for in `skipped_because`.
`finished` is true once every step is done.
Digests list their `cuts`, picks give the fragment they kept as `picked`, and ligations list their `junctions`.
A digest with a part gives it as `part` (`start`, `end` and the `fragment` that holds it).
A ligation lists the stretch of the product each input became as `parts` (`from`, `start`, `end`).
Positions count from 1.

| Code | Meaning |
|---|---|
| `ends.incompatible` | two ends that should join do not pair |
| `ends.no_phosphate` | neither end of a junction has a 5' phosphate |
| `digest.no_cut` | the enzymes do not cut the molecule (the fixes name enzymes that cut it once) |
| `digest.cuts_part` | a warning: an enzyme cuts inside the digest's part |
| `digest.no_part` | a warning: the digest's part names no feature or range of the molecule |
| `digest.blocked` | a warning: methylation blocks a site, so it is not cut (the fix suggests a PCR copy, which is unmethylated) |
| `digest.impaired` | a warning: methylation impairs a site's cut, or may block it |
| `pick.no_match` | no fragment meets the criteria |
| `pick.stale` | the pinned fragment no longer meets the criteria |
| `input.unknown` | an input names no step, or a step further down |
| `input.several` | a step that works on one molecule was given several |

## Verifying sequencing results

| Task | Command |
|---|---|
| Align one read to one map | `slopgene align map.dna read.gbk --json` |
| See the aligned bases too | `slopgene align map.dna read.gbk --show` |
| Draw the map with the read's differences | `slopgene align map.dna read.gbk --map read.svg` |
| Match many reads to many maps | `slopgene verify "results/*.gbk" --refs maps/ -r --jsonl` |

Reads can be whole-plasmid results (GenBank, FASTA, SnapGene) or Sanger reads.
They align in either orientation and through the origin of a circular map.
`verify` ranks every reference by the k-mers it shares with the read and aligns the read to the best one.
`--top N` (3 by default) lists N candidates: the aligned best and N - 1 others under `runners_up`.

Alignment works like SnapGene's "Align to Reference" and minimap2.
Shared 15-mers are chained in the order they occur along both sequences, and the stretches between them are aligned base by base.
A chain carries on across a deletion or insertion of any size, so a clone that lost or gained a whole element still aligns along its full length.
Changes of 50 bp or more are reported as structural changes under `events`, and identity is measured outside them.
A read that is nearly as long as a circular map, or that comes from a circular record, is treated as a whole plasmid: a change across the base where sequencing started is still found.

Each result has a `verdict`.

| Verdict | Meaning |
|---|---|
| `exact match` | No differences, and the read covers more than 98% of the map |
| `exact over the region read` | No differences, but the read covers only part of the map |
| `match (homopolymer indels only)` | At least 99% identity, and every difference is a homopolymer indel (see below) |
| `match with differences to inspect` | At least 99% identity, with at least one other difference |
| `structural change` | At least 99% identity outside one or more structural changes (see `events`) |
| `structural change in a sibling construct` | Structural changes, and 95% to 99% identity outside them |
| `near match: probably a sibling construct` | At least 95% identity, no structural change |
| `not this plasmid` | Below 95% identity, less than half the read aligned, or no shared sequence at all |

The other keys are `identity` (percent), `orientation`, `reference_start`, `reference_covered`, `reference_length`, `coverage` (percent), `read_length`, `read_aligned`, `whole_read`, `kmer_containment`, `events` and `differences`.

`align --show` writes out the aligned bases, 60 columns to a line (`--width N` changes this).
Each line of the reference sits over the matching line of the read.
A `|` between them marks bases that agree.
The positions at both ends of each line are 1-based.
A read on the reference's bottom strand is shown reverse-complemented.
Its positions count down.
A whole-plasmid read starts at the map's first base.
Each structural change takes one line in brackets in place of its bases.
With `--json`, `--show` adds `aligned`, which holds the two rows as strings with `-` for each gap: `aligned.reference` and `aligned.read`.
`differences` follow the alignment, so a whole-plasmid read lists them from the map's first base.

`align --map FILE` draws the reference map with the read as a ring around it, as SVG, or as PNG when the name ends in `.png`.
The ring is green where the read aligns, a thin red line where it lacks the map's bases, and absent beyond the read's ends.
Each difference is a mark on the ring in its kind's colour, labelled with its position and change.
Orange marks different bases, violet marks bases only the read has, and red marks bases the read lacks.
With `--json` the output then has `map`, the path written.

Every alignment's JSON has `method`, which names the method and its settings: `seed_length`, `chaining`, `alignment`, `scores` and `structural_min_bp`.
The text output states the method on its second line.

Each structural change in `events` has `type`, `length` (reference bases), `read_bases`, `features` and `description`.
It has `reference_start` and `reference_end`, or `after` for a change that adds bases at one point.

| Event type | Meaning |
|---|---|
| `deletion` | Reference bases the read lacks |
| `insertion` | Read bases found nowhere in the reference |
| `duplication` | Read bases that copy part of the reference, including a dimer that holds the whole plasmid twice |
| `inversion` | A reference stretch the read carries in the opposite orientation |
| `replacement` | A reference stretch the read has replaced by other sequence |

The description names the features a change touches, as in `1,500 bp deletion of reference 2,001..3,500 (removes WPRE)`.

Each difference has `type` (`substitution`, `insertion` or `deletion`), `position`, `reference`, `read`, `features`, `effects` and `homopolymer`.
`position` is the first reference base involved, or the base an insertion follows.
`effects` gives the protein change in each coding feature, such as `AmpR: L191V` or `lacZα: frameshift at codon 12`.
`homopolymer` is true for an insertion or deletion of one repeated base inside a run of four or more of that base.
Such indels are the usual nanopore sequencing error, so they deserve caution before a clone is rejected.

## Annotating

| Task | Command |
|---|---|
| Annotate one file from trusted maps | `slopgene annotate read.gbk --library maps/ -r -o read_annotated.dna` |
| Annotate many files in place | `slopgene annotate "results/*.gbk" --library maps/ -r --in-place` |
| Keep the library for reuse | `slopgene annotate read.gbk --library maps/ -r --save-library library.dna` |

The library is every feature of at least 30 bp (`--min-length`) in the library files.
A feature is added where it matches over its whole length with at least 95% identity (`--min-identity`).
The best-matching library feature wins where several match the same stretch.
A library feature of the same type and a different name on nearly the same stretch becomes an `also_called` qualifier of the winner.
A match is skipped where the target already has a feature of the same name over most of that stretch.
It is also skipped where the target has a feature of the same type over nearly the same bases.
Without `-o` or `--in-place`, the command only reports what it would add.

## Editing

Edits change the file in place unless `-o` names another file.
`--dry-run` reports the change without writing.
Writes are atomic, so a reader never sees a half-written file.
Features and primer binding sites move with every edit.

```
slopgene edit construct.dna insert --at 120 --seq GAATTC
slopgene edit construct.dna delete --range 120..125
slopgene edit construct.dna replace --range 300..300 --seq T
slopgene edit construct.dna add-feature --name "His tag" --type CDS --range 801..818 --strand forward
slopgene edit construct.dna update-feature "His tag" --color "#515faa"
slopgene edit construct.dna add-primer --name seqF --seq GTAAAACGACGGCCAGT
slopgene edit construct.dna origin --at 1500
slopgene edit construct.dna topology linear
slopgene edit construct.dna methylation "dam- dcm-"
```

A replacement of the same length is a point mutation, and nothing moves.
Feature names are matched exactly and without regard to case.
Use `#N` to pick the Nth feature when names repeat.

## Converting and drawing

```
slopgene convert "lab/*.dna" --to genbank --out-dir genbank/
slopgene convert result.gbk -o result.dna
slopgene render pUC19.dna -o pUC19.png
slopgene render pUC19.dna -o pUC19.svg --enzymes EcoRI,BamHI --dark
```

SnapGene files survive a read and write unchanged, byte for byte.
After an edit the writer drops SnapGene's history and cached data, which SnapGene rebuilds when it opens the file.

### Which format to save in

Save annotated work as SnapGene `.dna` (`.rna`, `.prot`).
It keeps everything SlopGene models, and SnapGene opens what SlopGene writes.
GenBank keeps features and their colours, but primers survive only as `primer_bind` features, and sticky ends are lost.
FASTA keeps only the name and the bases.
SlopGene JSON keeps everything but is read by SlopGene alone.

An edit written back into its own file never drops anything silently.
`slopgene edit` and `slopgene annotate --in-place` refuse with an error and write nothing when the file's format cannot hold the result, such as a feature added to a FASTA file.
Name a `.dna` output with `-o` instead.
In the desktop app, Save writes a `.dna` file beside such a file and leaves the original as it was.

## Scripts

Work that needs loops, conditions or several steps fits in a script.
SlopGene runs scripts in the Rhai language.

```
slopgene run report.rhai --arg dir=maps --json
slopgene run -e 'open("pUC19.dna").cut_counts("EcoRI,BamHI")' --json
slopgene run --api
```

`--arg key=value` sets `ARGS.key` inside the script.
With `--json`, the command prints `{"ok": true, "output": ..., "result": ...}`.
`output` holds what the script printed, and `result` holds the value of its last expression.
`slopgene run --api` lists every function with a one-line description.
[scripting.md](scripting.md) documents the language, the functions and complete examples.

## Updating

`slopgene update` says whether a newer version of SlopGene is released on GitHub, and `slopgene update --install` downloads it, checks it against the release's checksums and installs it the way this copy was installed.
With `--json` the answer has `current`, `latest`, `newer`, `published`, `page`, `installation` (`windows-installer`, `windows-portable`, `mac-bundle`, `deb`, `appimage`, `portable` or `development`) and `file`, the download for this system.
A `development` copy, built from the source tree, is not installed over: build it again.
The desktop app checks the same way when it starts, and its top bar then offers the update.

## The desktop app

`slopgene view construct.dna` opens files in the desktop app.
A protocol file opens in a protocol tab, which shows each step with its maps.
The protocol commands above make the same changes a person makes on that tab's bench.
The app reloads a file that changes on disk, as long as its window holds no unsaved edits.
A person can therefore watch the map while an agent edits the file with `slopgene edit`.

## Patterns for bulk work

Summarise a folder and keep the circular plasmids over 5 kb:

```
slopgene info lab/ -r --jsonl | jq -c 'select(.topology == "circular" and .length > 5000)'
```

Check which plasmids an enzyme pair cuts exactly once:

```
slopgene sites lab/ -r -e EcoRI,NotI --jsonl | jq -c '{name, cuts: [.enzymes[] | {(.name): .cuts}]}'
```

List the reads that need a closer look:

```
slopgene verify "results/*.gbk" --refs maps/ -r --jsonl | jq -c 'select(.verdict | startswith("exact") | not) | {read, reference, verdict}'
```

From Python, call the command and parse its JSON:

```python
import json, subprocess
out = subprocess.run(["slopgene", "show", "pUC19.dna", "--json"], capture_output=True, text=True, check=True)
doc = json.loads(out.stdout)
cds = [f for f in doc["features"] if f["type"] == "CDS"]
```
