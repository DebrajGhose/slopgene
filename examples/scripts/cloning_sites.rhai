// Which enzymes cut a plasmid exactly once, inside a chosen feature (for example the MCS)?
// Those are the enzymes that open the plasmid there without cutting it anywhere else.
//
//   slopgene run examples/scripts/cloning_sites.rhai --arg file=pUC19.dna --arg feature=MCS
//
// Add --arg set=6+ to consider every enzyme with a 6-base site, not just the common ones.

const TAB = "\t";
let doc = open(ARGS.file);
let target = doc.feature(ARGS.feature);
let set = if "set" in ARGS { ARGS.set } else { "common" };
print(`${doc.name}: enzymes cutting once, inside ${target.name} (${target.start}..${target.end})`);

// Whether a cut lies inside the feature, allowing for features that run through the origin.
fn inside(pos, start, end) {
    if start <= end { pos >= start && pos <= end } else { pos >= start || pos <= end }
}

let count = 0;
for cutter in doc.unique_cutters(set) {
    if inside(cutter.cut, target.start, target.end) {
        let e = enzyme(cutter.enzyme);
        print(`  ${cutter.enzyme}${TAB}${e.site_with_cuts}${TAB}cuts after ${cutter.cut}`);
        count += 1;
    }
}
print(`${count} enzymes`);
