// Summarise every sequence file in a folder as tab-separated text.
//
//   slopgene run examples/scripts/library_summary.rhai --arg dir=plasmids --arg recursive=true
//
// Columns: file, length, topology, GC %, CDS count, and the common enzymes that cut once.
// Backtick strings interpolate ${...} but do not interpret escapes, so TAB holds a tab.

const TAB = "\t";
let dir = if "dir" in ARGS { ARGS.dir } else { "." };
let recursive = "recursive" in ARGS && ARGS.recursive == "true";

print("file\tlength\ttopology\tgc\tcds\tunique_cutters");
for path in files(dir, recursive) {
    let doc = open(path);
    let cds = doc.features_of_type("CDS").len();
    let cutters = join(doc.unique_cutters("common").map(|c| c.enzyme), ",");
    print(`${file_name(path)}${TAB}${doc.length}${TAB}${doc.topology}${TAB}${doc.gc.to_int()}${TAB}${cds}${TAB}${cutters}`);
}
