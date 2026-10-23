// Check parts for internal Type IIS sites before Golden Gate assembly ("domestication").
// A part must not contain the assembly enzyme's site except where it is meant to be.
//
//   slopgene run examples/scripts/golden_gate_check.rhai --arg dir=parts --arg enzyme=BsaI

let dir = if "dir" in ARGS { ARGS.dir } else { "." };
let enzyme_name = if "enzyme" in ARGS { ARGS.enzyme } else { "BsaI" };

let flagged = 0;
for path in files(dir, true) {
    let doc = open(path);
    let sites = doc.sites(enzyme_name);
    if sites.len() > 0 {
        flagged += 1;
        let places = join(sites.map(|s| `${s.start}..${s.end} (${s.strand})`), ", ");
        print(`${file_name(path)}: ${sites.len()} ${enzyme_name} site(s) at ${places}`);
    }
}
print(`${flagged} file(s) carry ${enzyme_name} sites`);
