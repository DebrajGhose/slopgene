"""Fill the methylation columns of SlopGene's enzyme table from REBASE.

REBASE sums up, for every enzyme, what Dam, Dcm, CpG, EcoBI and EcoKI methylation do to it when
a methylated base lies in its recognition site (https://rebase.neb.com/cgi-bin/damlist?eEcoRI).
This script reads that summary for each enzyme of the table and writes the five columns
`dam dcm cpg ecobi ecoki`, with the codes the table's header explains: c = cut as usual,
b = blocked, i = impaired, v = varies between sites, r = required for cutting, - = the methylase's
motif cannot overlap the site, ? = no data.

    py scripts/rebase-methylation.py crates/slopgene-core/data/enzymes.tsv --cache /tmp/rebase

The cache folder keeps one page per enzyme, so a second run downloads nothing. REBASE is a
free resource of New England Biolabs (Roberts et al., Nucleic Acids Res.); leave a pause between
requests, as the script does.
"""
import argparse
import pathlib
import re
import time
import urllib.parse
import urllib.request

COLUMNS = ["dam", "dcm", "cpg", "ecobi", "ecoki"]
# REBASE's words for an overlap, and the table's code for each.
WORDS = {"cut": "c", "blocked": "b", "impaired": "i", "variable": "v", "some blocked": "v", "some impaired": "i"}
# Enzymes that cut only methylated DNA. REBASE's overlap summary says that methylated sites are
# cut, but not that unmethylated ones are not, so these come from the enzymes' own entries.
REQUIRED = {"DpnI": {"dam": "r"}, "MalI": {"dam": "r"}}
# Enzymes whose REBASE summary page is empty: the Dam block of BclI and DpnII is textbook
# knowledge, written in every supplier's table.
KNOWN = {"BclI": {"dam": "b"}, "DpnII": {"dam": "b"}}
# Each methylase's motif with the indices of the bases it methylates, on the top strand and
# opposite the bottom strand, in the order of COLUMNS.
MOTIFS = [("GATC", [1], [2]), ("CCWGG", [1], [3]), ("CG", [0], [1]), ("TGANNNNNNNNTGCT", [2], [11]), ("AACNNNNNNGTGC", [1], [10])]
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
PAIR = {"A": "T", "C": "G", "G": "C", "T": "A", "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
        "B": "V", "V": "B", "D": "H", "H": "D", "N": "N"}
CELL = re.compile(r"<td><font size=1>(.*?)</font></td>", re.S)
# The table's header, the same lines the program writes (EnzymeDb::to_tsv).
HEADER = """# SlopGene restriction enzymes: commercially available enzymes from {source}
# REBASE, the Restriction Enzyme Database, http://rebase.neb.com (Roberts RJ, Vincze T, Posfai J, Macelis D. Nucleic Acids Res.)
# Cut positions are boundaries counted from the first base of the site (0 = just before it,
# negative = upstream), on the top strand (top) and bottom strand (bottom), in top-strand coordinates.
# The last five columns say what Dam, Dcm, CpG, EcoBI and EcoKI methylation do to the enzyme when a
# methylated base lies in its site, from REBASE's methylation sensitivity pages: c = cut as usual,
# b = blocked, i = impaired, v = varies between sites, r = required for cutting, - = the methylase's
# motif cannot overlap the site, ? = no data.
# Regenerate with: slopgene enzymes import-rebase --emboss-e link_emboss_e --emboss-r link_emboss_r,
# and scripts/rebase-methylation.py for the methylation columns."""


def fetch(name, cache):
    """The summary page of one enzyme, from the cache or from REBASE."""
    page = cache / (name + ".html")
    if not page.exists() or page.stat().st_size < 500:
        url = "https://rebase.neb.com/cgi-bin/damlist?e" + urllib.parse.quote(name)
        request = urllib.request.Request(url, headers={"User-Agent": "SlopGene rebase-methylation.py"})
        with urllib.request.urlopen(request, timeout=60) as response:
            page.write_bytes(response.read())
        time.sleep(0.5)
    return page.read_text(encoding="latin-1")


def codes_of(html):
    """The five codes of a summary page, or None when the page has no summary."""
    overlaps_at = html.find("overlaps?")
    sensitivity_at = html.find("sensitivity?")
    if overlaps_at < 0 or sensitivity_at < 0:
        return None
    overlaps = CELL.findall(html[overlaps_at:sensitivity_at])[:5]
    effects = CELL.findall(html[sensitivity_at:sensitivity_at + 4000])[:5]
    if len(overlaps) != 5 or len(effects) != 5 or any(o not in ("y", "n") for o in overlaps):
        return None
    codes = []
    for overlap, effect in zip(overlaps, effects):
        if overlap == "n":
            codes.append("-")
        else:
            codes.append(WORDS.get(effect.strip().lower(), "?"))
    return codes


def can_overlap(motif, site):
    """Whether some sequence puts a base the methylase modifies inside the site, with the motif in
    either orientation."""
    pattern, top, bottom = motif
    forms = [(pattern, list(top) + list(bottom))]
    reverse = "".join(PAIR[b] for b in reversed(pattern))
    if reverse != pattern:
        n = len(pattern)
        forms.append((reverse, [n - 1 - k for k in list(top) + list(bottom)]))
    for text, marks in forms:
        for offset in range(-(len(text) - 1), len(site)):
            fits = all(not 0 <= offset + k < len(site) or set(IUPAC[text[k]]) & set(IUPAC[site[offset + k]])
                       for k in range(len(text)))
            if fits and any(0 <= offset + k < len(site) for k in marks):
                return True
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("table", type=pathlib.Path)
    parser.add_argument("--cache", type=pathlib.Path, required=True, help="folder for the downloaded pages")
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    lines = args.table.read_text(encoding="utf-8").splitlines()
    source = next((line.split(" from ", 1)[1] for line in lines if line.startswith("#") and " from " in line), "REBASE")
    out = [HEADER.format(source=source.strip())]
    missing = []
    for line in lines:
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[0] == "name":
            out.append("\t".join(fields[:7] + COLUMNS))
            continue
        name = fields[0]
        codes = codes_of(fetch(name, args.cache))
        if codes is None:
            # No summary on REBASE: say "cannot overlap" where the geometry rules it out, and
            # "no data" elsewhere.
            codes = ["-" if not can_overlap(motif, fields[1]) else "?" for motif in MOTIFS]
            missing.append(name)
        for column, code in {**KNOWN.get(name, {}), **REQUIRED.get(name, {})}.items():
            codes[COLUMNS.index(column)] = code
        out.append("\t".join(fields[:7] + codes))
    args.table.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {args.table}; no REBASE summary for: {', '.join(missing) or 'none'}")


if __name__ == "__main__":
    main()
