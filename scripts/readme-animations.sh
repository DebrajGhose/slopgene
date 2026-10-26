#!/usr/bin/env bash
# Record the README's animations with the desktop app, from the made-up plasmids in examples/maps.
#
#   bash scripts/readme-animations.sh
#
# The app plays each script of clicks into its own window and saves the frames as an animated GIF
# in docs/images. CONTRIBUTING.md describes --play and --record under "Checking the app without a
# person at the screen". The points in the scripts are window coordinates for a 1000 by 760 point
# window at one pixel per point, read off screenshots taken with --screenshot, so a change to the
# layout means reading them again. A recording takes about as long as the animation plays.
set -euo pipefail
cd "$(dirname "$0")/.."

app=${SLOPGENE_GUI:-target/release/slopgene-gui}
if [ ! -x "$app" ] && [ ! -x "$app.exe" ]; then
    cargo build --release -p slopgene-gui
fi
size=1000x760
mkdir -p docs/images

# Cutting GFP out of pGFP and putting it into pDoom. The gene is clicked on its map, Clone opens
# the bench, the cut zooms into the map, the two enzymes nearest the gene are picked in the panel,
# and the join shows the new plasmid with its bases underneath, where a stretch is selected.
"$app" examples/maps/pDoom.dna examples/maps/pGFP.dna --view map --size "$size" \
    --record docs/images/clone-by-cutting.gif --play "
    wait 1.2
    click 508,345        # the GFP gene on the map
    wait 1.5
    click 670,21         # Clone, in the top bar
    wait 2.5
    click 373,385        # the cut of pGFP, which zooms into its map
    wait 2.2
    click 543,505        # EcoRI, before GFP
    wait 1.5
    click 553,587        # BamHI, after GFP
    wait 2.5
    key Escape           # back to the bench
    wait 3.0
    click 245,610        # the join, which zooms into the new plasmid
    wait 2.5
    drag 95,600,385,600  # a stretch of the bases under the map
    wait 3.0
"

# Copying mCherry out of pCherry by PCR and putting it into pDoom. The lane is switched from a cut
# to a PCR from the card's menu, the Copy card zooms into the map with the primers, one enzyme is
# picked for each primer's tail, and the bench scrolls down to the join.
"$app" examples/maps/pDoom.dna examples/maps/pCherry.dna --view map --size "$size" \
    --record docs/images/clone-by-pcr.gif --play "
    wait 1.2
    click 508,400        # the mCherry gene on the map
    wait 1.5
    click 670,21         # Clone, in the top bar
    wait 2.5
    click 442,303        # the menu of the pCherry cut
    wait 1.5
    click 544,381        # Copy the part out by PCR instead
    wait 2.5
    click 373,385        # the Copy card, which zooms into the map with the primers
    wait 2.5
    click 626,593        # EcoRI for the forward primer's tail
    wait 2.0
    scroll 900,150,300   # the panel, down to the reverse primer's enzymes
    wait 1.0
    click 543,584        # BamHI for the reverse primer's tail
    wait 2.5
    key Escape           # back to the bench
    wait 2.5
    scroll 250,400,400   # the bench, down to the join
    wait 1.5
    click 245,625        # the join, which zooms into the new plasmid
    wait 3.5
"
ls -l docs/images/*.gif
