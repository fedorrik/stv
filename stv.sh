#!/bin/bash
# Pipeline creates StV bed file and StV stats table from monomeric AS-HOR bed file

if [[ $# -ne 1 ]]
then
    echo "usage: ./stv <monomeric track>"
    exit 1
fi

# put here path to dir with program
path='/storage/Bioinformatic/centromere/AnVIL/fedorrik/stv'

# take live HORs only
python3 $path/scripts/live_HORs_filter.py $1 > AS_liveHORs.bed

# monomers to StVs
python3 $path/scripts/mon2stv.py AS_liveHORs.bed > stv_row.bed

# stats file
python3 $path/scripts/bed2stat.py stv_row.bed > stats.tsv

# coloring
python3 $path/scripts/coloring.py stv_row.bed > stv_colored.bed

# numbering
python3 $path/scripts/numbering.py stv_colored.bed > stv.bed

# add descriotion line
sed -i "1 i\track name=\"ASat StV\" description=\"ASat HORs Structural Variants\" itemRgb=\"On\" visibility=\"1\"" stv.bed

rm stv_colored.bed
rm AS_liveHORs.bed
