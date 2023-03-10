# StV

HOR Structural Variant (StV) prediction using HOR-monomer annotation

Usage: ./stv.sh input.bed

Note: input.bed must be sorted (bedtools sort) and must contain the HOR monomers. You can generate it from sequence using [HumAS-HMMER](https://github.com/enigene/HumAS-HMMER) for example.

Config: insert path to dir with program into stv.sh
___

Output files:

• stv.bed - resulting file which can be put in the browser
  
• stv_row.bed - same as stv.bed but doesn't contain stv numbering, colors, the first description line

• stv_stats.tsv - contains the number of each stv in each chr/contig
___

Method:

1. Extract only "live HORs" from the input bed.
2. Iterate through AS_liveHORs.bed merging monomers. StV is cut after the last monomer (the one with maximal number in a given array) or before the first monomer (or vice versa in the case of a "-" strand).
3. Squeeze huge "6/4_5" dimer stretches in the human cen1/cen5/cen19.
4. Count StV stats.
5. Add numbering and colors. 

___

[Here](https://github.com/fedorrik/stv_chm13) is special version for T2T-CHM13 assembly which was used in the ["Complete genomic and epigenetic maps of human centromeres"](https://www.science.org/doi/10.1126/science.abl4178) article.
