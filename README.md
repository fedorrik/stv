# StV

HOR Structrual Variant (StV) prediction using HOR-monomer annotation

Usage: ./stv.sh AS-HOR-monomers.bed

Config: insert path to dir with program into stv.sh
___

Output files:

• stv.bed - resulting file which can be put in the browser
  
• stv_row.bed - same as stv.bed but doesn't contain stv numbering, colors, the first description line

• stv_stats.tsv - contains the number of each stv in each chr/contig
___

Method:

1. Extract only "live HORs" from input bed.
2. Iterate through AS_liveHORs.bed merging monomers. StV is cut after the last monomer (the one with maximal number in a given array) or before the first monomer (or vice versa in the case of a "-" strand).
3. Squizze huge "6/4_5" dimer stretches in the cen1/cen5/cen19.
4. Count StV stats.
5. Add numbering and colors. 

___

[Here](https://github.com/fedorrik/stv_chm13) is special version for T2T-CHM13 assembly which was used in the ["Complete genomic and epigenetic maps of human centromeres"](https://www.biorxiv.org/content/10.1101/2021.07.12.452052v1) artical.
