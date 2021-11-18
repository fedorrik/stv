# find the most frequent live HOR in each contig and print only mons of this HOR
# Usage: python3 live_HORs_filter.py <AS-HOR>.bed > AS_liveHORs.bed
from collections import Counter
from sys import argv


input_bed_path = argv[1]

# parse input bed
input_bed = []
with open(input_bed_path) as bed:
    for line in bed:
        if line[:5] != 'track': # skip header
            input_bed.append(line.split())

contig = input_bed[0][0]
live_HORs = {}
contig_live_HORs = []
for line in input_bed:
    if line[0] == contig: # current contig
        
        # take all live HOR in list
        HOR_name = line[3].split('.')[0]
        if HOR_name[-1] == 'L':
            contig_live_HORs.append(HOR_name)

    else: # next contig
        # deal with previous cintig
        if len(contig_live_HORs) > 0:
            live_HORs[contig] = Counter(contig_live_HORs).most_common(1)[0][0]
        # start procces next contig
        contig = line[0]
        contig_live_HORs = []
        # procces first line
        HOR_name = line[3].split('.')[0]
        if HOR_name[-1] == 'L':
            contig_live_HORs.append(HOR_name)
live_HORs[contig] = Counter(contig_live_HORs).most_common(1)[0][0]

# print chosen HORs
for line in input_bed:
    if line[0] in live_HORs and line[3].split('.')[0] == live_HORs[line[0]]:
        print('\t'.join(line))
