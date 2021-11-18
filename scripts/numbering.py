# Put numer in the name for StVs in each chr
# Usage: python3 numbering.py <input>.bed > <output>.bed
from sys import argv


# parse input bed
input_bed = []
with open(argv[1]) as bed:
    for line in bed:
        input_bed.append(line.split())

current_chr = input_bed[0][0]
current_bed = []
cnt = 0
for line in input_bed:
    if line[0] == current_chr:
        cnt += 1
        line[3] = '#{}:{}'.format(cnt, line[3])
        print('\t'.join(line))
    else:
        cnt = 0
        current_chr = line[0]
        print('\t'.join(line))
