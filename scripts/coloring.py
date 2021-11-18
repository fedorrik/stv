# Fill cillumn 8 in the stv bed file. "Clever" random coloring
# Usage: python3 bed2stat.py <input>.bed > <output>.bed
from collections import Counter
from random import choice, randint
from re import sub
from sys import argv, exit


def count_uniq_stv(stv_bed):
    all_stv = [sub(r'{\d+}', '{n}', line[3]) for line in stv_bed]
    return dict(Counter(all_stv).most_common())

def stv_clever_coloring(stv_bed, stv_stat_dict):
    bright_colors = ['255,0,0', '0,255,0', '0,0,255', '255,255,0', '0,255,255', '255,0,255', '130,255,0', '0,130,255', '130,0,255', '255,130,0', '0,255,130', '255,0,130',  '0,0,0']
    colors = {}
    for stv in stv_stat_dict:
        if stv_stat_dict[stv] > len(stv_bed) / 10:
            colors[stv] = '{},{},{}'.format(randint(150, 220), randint(150, 220), randint(150, 220))
        else:
            colors[stv] = choice(bright_colors)
    for line in stv_bed:
        name = sub(r'{\d+}', '{n}', line[3])    # (_6/4_5){3} -> (_6/4_5){n}
        line[8] = colors[name]
    return stv_bed

def name_reverser(name):
    name = name.replace('(5_6/4_)', '(X)')
    name_splitted = name.split('.')
    name_digits = '.'.join(name_splitted[1:])
    # Braindead line. Split name first by '_' and than by '-'. Than reverse both. So S1C12H1L.8-4_7-1 --> [['1', '7'], ['4', '8']]
    items = list(reversed([list(reversed(i.split('-'))) for i in name_digits.split('_')]))
    # shit with cen1 inversion
    for i in range(len(items)):
        for j in range(len(items[i])):
            if 'X' in items[i][j]:
                items[i][j] = '{}{}{}'.format(items[i][j].split('}')[1], items[i][j].split('}')[0], '}')
    # And now join this shit back
    reversed_name = '_'.join(['-'.join(i) for i in items])
    reversed_name = reversed_name.replace('(X)', '(_6/4_5)')
    return '.'.join([name_splitted[0], reversed_name])


# parse input bed
input_bed = []
with open(argv[1]) as bed:
    for line in bed:
        input_bed.append(line.split())

# cen1 separatly because inversion
chr1_bed = [line for line in input_bed if line[0] == 'chr1']
inversion_names_backup = []
for line in chr1_bed:
    if line[5] == '-':
        inversion_names_backup.append(line[3])
        line[3] = name_reverser(line[3])
# print coloured chr1
stv_stat = count_uniq_stv(chr1_bed)
cnt = -1
for stv in stv_clever_coloring(chr1_bed, stv_stat):
    if stv[5] == '-':
        cnt += 1
        stv[3] = inversion_names_backup[cnt]
    print('\t'.join(stv))
    
input_bed = [line for line in input_bed if line[0] != 'chr1']
current_chr = input_bed[0][0]
current_bed = []
for line in input_bed:
    chr = line[0]
    strand = line[5]
    if chr == current_chr:
        current_bed.append(line)
    else:
        stv_stat = count_uniq_stv(current_bed)
        for stv in stv_clever_coloring(current_bed, stv_stat):
            print('\t'.join(stv))
        current_chr = chr 
        current_bed = [line]
# for the last chr
stv_stat = count_uniq_stv(current_bed)
for stv in stv_clever_coloring(current_bed, stv_stat):
    print('\t'.join(stv))
