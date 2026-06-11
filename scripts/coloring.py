#!/usr/bin/env python3
# Fill column 8 in the STV bed file. "Clever" random coloring.
import argparse
from collections import Counter
from random import choice, randint
from re import sub


def count_uniq_stv(stv_bed):
    all_stv = [sub(r'{\d+}', '{n}', line[3]) for line in stv_bed]
    return dict(Counter(all_stv).most_common())

def stv_clever_coloring(stv_bed, stv_stat_dict, highlight_common=False):
    bright_colors = ['255,0,0', '0,255,0', '0,0,255', '255,255,0', '0,255,255', '255,0,255', '130,255,0', '0,130,255', '130,0,255', '255,130,0', '0,255,130', '255,0,130',  '0,0,0']
    colors = {}
    for stv in stv_stat_dict:
        is_common = stv_stat_dict[stv] > len(stv_bed) / 10
        if is_common == highlight_common:
            colors[stv] = choice(bright_colors)
        else:
            colors[stv] = '{},{},{}'.format(randint(150, 220), randint(150, 220), randint(150, 220))
    for line in stv_bed:
        name = sub(r'{\d+}', '{n}', line[3])    # (_6/4_5){3} -> (_6/4_5){n}
        line[8] = colors[name]
    return stv_bed

def name_reverser(name):
    name = name.replace('(5_6/4_)', '(X)')
    name_splitted = name.split('.')
    name_digits = '.'.join(name_splitted[1:])
    # Split name first by '_' and than by '-'. Than reverse both. So S1C12H1L.8-4_7-1 --> [['1', '7'], ['4', '8']]
    items = list(reversed([list(reversed(i.split('-'))) for i in name_digits.split('_')]))
    # cen1 inversion issues 
    for i in range(len(items)):
        for j in range(len(items[i])):
            if 'X' in items[i][j]:
                items[i][j] = '{}{}{}'.format(items[i][j].split('}')[1], items[i][j].split('}')[0], '}')
    # join it back
    reversed_name = '_'.join(['-'.join(i) for i in items])
    reversed_name = reversed_name.replace('(X)', '(_6/4_5)')
    return '.'.join([name_splitted[0], reversed_name])


def parse_args():
    parser = argparse.ArgumentParser(
        description='Add STV colors to BED column 9 using one color system across the whole file.'
    )
    parser.add_argument('input_bed', help='Input STV BED file.')
    color_group = parser.add_mutually_exclusive_group()
    color_group.add_argument(
        '--highlight_rare',
        action='store_false',
        dest='highlight_common',
        help='Use bright colors for rare STVs and dim colors for common STVs. This is the default.',
    )
    color_group.add_argument(
        '--highlight_common',
        action='store_true',
        dest='highlight_common',
        help='Use bright colors for common STVs and dim colors for rare STVs.',
    )
    parser.set_defaults(highlight_common=False)
    return parser.parse_args()


args = parse_args()

# parse input bed
input_bed = []
with open(args.input_bed) as bed:
    for line in bed:
        input_bed.append(line.split())

# cen1 separately because of inversion
chr1_bed = [line for line in input_bed if line[0] == 'chr1']
inversion_names_backup = []
for line in chr1_bed:
    if line[5] == '-':
        inversion_names_backup.append(line[3])
        line[3] = name_reverser(line[3])

# Count and color all chromosomes together so the same STV in column 4
# always gets the same color, and the brightness rule is file-wide.
stv_stat = count_uniq_stv(input_bed)
stv_clever_coloring(input_bed, stv_stat, args.highlight_common)

cnt = -1
for stv in input_bed:
    if stv[0] == 'chr1' and stv[5] == '-':
        cnt += 1
        stv[3] = inversion_names_backup[cnt]
    print('\t'.join(stv))
