# Convert monomeric bed file into the HOR bed file
# Usage: python3 hor2stv.py <ASat live HORs>.bed > <StV>.bed
from collections import Counter
from random import choice, randint
from re import search, sub
from sys import argv


input_bed_path = argv[1]


def get_max_mon(live_mons_bed):
    numbers = []
    for mon in live_mons_bed:
        n = mon[3].split('.')[1]
        # problematic S2C8H1L.6/7s, so skip it 
        if 's' in n:
            numbers.append(7)
            continue
        # split hybrids on separate numbers
        if '/' in n:
            for hybrid_part in map(int, n.split('/')):
                numbers.append(hybrid_part)
            continue
        numbers.append(int(n))
    max_mon = str(max(numbers))
    return max_mon

def stv_namer(live_stv_name, mons_numbers, strand):
    if strand == '+':
        stv_name = mons_numbers[0]
        if mons_numbers[0].isdigit():
            i_prev = mons_numbers[0]
        elif '/' in mons_numbers[0] and 'S' not in mons_numbers[0]: # first mon is hybrid
            i_prev = 25
        else: # 8&12 in chr18
            i_prev = 25
        status = 'Closed'
        for i in mons_numbers[1:]:
            if i.isdigit():
                if int(i) == int(i_prev) + 1:
                    i_prev = i
                    status = 'ToBeClosed'
                else:
                    if status == 'ToBeClosed':
                        stv_name += '-{}_{}'.format(i_prev, i)
                    else:
                        stv_name += '_{}'.format(i)
                    i_prev = i
                    status = 'Closed'
            # hybrid like '4/7'
            elif '/' in i and 'S' not in i: 
                if status == 'ToBeClosed':
                    stv_name += '-{}_{}'.format(i_prev, i)
                else:
                    stv_name += '_{}'.format(i)
                i_prev = 25
                status = 'Closed'
            # 8&12 in chr18
            else: 
                if status == 'ToBeClosed':
                    stv_name += '-{}_{}'.format(i_prev, i)
                else:
                    stv_name += '_{}'.format(i)
                status = 'Closed'
                i_prev = '20'
        if status == 'ToBeClosed':
            stv_name += '-{}'.format(i)
    # reversed. strand == '-'
    else: 
        i_prev = 25
        stv_name = mons_numbers[0]
        if mons_numbers[0].isdigit():
            i_prev = mons_numbers[0]
        elif '/' in mons_numbers[0] and 'S' not in mons_numbers[0]: # first mon is hybrid
            i_prev = 25
        status = 'Closed'
        for i in mons_numbers[1:]:
            if i.isdigit():
                if int(i) == int(i_prev) - 1:
                    i_prev = i
                    status = 'ToBeClosed'
                else:
                    if status == 'ToBeClosed':
                        stv_name += '-{}_{}'.format(i_prev, i)
                    else:
                        stv_name += '_{}'.format(i)
                    i_prev = i
                    status = 'Closed'
            # hybrid like '4/7'
            elif '/' in i and 'S' not in i: 
                if status == 'ToBeClosed':
                    stv_name += '-{}_{}'.format(i_prev, i)
                else:
                    stv_name += '_{}'.format(i)
                i_prev = 25
                status = 'Closed'
            else:
                if status == 'ToBeClosed':
                    stv_name += '-{}_{}'.format(i_prev, i)
                else:
                    stv_name += '_{}'.format(i)
                status = 'Closed'
                i_prev = '20'
        if status == 'ToBeClosed':
            stv_name += '-{}'.format(i)
    stv_name = '{}.{}'.format(live_stv_name, stv_name)
    return(stv_name)

def cen1_5_19_compressor(stv_bed):
    # func is so fat because there are two case: normal (_6/4_5) and inversion in cen1 (5_6/4)
    corrected_bed = []
    # first loop _6/4_5_6/4_5... --> (X){n}, n > 1
    for line in stv_bed:
        if line[0] == 'chr1' and line[5] == '-':
            if search(r'[-_](5_6/4_){2,}', line[3]):
                old_name = line[3]
                new_name = inv_repeat_squizzer(old_name)
                line[3] = new_name
                corrected_bed.append(line)
            else:
                corrected_bed.append(line)
        else:
            if search(r'(_6/4_5){2,}[-_]', line[3]):
                old_name = line[3]
                new_name = repeat_squizzer(old_name)
                line[3] = new_name
                corrected_bed.append(line)
            else:
                corrected_bed.append(line)
    # one more loop to _6/4_5 --> (X){1}
    recorrected_bed = []
    for line in corrected_bed:
        if line[0] == 'chr1' and line[5] == '-':
            search_result = search(r'([-_])5_6/4_', line[3])
            if search_result:
                connector = search_result.groups()[0]
                old_name = line[3]
                new_name = sub(r'[-_]5_6/4_', '{}(X){}1{}'.format(connector, '{', '}'), old_name)
                line[3] = new_name
                recorrected_bed.append(line)
            else:
                recorrected_bed.append(line)
        else:
            search_result = search(r'_6/4_5([-_])', line[3])
            if search_result:
                connector = search_result.groups()[0]
                old_name = line[3]
                new_name = sub(r'_6/4_5[-_]', '(X){}1{}{}'.format('{', '}', connector), old_name)
                line[3] = new_name
                recorrected_bed.append(line)
            else:
                recorrected_bed.append(line)
    # (X) --> (_6/4_5)
    for line in recorrected_bed:
        if line[0] == 'chr1' and line[5] == '-':
            line[3] = line[3].replace('(X)', '(5_6/4_)')
        else:
            line[3] = line[3].replace('(X)', '(_6/4_5)')
    return recorrected_bed

def repeat_squizzer(name):
    result = search(r'(_6/4_5){2,}[-_]', name)
    all_repeats = result.group()
    n_repeats = int((len(all_repeats)-1)/6)
    connector = all_repeats[-1]
    name = name.replace(all_repeats, '(X){}{}{}{}'.format('{', n_repeats, '}', connector), 1)
    if search(r'(_6/4_5){2,}[-_]', name):
        name = repeat_squizzer(name)
    return name

def inv_repeat_squizzer(name):
    result = search(r'[-_](5_6/4_){2,}', name)
    all_repeats = result.group()
    n_repeats = int((len(all_repeats)-1)/6)
    connector = all_repeats[0]
    name = name.replace(all_repeats, '{}(X){}{}{}'.format(connector, '{', n_repeats, '}'), 1)
    if search(r'[-_](5_6/4_){2,}', name):
        name = inv_repeat_squizzer(name)
    return name

def print_bed(bed):
    for line in bed:
        print('\t'.join(line))


# parse input bed
input_bed = []
with open(input_bed_path) as bed:
    for line in bed:
        if line[:5] != 'track': # skip header
            input_bed.append(line.split())

# list of contigs
contigs = []
for line in input_bed:
    if line[0] not in contigs:
        contigs.append(line[0])
# MAIN LOOP
stv_bed = []
for contig in contigs:
    #print(contig)
    live_mons = []
    for line in input_bed:
        if line[0] == contig:
            live_mons.append(line)
    # get max mon
    max_mon = get_max_mon(live_mons)
    only_live_mons = []

    # fill stv list
    stvs = []
    mons_numbers = []
    n_prev = max_mon
    start = live_mons[0][1]
    end = live_mons[0][2]
    is_prev_max = True
    strand_prev = live_mons[0][5]
    for line in live_mons:
        #print(line)
        name, n = line[3].split('.')
        strand = line[5]
        mons_numbers.append(n)
        if strand == '+':

            # cut before mon which is after gap
            if int(line[1])-int(end) > 160:
                mons_numbers.pop()
                # check if mon isn't lonly mon
                if len(mons_numbers) > 0:
                    stv_name = stv_namer(name, mons_numbers, strand_prev)
                    stvs.append([contig, start, end, stv_name, '0', strand_prev, start, end, '0,0,0'])
                start = line[1]
                stv_name = []
                mons_numbers = [n]                

            # max mon (or max mon last in hybrid) THAN cut after it
            if n == max_mon or n[-2:] == '/{}'.format(max_mon) or n[-3:] == '/{}'.format(max_mon):
                end = line[2]
                stv_name = stv_namer(name, mons_numbers, strand_prev)
                stvs.append([contig, start, end, stv_name, '0', strand_prev, start, end, '0,0,0'])
                start = line[2]
                stv_name = []
                mons_numbers = []
                n_prev = n
                is_prev_max = True
                continue
            # first mon (or first mon is 1st in hybrid) AND previous wasn't max (or hybrid with max)
            elif ((n == '1' or n[:2] == '1/') and is_prev_max == False):
                mons_numbers.pop()
                # check if mon isn't lonly mon
                if len(mons_numbers) > 0:
                    stv_name = stv_namer(name, mons_numbers, strand)
                    stvs.append([contig, start, end, stv_name, '0', strand, start, end, '0,0,0'])
                start = line[1]
                stv_name = []
                mons_numbers = [n]
            end = line[2]
            n_prev = n
            is_prev_max = False
            strand_prev = strand


        else: # strand == '-'

            # cut before mon which is after gap
            if int(line[1])-int(end) > 160:
                # check if mon isn't lonly mon
                mons_numbers.pop()
                if len(mons_numbers) > 0:
                    stv_name = stv_namer(name, mons_numbers, strand_prev)
                    stvs.append([contig, start, end, stv_name, '0', strand_prev, start, end, '0,0,0'])
                start = line[1]
                stv_name = []
                if n == max_mon:
                    mons_numbers = [max_mon]
                else:
                    mons_numbers = [n]

            # first mon (or 1st is the last in hybrid) THAN cut after it
            if n == '1' or n[-2:] == '/1':
                end = line[2]
                stv_name = stv_namer(name, mons_numbers, strand_prev)
                stvs.append([contig, start, end, stv_name, '0', strand_prev, start, end, '0,0,0'])
                start = line[2]
                stv_name = []
                mons_numbers = []
                n_prev = n
                is_prev_max = True
                continue
            # fix BUG with cen1 inversion 
            elif contig == 'chr1' and n == max_mon and is_prev_max == False:
                mons_numbers.pop()
                stv_name = stv_namer(name, mons_numbers, strand)
                stvs.append([contig, start, end, stv_name, '0', strand, start, end, '0,0,0'])
                start = line[1]
                stv_name = []
                if n == max_mon:
                    mons_numbers = [max_mon]
                else:
                    mons_numbers = [n]
            # max mon (or max is 1st in hybrid) AND previous wasn't the first
            elif ((n == max_mon or n[:2] == '{}/'.format(max_mon) or n[-3:] == '{}/'.format(max_mon)) and is_prev_max == False and contig != 'chr1'):
                # check if mon isn't lonly mon
                mons_numbers.pop()
                if len(mons_numbers) > 0:
                    stv_name = stv_namer(name, mons_numbers, strand)
                    stvs.append([contig, start, end, stv_name, '0', strand, start, end, '0,0,0'])
                start = line[1]
                stv_name = []
                if n == max_mon:
                    mons_numbers = [max_mon]
                else:
                    mons_numbers = [n]
            end = line[2]
            n_prev = n
            is_prev_max = False
            strand_prev = strand

    if len(mons_numbers) > 0:
        stv_name = stv_namer(name, mons_numbers, strand)
        stvs.append([contig, start, end, stv_name, '0', strand, start, end, '0,0,0'])

    stvs_corrected = cen1_5_19_compressor(stvs)
    stv_bed += stvs_corrected

print_bed(stv_bed)
