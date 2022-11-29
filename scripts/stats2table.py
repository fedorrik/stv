from sys import argv
import pandas as pd


f = open(argv[1]).readlines()
dfs = {}
# fill dfs of each contig in the dict 
for line in f:
    line = line.strip()
    if '\t' not in line: # header
        contig = line.split(':')[0]
        dfs[contig] = pd.DataFrame(columns=[contig])
    elif 'total' in line: # line with total
        continue
    else: # common
        stv, cnt = line.split('\t')
        dfs[contig].loc[stv, contig] = cnt
# concat dfs
df = pd.concat(dfs.values(), axis=1).astype('float').fillna(0)
df.index.name = 'stv'
# sort values
df['sum'] = df.sum(axis=1)
df = df.sort_values('sum', ascending=False).drop('sum', axis=1)
# write
df.to_csv('stv_stats.table', sep='\t')

