import pandas as pd

path_blast = 'results/Ntab/blast_families_Ntab.csv'
path_blast_filtered = 'results/ata/blast_families_Ntab.filtered.csv'
params = {'min_len':50,'max_len':False,'min_distance':10,'max_q':1.2,'min_q':0.7,'min_pident':80,'min_qcov':80}

df = pd.read_csv(path_blast, sep='\t', header=None)
df.columns = ['qseqid','sseqid','qstart','qend','sstart','send','score','length','mismatch','gaps','gapopen','nident','pident','evalue','qlen','slen','qcovs']
print('initial:',len(df.index))
initial = len(df.index)

#filter by length
if(params['min_len']):
    df = df[df.qlen > params['min_len']]
print('Min len: ' + str(len(df.index)))
min_length = str(len(df.index))

if(params['max_len']):
    df = df[df.qlen < params['max_len']]
print('Max len: ' + str(len(df.index)))    
max_length = str(len(df.index))

#filter by query / subject length threshold
df = df[((df.length / df.qlen) >= params['min_q'])]
print('min treshold:',len(df.index))
min_treshold = str(len(df.index))

df = df[((df.length / df.qlen) <= params['max_q'])]
print('max treshold:',len(df.index))
max_treshold = str(len(df.index))

#filter by pident
df = df[(df.pident >= params['min_pident'])]
print('Min_pident: ' + str(len(df.index)))
min_pident = str(len(df.index))

#filter by qcov
df = df[(df.qcovs >= params['min_qcov'])]
print('Min qcov: ' + str(len(df.index)))
min_qcov = str(len(df.index))

#order sstart and send
df['new_sstart'] = df[['sstart','send']].min(axis=1)
df['new_ssend'] = df[['sstart','send']].max(axis=1)
df['sstart'] = df['new_sstart']
df['send'] = df['new_ssend']
df = df.drop('new_sstart',axis=1).drop('new_ssend',axis=1)
df = df.sort_values(by=['sseqid','sstart', 'send'])
df = df.reset_index(drop=True)
# sep by chr
dfs = {}
for seq in df.sseqid.unique():
    dfs[seq] = df[df.sseqid == seq]
print('done ordering sstart and send')

# filter overlapped 
rows = []
discard = []
total = len(df.index)
count = 0
curr = 0
for index, row in df.iterrows():
    count += 1
    curr_new = int(count * 100 * 1.0 / (total * 1.0))
    if curr_new != curr:
        curr = curr_new
        if curr_new % 1 == 0:
            print(curr_new)
    if index in discard:
        continue
    for k2, v2 in df.loc[index:,].iterrows():
        if abs(v2.sstart - row.sstart) > params['min_distance']:
            break
        if abs(v2.sstart - row.sstart) <= params['min_distance'] and abs(v2.send - row.send) <= params['min_distance']:
            discard.append(k2)
    rows.append(row)
print('done filtering overlap')

df = pd.DataFrame(rows)
print('Non overlapped: ' + str(len(df.index)))
non_overlapped = str(len(df.index))

df.to_csv(path_blast_filtered, index=None, sep='\t')

print('Initial: ' + str(initial))
print('Min len: ' + str(min_length))
print('Max len: ' + str(max_length))
print('Min treshold: ' + str(min_treshold))
print('Max treshold: ' + str(max_treshold))
print('Min pident: ' + str(min_pident))
print('Min qcov: ' + str(min_qcov))
print('Non overlapped: ' + str(non_overlapped))
print('Saved: ' + path_blast_filtered)
