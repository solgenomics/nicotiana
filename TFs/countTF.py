import re
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(
    description='Path to the output from hmmscan formatted as parsable table')
parser.add_argument('-f', action="store", dest="f", type=str, required=True)
args = parser.parse_args()

res = []
p = re.compile(r"\s+")
with open(args.f) as f:
    for line in f.read().splitlines():
        if not line.startswith('#'):
            res.append(p.split(line, maxsplit=18))

df = pd.DataFrame(columns=['targetName', 'targetAccession', 'queryName', 'queryAccession',
                           'E-value1', 'score1', 'bias1', 'E-value2', 'score2', 'bias2',
                           'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'targetDescription'],
                  data=res)

TFs = set(df['targetName'])


output = []
for tf in TFs:
    num = np.sum(np.array(df['targetName'] == tf))
    listQuery = df[df['targetName'] == tf]['queryName'].tolist()
    #print(listQuery)
    #out = out.append([tf, num, listQuery])
    output.append([tf, num, ','.join(listQuery)])

out = pd.DataFrame(columns=['TF family name', 'number of TFs','list of TFs'], data=output)
out.to_csv(f"./{args.f}.summary.csv", sep=',', header=True)



