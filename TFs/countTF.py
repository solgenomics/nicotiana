import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(
        description='Path to the output from hmmscan formatted as parsable table')
parser.add_argument('-f', action="store", dest="f", type=str, required=True)
args = parser.parse_args()


df = pd.read_csv(args.f, header=None, comment='#', sep='\s{18}')
print(df)
df.columns = ['targetName','targetAccession','queryName','queryAccession',
              'E-value1','score1','bias1','E-value2', 'score2', 'bias2',
              'exp','reg','clu','ov','env','dom','rep','inc','targetDescription']

TFs = set(df['targetName'])

out = df.DataFrame(columns=['TF family name', 'number of TFs','list of TFs'])
for tf in TFs:
    num = np.sum(np.array(df['targetName'] == tf))
    listQuery = df[df['targetName' == tf]]['queryName']
    out = out.append([tf, num, listQuery])

out.to_csv(f"./{args.f}.summary.csv", sep=',', header=True)



