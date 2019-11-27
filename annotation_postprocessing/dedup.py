# retain only 1 protein fasta entry for each gene (remove isoforms)

import sys

fasta = sys.argv[1]
protDict = {}

with open(fasta) as f:
    line = f.readline()
    currProt = ''
    while line:
        if line.startswith('>'):
            proteinID = line.strip()[1:]
            currProt = proteinID
            protDict[currProt] = ''
        else:
            protDict[currProt] += (line.strip()+'\n')
        line = f.readline()

with open(f'{fasta}.dedup','w') as output:
    for key,value in protDict.items():
        value = value.strip()
        if not 'isoform' in key:
            output.write(f'>{key}\n')
            output.write(f'{value}\n')
        elif 'isoform 1' in key or 'isoform X1' in key:
            output.write(f'>{key}\n')
            output.write(f'{value}\n')
        else:
            continue
