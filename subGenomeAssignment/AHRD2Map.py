### Convert the output from AHRD to the gene2GO mapping as required by topGO
import sys

with open(sys.argv[1]) as ahrd:
    with open('gene2go.ahrd','w') as out:
        line = ahrd.readline()
        while line:
            if not line.startswith('#'):
                temp = line.strip().split('\t')
                if len(temp) == 3:
                    protID, _, go = temp
                    go = go.split(' ,')
                    go = ",".join(go)
                    out.write(f'{protID[:-2]}\t{go}\n')

            line = ahrd.readline()
