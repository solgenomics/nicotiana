# return a list of homeologous gene pairs that are not biased in all tissues given.
import sys
from os import walk
import pandas as pd
from functools import reduce
from collections import defaultdict

dict = defaultdict(lambda:[])
fileDir = sys.argv[1]
for (dirpaths, dirnames, filenames) in walk(fileDir):
    for file in filenames:
        print(file)
        df = pd.read_csv(f"./{fileDir}/{file}", sep=",", header=0)
        for s,t,q in zip(df.Nsyl, df.Ntom, df.q_value):
            dict[(s,t)].append(q <= 0.05)

    #for file in [file for file in filenames if 'Nsyl' in file]:
    #    df_Nsyl = pd.read_csv(f"./{fileDir}/{file}", sep=",", header=0)
    #    file_sister = file.replace("Nsyl","Ntom")
    #    df_Ntom = pd.read_csv(f"./{fileDir}/{file_sister}",sep=",",header=0)
    #    # subset the dataframe such that only non-sig gene pairs are retained
    #    df_Nsyl = df_Nsyl[df_Nsyl.q_value > 0.05]
    #    df_Ntom = df_Ntom[df_Ntom.q_value > 0.05]

    #    non_sig_Nsyl = df_Nsyl.Nsyl + df_Ntom.Nsyl
    #    non_sig_Ntom = df_Nsyl.Ntom + df_Ntom.Ntom
    #    for i in range(len(non_sig_Nsyl)):
    #        dict[(non_sig_Nsyl[i], non_sig.Ntom[i])].append()

    #    nonSig.append(frozenset(map(tuple, zip(df_Nsyl.Nsyl + df_Ntom.Nsyl, df_Nsyl.Ntom + df_Ntom.Ntom))))

#print(len(dict))
with open('80%.S-biased.txt','w') as output:
    for pair,indicator in dict.items():
        #print(indicator)
        if indicator.count(True)/len(indicator) >= 0.8:
            output.write(f"{pair[0]}\t{pair[1]}\n")


#print(len(nonSig))
#nonSigConsensus = reduce((lambda set1, set2 : set1 & set2), nonSig)
#print(len(nonSigConsensus))
#with open('non.sig.consensus.txt','w') as output:
#    for pair in nonSigConsensus:
#        output.write(f'{pair[0]}\t{pair[1]}\n')

