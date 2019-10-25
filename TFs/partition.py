'''
Given a protein fasta file, produce, for each transcription family, a separate fasta file containing sequence
of proteins belonging to each family
'''

import argparse
import pandas as pd
parser = argparse.ArgumentParser(
    description='produce separate fasta files for protein sequences belonging to different TF families')
parser.add_argument('-f', action="store", dest="f", type=str, required=True,
                    help='path to the protein fasta file')
parser.add_argument('-t',action="store", dest='t',type=str, required=True,
                    help='path to the .csv file processed from hmmer output by countTF.py')
parser.add_argument('--prefix',action="store",dest='p',type=str,required=True,
                    help='path to the directory to which output files should be written')
args = parser.parse_args()

#first produce a dictionary of from ID to protein fasta sequence
seqDict = {}
with open(args.f) as fasta:
    line = fasta.readline()
    while line:
        if line.startswith('>'):
            curr = line.strip()[1:]
            seqDict[curr] = ''
        else:
            seqDict[curr] += line.strip()
        line = fasta.readline()

#print(args.t)
df = pd.read_csv(args.t, sep=',', header=0)
for tf in df["TF family name"]:
    geneList = df[df['TF family name'] == tf]['list of TFs'].tolist()[0].split(',')
    #print(f"geneList is {geneList}")
    with open(f'{args.p}/{tf}.fasta','w') as out:
        for gene in geneList:
            out.write(f'>{gene}\n')
            out.write(f'{seqDict[gene]}\n')

