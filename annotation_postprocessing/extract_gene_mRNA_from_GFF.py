#!/usr/bin/env python

#for Nitab maker genes - extracting genes using maker ID with gffutils
#This script will filter out only gene and mRNA based on their ID from input file. 
#Mulitple isoforms will be discarded if it is not in the ID file.
import gffutils
import numpy as np
import argparse
from sys import stdout, stderr
from pathlib import Path

parser = argparse.ArgumentParser(
    description='Script to filter genes and their isoform from input annotation GFF file',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("input", help="Path of input GFF3 file.", type=Path)
parser.add_argument(
    "-o", "--output",
    help="Path of output GFF3 file. Prints to STDOUT if not provided.",
    type=Path, default=None, required=False)
parser.add_argument(
	"-f", "--filter_list",
	help="Path of file with tab delimited list of gene and transcript ID's",
	type=Path, default=None, required=True)
args = parser.parse_args()

assert args.input.exists()

db_save = "%s.db" % (args.input)
config = Path(db_save)

if config.is_file():
	print("gffutils (%s) database already created" % db_save)
else:
	print("------ Saving file %s ------\n" % db_save)
	#create gffutils db - save db to file as it takes long time depending on size of GFF
	fn = gffutils.example_filename(args.input)
	#create gff db - may need to change options based on the input gff format. see ID_spec from gffutils http://daler.github.io/gffutils/database-ids.html
	#create_db takes long time based on the size of input GFF. Therefore saving it to file with extension .db for future use. comment it out once you have database created for future runs.
	gffutils.create_db(fn, db_save, merge_strategy="merge")


#import db for GFF using FeatureDB
db = gffutils.FeatureDB(db_save)

#import tab delimited file containing gene\tmRNA ID's
#e.g. 
#grep -P "\tmRNA\t" annotation.gff | cut -f 9 | awk -F ";" '{print $2,$1}' | sed -e 's/Parent=//' -e 's/ ID=/\t/' > gene-mRNA.ID
#pick first (which is also longest from Mikado isoforms) transcript from gff with mRNA ID ending in mRNA-1
#sort -k1,1n -k 2,2n gene-mRNA.ID | awk '{if(!seen[$1]++) print $0}' > longest_isoform.gene-mRNA.ID

#Use numpy array to import gene and mRNA ID's from tab delimited text file
gene_isoform = np.loadtxt(args.filter_list, dtype="str", delimiter="\t")
mRNA_count = len(gene_isoform)
print("%s transcripts will be filtered \n" % mRNA_count)
#Output only genes and mRNA from provided list

with open(args.output, 'w') as f:
    i = 0
    while i < len(gene_isoform):
        print(db[gene_isoform[i, 0]], file=f)
        print(db[gene_isoform[i, 1]], file=f)
        for x in db.children(gene_isoform[i, 1]):
            print(x, file=f)
        i += 1



# #create gffutils db - save db to file as it takes long time depending on size of GFF
# fn = gffutils.example_filename('/Users/prashant/Documents/Nicotiana/Nitab_v5/Annotation/maker_run/x4-nitab5-v1.all.maker.fixed.gff')
# #create gff db - may need to change options based on the input gff format. see ID_spec from gffutils http://daler.github.io/gffutils/database-ids.html
# gffutils.create_db(fn, '/Users/prashant/Documents/Nicotiana/Nitab_v5/Annotation/maker_run/x4-nitab5-v1.all.maker.fixed.db', merge_strategy="merge")

# #import db for GFF using FeatureDB
# db = gffutils.FeatureDB('/Users/prashant/Documents/Nicotiana/Nitab_v5/Annotation/maker_run/x4-nitab5-v1.all.maker.fixed.db')

# #import tab delimited file containing gene\tmRNA ID's
# #e.g. from 
# #grep -P "\tmRNA\t" annotation.gff | cut -f 9 | awk -F ";" '{print $2,$1}' | sed -e 's/Parent=//' -e 's/ ID=/\t/' > gene-mRNA.ID
# #pick first (which is also longest from MIkado) transcript from gff with mRNA ID ending in mRNA-1
# #sort -k1,1n -k 2,2n gene-mRNA.ID | awk '{if(!seen[$1]++) print $0}' > longest_isoform.gene-mRNA.ID

# #Use numpy array to import in gene and mRNA ID's
# gene_isoform = np.loadtxt('/Users/prashant/Documents/Nicotiana/Nitab_v5/Annotation/maker_run/test.gene-mRNA.ID', dtype="str", delimiter="\t")

# #Output only genes and mRNA from provided list

# with open('/Users/prashant/Documents/Nicotiana/Nitab_v5/Annotation/maker_run/test10.gff', 'w') as f:
#     i = 0
#     while i < len(gene_isoform):
#         print(db[gene_isoform[i, 0]], file=f)
#         print(db[gene_isoform[i, 1]], file=f)
#         for x in db.children(gene_isoform[i, 1]):
#             print(x, file=f)
#         i += 1


