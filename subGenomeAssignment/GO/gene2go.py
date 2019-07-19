# Given a gff file and a tsv file (output from interproscan)
# output a gene to GO mapping in the format as required by topGO R package.
# gene_id<tab>GO_id1,GO_id2,...
from optparse import OptionParser
import re

parser = OptionParser()
parser.add_option("-g","--gff",action="store",type="string",dest="gff",
                  help="path to the .gff file")
parser.add_option("-t","--tsv",action="store",type="string",dest="tsv",
                  help="Path to the .tsv file")

(options,args) = parser.parse_args()

def featureDict(attribute):
    dict = {}
    attr_splitted = attribute.split(';')
    for p in attr_splitted:
        tag,value = p.split('=')
        dict[tag]=value
    return dict


mRNA2gene = {}
with open(options.gff) as gff:
    line = gff.readline()
    while line:
        if not line.startswith('#'):
            seqid, source, feature, start, end, score, strand, phase, attributes = line.strip().split('\t')
            if feature == 'mRNA':
                dict = featureDict(attributes)
                mRNA2gene[dict['Name']] = dict['Parent']
        line = gff.readline()

# Extract GO terms for each mRNA, and associate them with genes.
gene2go = {}
with open(options.tsv) as tsv:
    line = tsv.readline()
    while line:
        s = line.split('\t')
        mRNA = s[0]
        gene = mRNA2gene[mRNA]
        for item in s:
            if item.startswith('GO:'):
                if gene in gene2go:
                    gene2go[gene].update(go for go in item.split('|'))
                else:
                    gene2go[gene] = {go for go in item.split('|')}
        line = tsv.readline()

with open('Ntab.gene2go','w') as op:
    for gene, goTerms in gene2go.items():
        go_cat = ','.join(goTerms)
        op.write(f'{gene}\t{go_cat}\n')


#print(gene2go)
#print(f'number of mRNA {len(mRNA2gene)}')
#print(f'number of genes {len(set(mRNA2gene.values()))}')
#print(f'number of genes that have a GO term {len(gene2go)}')