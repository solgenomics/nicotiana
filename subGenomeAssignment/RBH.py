# find the reciprocal best hits given two .blast files

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-b",action="store",type="string",dest="blast",
                  help="path to the two blast files, delimited by comma")
#parser.add_option("-g",action="store",type="string",dest="geneSet",
#                  help="path to the two files that storing genes from different subgenomes. FASTA format")
(options,args) = parser.parse_args()

blastFile = options.blast.strip().split(',')
blast1, blast2 = blastFile[0], blastFile[1]
#geneFile = options.geneSet.strip().split(',')
#geneF1, geneF2 = geneFile[0], geneFile[1]


def find_best_hit(file, dict):
    with open(file) as f:
        line = f.readline()
        prev_query_gene = ''
        while line:
            # the blast file is sorted by descending bit-score
            blast_info = line.strip().split('\t')
            query_gene, ref_gene, seq_iden = blast_info[0], blast_info[1], float(blast_info[2])
            if prev_query_gene != query_gene and seq_iden > 70.0:
                dict[query_gene] = ref_gene
            line = f.readline()


def build_gene_set(geneF, geneSet):
    with open(geneF) as f:
        line = f.readline()
        while line:
            if line.startswith('>'):
                geneSet.add(line.strip()[1:])
            line = f.readline()

def diff_origin(gene1, gene2, geneSet1, geneSet2):
    return (gene1 in geneSet1 and gene2 in geneSet2) or (gene1 in geneSet2 and gene2 in geneSet1)

rbh_dict1 = {}
rbh_dict2 = {}
geneSet1 = set()
geneSet2 = set()
find_best_hit(blast1, rbh_dict1)
find_best_hit(blast2, rbh_dict2)
#build_gene_set(geneF1, geneSet1)
#build_gene_set(geneF2, geneSet2)

for q,r in rbh_dict1.items():
    if r in rbh_dict2 and rbh_dict2[r] == q:       
        #if diff_origin(q, r, geneSet1, geneSet2):
            print('{}\t{}'.format(q,r))


   