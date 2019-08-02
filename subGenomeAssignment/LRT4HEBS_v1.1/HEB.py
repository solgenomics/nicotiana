# Homeologous Gene Expression Bias Analysis
import gffOp
import matlab.engine
from optparse import OptionParser
import numpy as np
from homeolog import genePair, experiment, compExpr
import itertools

parser = OptionParser()
parser.add_option("-f","--file",action="store",type="string",dest="gff",
                  help="path to the GFF file")
parser.add_option("-i","--ignore",action="store",type="string",dest="garbage",
                  help="A list of scaffolds that you'd like to ignore in downstream analysis")
parser.add_option("-r","--rna",action="store",type="string",dest="RNA",
                  help="path to htseq-count output of RNA-seq alignment. Replicates should be delimited by comma.\
                  You can provide all htseq outputs from different experiment settings, but make sure specify the\
                  number of replicates each experiment has.")
parser.add_option('--name', action="store",type="string",dest="names",
                  help="Name for each experiment")
# In our file, each line is a homeologous gene pair. 
# In each line, the first entry is the gene coming from Nsyl, and the second entry is the gene from Ntom 
parser.add_option("-p","--pair",action="store",type="string",dest="homo_pair",
                  help="path to the file specifying homeologous gene pairs.")

(options,args) = parser.parse_args()

ignore_list = options.garbage.split(',')
htseq_count_list = options.RNA.split(',')
# you may want to check if we have recorded the exon length for every gene in the dictionary
# also note that htseq contains some genes that are on chromosomes U,S,T and we don't include them here
exonLen = gffOp.calcExonLen(options.gff, ignore_list)
genePairDict = gffOp.getHomeoGenePair(options.homo_pair)

# Next we want to perform a likelihood ratio test for each homeologous gene pair by calling the matlab functions
# Before we get started, let's first list what information we need for each gene pairs (g1, g2)
# read counts for g1, g2 in different replicates
# length of exonic region of g1, g2
# a row vector of aggregation parameter for each replicate
# a row vector of the total number of mapped reads in each replicate

# so clearly, we need to estimate the aggregation parameter for each replicate first
# aggregation parameter estimate function: get_R(data,plot)
# data: an array whose rows are genes and columns are total number of mapped reads for each replicate. 
# Data should include all genes even if some of them will not be tested for homeologous expression bias.

eng = matlab.engine.start_matlab()
htseq_Dirs = (htseq_dir for htseq_dir in options.RNA.split(','))
names = [name for name in options.names.split(',')]
experiments = []
print(f'options.names{options.names}')
print("name list is {}".format(names))

for expr_name in names:
    htseq_dir = next(htseq_Dirs)
    gc_dict = gffOp.geneCount(htseq_dir)
    gclist = [count for gene,count in gc_dict.items()]
    r, D = eng.get_R(matlab.int64(gclist), nargout=2)
    print("Aggregation parameter for each replicate is {}".format(r))
    print("Total sequencing depth for each replicate is {}".format(D))
    # D is the number of total mapped reads, in millions
    # which is used in calculting FPKM in later steps
    r,D = r[0],D[0]
    expr = experiment(r, D, genePairDict, gc_dict, exonLen, expr_name)
    expr.fdr(0.05)
    experiments.append(expr)
    

# Draw Venn Diagram
if len(experiments) == 3:
    gffOp.draw_venn3(experiments)

if len(experiments) == 2:
    gffOp.draw_venn2(experiments)
    #compExpr(experiments[0], experiments[1])

if len(experiments) >= 2:
    # print the intersection of biased genes in all given experiments
    gffOp.intersection(experiments)