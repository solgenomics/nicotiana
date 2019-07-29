from repeats import TE, chrom, gene, repeatCal
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import math
import gffOp

# Explore relationship between TE density around a gene and its expression level
parser = OptionParser()
parser.add_option("--gff",action="store",type="string",dest="gff",
                  help="path to the GFF file")
parser.add_option("-r","--repeats",action="store",type="string",dest="repeats",
                  help="Path to a repeat bedgraph file")
parser.add_option('-w', action="store",type="int",dest="winSize", default=5000,
                  help="Window size. Default 5kb")
parser.add_option('-s','--size',action='store',type='string',dest='chromSize',
                  help='Path to the chromosome size file')
parser.add_option("-i","--ignore",action="store",type="string",dest="garbage",
                  help="A list of scaffolds that you'd like to ignore in downstream analysis")
parser.add_option('-h','--htseq',action="store",type="string",dest="htseq_dir",
                  help="Path to a directory that contains htseq count output for all biological replicates of one experimental condition")
(options,args) = parser.parse_args()

# We can query the TE density from this object with gene ID
gDen = repeatCal(options.winSize, options.chromSize, options.repeats, options.gff).gene_Dict
gc_dict = gffOp.geneCount(options.htseq_dir)
exon_dict = gffOp.calcExonLen(options.gff, options.garbage.split(','))

# sort gene ID by increasing TE density
gList = [g for g,d in gDen.items()]
gList.sort(key=lambda g:gDen[g])
exonList = [exon_dict[g] for g in gList]
gcList = [gc_dict[g] for g in gList] # now every thing is in the same order, we can do vector calc below
# transpose gcList so that, each row is the count for every gene in that biological replicate
gcList = gcList.T
T = np.sum(gcList/exonList, axis=1) # sum by row, that gives the scaling T for each replicate
# Now calculate TPM for each gene in each biological replicate
# +1 as pseudocounts to avoid calculatioin as log 0
# use log and exp operation to avoid overflow
# formula taken from Wagner et.al
TPM = np.exp(np.log(gcList+1)+np.log(1e6)-np.log(exonList)-np.log(T))
TPM = np.mean(TPM, axis=0) # take the mean column wise
print(TPM)



