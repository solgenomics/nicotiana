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
parser.add_option('--htseq',action="store",type="string",dest="htseq_dir",
                  help="Path to a directory that contains htseq count output for all biological replicates of one experimental condition")
parser.add_option('-b','--bin',action="store",type="float",dest="binSize",default=0.1,
                  help="bin sizes for grouping genes with similar TE density")

(options,args) = parser.parse_args()

# We can query the TE density from this object with gene ID
# repeatCal alreadys ignores genes from S,T,U
gDen = repeatCal(options.winSize, options.chromSize, options.repeats, options.gff).geneDict
exon_dict = gffOp.calcExonLen(options.gff, options.garbage.split(','))
gList = [g for g,d in gDen.items()]
gList.sort(key=lambda g:gDen[g])
exonList = [exon_dict[g] for g in gList]

bins = []
r = 0
index = 0
for g in gList:
    if gDen[g] > r:
        bins.append(index)
        r += options.binSize
    index += 1
bins.append(index+1)


def TPM_mean(gList, htseq_dir, exonList, bins):
    gc_dict = gffOp.geneCount(htseq_dir)
    # sort gene ID by increasing TE density
    gcList = [gc_dict[g] for g in gList] # now every thing is in the same order, we can do vector calc below
    # transpose gcList so that, each row is the count for every gene in that biological replicate
    gcList = np.array(gcList).T
    T = np.sum(gcList/exonList, axis=1)[:,np.newaxis] # sum by row, that gives the scaling T for each replicate
    #print(T)
    # Now calculate TPM for each gene in each biological replicate
    # +1 as pseudocounts to avoid calculatioin as log 0
    # use log and exp operation to avoid overflow
    # formula taken from Wagner et al.
    #print((gcList+1)[1:100])
    #print(exonList[1:100])
    TPM = np.exp(np.log(gcList+1)+np.log(1e6)-np.log(exonList)-np.log(T))
    #print(np.mean(TPM,axis=1))
    TPM = np.mean(TPM, axis=0) # take the mean column wise
    
    TPM_mean = []
    for i in range(len(bins)-1):
        TPM_mean.append(np.mean(TPM[bins[i]:bins[i+1]]))
    return TPM_mean
#print(TPM_mean)

dirs = options.htseq_dir.split(',')
TPM_mean_Root = TPM_mean(gList,dirs[0],exonList,bins)
TPM_mean_Shoot = TPM_mean(gList,dirs[1],exonList,bins)
TPM_mean_ShootApex = TPM_mean(gList,dirs[2],exonList,bins)

fig = plt.figure()
plt.plot(np.arange(0,1,options.binSize), np.log(TPM_mean_Root),'b', label='root')
plt.plot(np.arange(0,1,options.binSize), np.log(TPM_mean_Shoot),'r', label='shoot')
plt.plot(np.arange(0,1,options.binSize), np.log(TPM_mean_ShootApex),'m', label='shoot-apex')
plt.xlabel("TE Density")
plt.ylabel("mean TPM($\log_2$)")
plt.title("TE Density and Gene Expression Level")
plt.legend(loc='lower left', fontsize='large')
plt.savefig('TE.Gene.Expr.jpg', dpi=300, format='jpg', quality=95)





