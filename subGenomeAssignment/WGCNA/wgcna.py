import pandas as pd
import numpy as np
from os import walk
from optparse import OptionParser
import sys
sys.path.append("../LRT4HEBS_v1.1")
import gffOp
#import rpy2.robjects as ro


parser = OptionParser()

parser.add_option('--htseq',action="store",type="string",dest="htseq_dir",
                  help="Path to a directory that contains htseq count output for all biological replicates of one experimental condition")
parser.add_option('--gff',action="store",type="string",dest="gff",
                  help="path to gff file")
(options,args) = parser.parse_args()

def TPM(htseq_dir, gff):
    gc_dict = gffOp.geneCount(htseq_dir)
    exonDict = gffOp.calcExonLen(gff, [])
    gcList = []
    exonList = []
    gList = []
    for g,c in gc_dict.items():
        gcList.append(c)
        exonList.append(exonDict[g])
        gList.append(g)
    # now every thing is in the same order, we can do vector calc below
    # transpose gcList so that, each row is the count for every gene in that condition

    # a matrix of raw counts
    gcMatrix = np.array(gcList)
    gcList = np.array(gcList).T
    T = np.sum(gcList/np.array(exonList), axis=1)[:,np.newaxis] # sum by row, that gives the scaling T for each condition
    # Now calculate TPM for each gene in each biological replicate
    # +1 as pseudocounts to avoid calculatioin as log 0
    # use log and exp operation to avoid overflow
    # formula taken from Wagner et al.
    TPM = np.exp(np.log(gcList+1)+np.log(1e6)-np.log(exonList)-np.log(T))
    TPM = TPM.T
    condition = [filename[:filename.find('.htseq')] 
                 for dirpath, dirnames, filenames in walk(htseq_dir)
                 for filename in sorted(filenames)]
    pd_TPM = pd.DataFrame(TPM, columns=condition, index=gList)
    pd_raw = pd.DataFrame(gcMatrix, columns=condition, index=gList)
    # filter out lowly expressed genes since they are probably just noises
    pd_raw = pd_raw[pd_TPM.mean(axis=1)>1]
    pd_raw.to_csv('rawCounts.csv')

TPM(options.htseq_dir, options.gff)