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
(options,args) = parser.parse_args()

# We can query the TE density from this object with gene ID
obj = repeatCal(options.winSize, options.chromSize, options.repeats, options.gff)

