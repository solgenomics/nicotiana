## Draw venn diagram, showing the overlap between DE and HEB.
## For the file lists of HEB, we assume that the first file contains genes biased towards Nsyl, and the second towards Ntom

from optparse import OptionParser
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import csv

parser = OptionParser()
parser.add_option("--HEB",action="store",type="string",dest="HEB",
                  help="Path to two files containing HEB genes. Delimited by comma.")
parser.add_option("--fdr",action="store",type="float",dest="fdr",default=0.05,
                  help="FDR")
parser.add_option("--csv",action="store",type="string",dest="csv",
                  help="Path to the csv file as produced by DESeq2")
parser.add_option("-n","--name",action="store",type="string",dest="experimental condition",
                  help="Such as meJA_T2, etc.")
(options,args) = parser.parse_args()

def getGeneSetFromHEB(file, col):
    gl = {}
    with open(file) as f:
        line = f.readline()
        while line:
            gl.append(line.strip().split('\t')[col])
            line = f.readline()
    return gl

def getGeneSetFromCSV(file, fdr):
    with open(file, newline='') as csvfile:
        r = csv.reader(csvfile, delimiter=',')
        for row in r:
            print(row)

toNsyl,toNtom = options.HEB.split(',')
g2Nsyl = getGeneSetFromHEB(toNsyl,0)
g2Ntom = getGeneSetFromHEB(toNtom,1)
DEGene = getGeneSetFromCSV(options.csv, options.fdr)
# now process the given csv file (output from DESeq2)



