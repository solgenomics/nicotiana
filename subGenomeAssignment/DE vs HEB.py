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
parser.add_option("-n","--name",action="store",type="string",dest="name",
                  help="Such as meJA_T2, etc.")
(options,args) = parser.parse_args()

def getGeneSetFromHEB(file, col):
    gl = set()
    with open(file) as f:
        line = f.readline()
        while line:
            gl.add(line.strip().split('\t')[col])
            line = f.readline()
    return gl

def getGeneSetFromCSV(file, fdr):
    gl = set()
    with open(file, newline='') as csvfile:
        r = csv.reader(csvfile, delimiter=',')
        first = True
        for row in r:
            if not first: # omit the header line
                gene,mean,LFC,lfcSE,stat,p,padj = row
                if padj != 'NA' and float(padj) < fdr:
                    gl.add(gene)
            else:
                first = False
    return gl

toNsyl,toNtom = options.HEB.split(',')
g2Nsyl = getGeneSetFromHEB(toNsyl,0)
g2Ntom = getGeneSetFromHEB(toNtom,1)
DEGene = getGeneSetFromCSV(options.csv, options.fdr)
fig = plt.figure()
v = venn3([g2Nsyl,g2Ntom,DEGene], set_labels=("towards Nsyl", "towards Ntom", "DE"))
plt.title(f'Venn Diagram of HEB and DE Genes for {options.name}')
plt.savefig(f'Venn Diagram of HEB and DE Genes for {options.name}.jpg',dpi=300, format='jpg',quality=95)




