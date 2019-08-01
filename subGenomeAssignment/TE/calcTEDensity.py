# required input:
# a gff file such we know the coordinate of each gene
# a bedgraph file for repeats
# a list of genes belonging to subgenome A
# a list of genes belonging to subgenome B
# a chromosome size file

### First, parse the chromosome size file, build chrom objects 
### Second, parse the repeats bedgraph file, add TEs to the corresponding chromosomes
### For each gene in the gff file, calculate their TE density, store the result in a dictionary
### with key being the gene ID and value being the density
### For gene in the given list (eg, subgenome A/B), retrive their respective TE density and plot.

from repeats import TE, chrom, gene, repeatCal
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import math

parser = OptionParser()
parser.add_option("--gff",action="store",type="string",dest="gff",
                  help="path to the GFF file")
parser.add_option("-r","--repeats",action="store",type="string",dest="repeats",
                  help="Path to a repeat bedgraph file")
parser.add_option("-g","--gene",action="store",type="string",dest="geneList",
                  help="Path to file containing genes of interest. Multiple files should be delimited by commas")
parser.add_option('-w', action="store",type="int",dest="winSize", default=5000,
                  help="Window size. Default 5kb")
parser.add_option('-s','--size',action='store',type='string',dest='chromSize',
                  help='Path to the chromosome size file')
parser.add_option('-f',action="store_true",dest="isFasta", default=False,
                  help="If the gene list is given in fasta format, set -f")
parser.add_option('-p',action="store", type="string", dest="pair",
                  help="Path to files specifying homeologous gene pairs")
(options,args) = parser.parse_args()

obj = repeatCal(options.winSize, options.chromSize, options.repeats, options.gff)

def parseFasta(fastaFile):
    l = []
    with open(fastaFile) as f:
        line = f.readline()
        while line:
            if line.startswith('>'):
                l.append(line.strip()[1:])
            line = f.readline()
    return l


# Assume the first file is Nsyl, and the second file is Ntom
fasta1, fasta2 = options.geneList.split(',')
# list of genes belong to Nsyl and Ntom
Nsyl = parseFasta(fasta1)
Ntom = parseFasta(fasta2)
# list of TE density around each gene (unordered)
#obj.drawCDF(Nsyl, Ntom, 'Nsyl', 'Ntom')
obj.drawDotPlot(options.pair)
# For drawCDFByChrom, the first argument is chromsome identified as belonging to Nsyl
# and the second is identified as belonging to Ntom
#obj.drawCDFByChrom('Nitab01','Nitab23')
#obj.drawCDFByChrom('Nitab10','Nitab02')
#obj.drawCDFByChrom('Nitab16','Nitab12')
#obj.drawCDFByChrom('Nitab03','Nitab18')
#obj.drawCDFByChrom('Nitab06','Nitab04')
#obj.drawCDFByChrom('Nitab11','Nitab13')


#list1 = []
#list2 = []
#with open('rbh.v1.blastp') as pair:
#    line = pair.readline()
#    while line:
#        g1, g2 = line.strip().split('\t')
#        list1.append(g1)
#        list2.append(g2)
#        line = pair.readline()

#obj.drawCDF(list1, list2, 'Nsyl_homeolog', 'Ntom_homeolog')
