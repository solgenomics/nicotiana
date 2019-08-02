import matplotlib.pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles, venn2
from os import walk
from operator import attrgetter
import collections

# This module contains several auxiliary functions used in homeologous gene expression bias analysis

# Given a gff file, calculate exonic length of each gene.
# Return a dictionary. <gene, exonic length>
# You can also specify a list of chromosomes such that genes on these chromosomes will be ignored
def calcExonLen(gff, ignore_list):
    exonLen = {}
    with open(gff) as gff:
        line = gff.readline()
        while line:
            if line.startswith('##sequence-region') or line.startswith('##gff-version'):
                line = gff.readline()
                continue
            if line.startswith('###'):
                gene_line = gff.readline().strip().split('\t')
                chrom = gene_line[0]
                # the last line (which is preceded by a line of '###') of the gff file is empty, so we check chrom == '' here
                if not chrom == '' and not chrom in ignore_list:
                    pos_1 = gene_line[8].find('Name=') + 5
                    pos_2 = gene_line[8][pos_1:].find(';')
                    name = gene_line[8][pos_1:pos_2+pos_1] if pos_2 != -1 else gene_line[8][pos_1:]
                    
                    line = gff.readline()
                    exon = 0
                    while line and (not line.startswith('###')):
                        feature_line = line.strip().split('\t')
                        if feature_line[2] == 'exon':
                            exon += (int(feature_line[4])-int(feature_line[3]))
                        line = gff.readline()
                    exonLen[name] = exon
                else:
                    line = gff.readline()
                    while line and not line.startswith('###'):
                        line = gff.readline()
            else:
                line = gff.readline()
    return exonLen


# Return a dictionary with <g1, g2> such that g1 comes from one subgenome and g2 coming from the other, and they
# are homeologous genes.
def getHomeoGenePair(file):
    dict = collections.OrderedDict()
    with open(file) as f:
        line = f.readline()
        while line:
            genes = line.strip().split('\t')
            dict[genes[0]] = genes[1]
            line = f.readline()
    return dict


# Given a list of files of htseq-count output, return a dictionary <g, count_list> such that count_list is
# a list of read counts in each replicate
# such dictionary will be used in construting the matrix that will be passed into get_R()
# and also used as part of the input that feeds into LRT_NB_HEB
def geneCount(fileDir):
    # Fetch all files in the given directory
    f = []
    for (dirpath, dirnames, filenames) in walk(fileDir):
        f.extend(filenames)
        break

    countDict = {}
    first = True
    for file in f:
        with open(f'{fileDir}/{file}') as htseq_count:
            line = htseq_count.readline()
            while line:
                if not line.startswith('__'): # the last few lines that start with __ in htseq-count output are summary stats
                    content = line.strip().split('\t')
                    gene, count = content[0], int(content[1])
                    if first:
                        countDict[gene] = [count]
                    else:
                        countDict[gene].append(count)
                line = htseq_count.readline()
        first = False
    return countDict

def printSet(inSet, dir):
    with open(f'Intersection of genes towards {dir}.txt','w') as out:
        (out.write(f'{gPair.h1}\t{gPair.h2}\n') for gPair in inSet)


def intersection(experiments):
    Nsylset = []
    Ntomset = [] # a list of sets
    Nsylset = (set(gPair for gPair in expr.gPairList if gPair.isSig and gPair.HEB > 0)
                for expr in experiments)
    Ntomset = (set(gPair for gPair in expr.gPairList if gPair.isSig and gPair.HEB < 0)
                   for expr in experiments)
    consistent_Nsyl = Nsylset[0].intersection((Nsylset[i] for i in range(1,len(Nsylset))))
    consistent_Ntom = Ntomset[0].intersection((Ntomset[i] for i in range(1,len(Ntomset))))
    printSet(consistent_Nsyl,'Nsyl')
    printSet(consistent_Ntom,'Ntom')


def draw_venn3(experiments):
    if len(experiments)==3:
        set1, set2, set3 = (set(gPair for gPair in expr.gPairList if gPair.isSig)
                   for expr in experiments)
        expr1, expr2, expr3 = (attrgetter('name')(expr) for expr in experiments)
        fig = plt.figure()
        v = venn3([set1, set2, set3], set_labels=(expr1, expr2, expr3))
        plt.title(f'Venn Diagram of HEB genes in {expr1}, {expr2}, {expr3}')
        plt.savefig(f'Venn Diagram of HEB genes in {expr1}, {expr2}, {expr3}.jpg',dpi=300, format='jpg',quality=95)
        #plt.show()

        set1, set2, set3 = (set(gPair for gPair in expr.gPairList if gPair.isSig and gPair.HEB > 0)
                   for expr in experiments)
        fig = plt.figure()
        v = venn3([set1, set2, set3], set_labels=(expr1, expr2, expr3))
        plt.title(f'Venn Diagram of HEB genes in towards Nsyl')
        plt.savefig(f'Venn Diagram of HEB genes biased towards Nsyl.jpg',dpi=300, format='jpg',quality=95)

        set1, set2, set3 = (set(gPair for gPair in expr.gPairList if gPair.isSig and gPair.HEB < 0)
                   for expr in experiments)
        fig = plt.figure()
        v = venn3([set1, set2, set3], set_labels=(expr1, expr2, expr3))
        plt.title(f'Venn Diagram of HEB genes in towards Ntom')
        plt.savefig(f'Venn Diagram of HEB genes biased towards Ntom.jpg',dpi=300, format='jpg',quality=95)

    else:
        raise ValueError("Must have three experiments to use venn3")

def draw_venn2(experiments):
    if len(experiments)==2:
        set1, set2 = (set(gPair for gPair in expr.gPairList if gPair.isSig)
                   for expr in experiments)
        expr1, expr2 = (attrgetter('name')(expr) for expr in experiments)
        fig = plt.figure()
        v = venn2([set1, set2], set_labels=(expr1, expr2))
        plt.title(f'Venn Diagram of HEB genes in {expr1}, {expr2}')
        plt.savefig(f'Venn Diagram of HEB genes in {expr1}, {expr2}.jpg', dpi=300, format='jpg',quality=95)
        #plt.show()

        set1, set2 = (set(gPair for gPair in expr.gPairList if gPair.isSig and gPair.HEB > 0)
                   for expr in experiments)
        fig = plt.figure()
        v = venn2([set1, set2], set_labels=(expr1, expr2))
        plt.title(f'Venn Diagram of HEB genes in towards Nsyl')
        plt.savefig(f'Venn Diagram of HEB genes biased towards Nsyl.jpg', dpi=300, format='jpg',quality=95)

        set1, set2 = (set(gPair for gPair in expr.gPairList if gPair.isSig and gPair.HEB < 0)
                   for expr in experiments)
        fig = plt.figure()
        v = venn2([set1, set2], set_labels=(expr1, expr2))
        plt.title(f'Venn Diagram of HEB genes in towards Ntom')
        plt.savefig(f'Venn Diagram of HEB genes biased towards Ntom.jpg', dpi=300, format='jpg',quality=95)

    else:
        raise ValueError("Must have two experiments to use venn2")


# draw scatter plot for visualizing changes in HEB under two different experiment settings
def draw_scatter(gpTwoL, p_cut, name1, name2):
    # plot scatter plot of bias for each gene pair in two different experimental conditions
    gpTwo_nonSig = [gptwo for gptwo in gpTwoL if gptwo.pValue > p_cut]
    gpTwo_Sig = [gptwo for gptwo in gpTwoL if gptwo.pValue <= p_cut]
    fig = plt.figure()
    plt.scatter([gptwo.HEB1 for gptwo in gpTwo_nonSig], [gptwo.HEB2 for gptwo in gpTwo_nonSig],
                c='#A9A9A9', label=f'$\Delta$HEB not significant(n={len(gpTwo_nonSig)})',
                s=(plt.rcParams['lines.markersize'] ** 2)/16)
    x_delta_sig_toNsyl = [gptwo.HEB1 for gptwo in gpTwo_Sig if gptwo.delta_HEB > 0]
    y_delta_sig_toNsyl = [gptwo.HEB2 for gptwo in gpTwo_Sig if gptwo.delta_HEB > 0]
    x_delta_sig_toNtom = [gptwo.HEB1 for gptwo in gpTwo_Sig if gptwo.delta_HEB < 0]
    y_delta_sig_toNtom = [gptwo.HEB2 for gptwo in gpTwo_Sig if gptwo.delta_HEB < 0]
    plt.scatter(x_delta_sig_toNsyl, y_delta_sig_toNsyl, c='#0000CD', 
                label=f'Significant $\Delta$HEB towards Nsyl($n_0$={len(x_delta_sig_toNsyl)})',
                s=(plt.rcParams['lines.markersize'] ** 2)/16)
    plt.scatter(x_delta_sig_toNtom, y_delta_sig_toNtom, c='#FF0000', 
                label=f'Significant $\Delta$HEB towards Ntom($n_1$={len(y_delta_sig_toNtom)})',
                s=(plt.rcParams['lines.markersize'] ** 2)/16)
    plt.xlabel(f'HEB in {name1}')
    plt.ylabel(f'HEB in {name2}')
    plt.title(f'$\Delta$HEB between {name1} and {name2}')
    plt.legend(loc='lower right', fontsize='small')
    plt.savefig(f'Changes of HEB between {name1} and {name2}.jpg',dpi=300,format='jpg', quality=95)

def draw_hist_changeHEB(gpTwoL, p_cut, name1, name2):
    # plot histogram of delta_bias
    delta_bias_all = [gptwo.delta_HEB for gptwo in gpTwoL]
    delta_bias_sig = [gptwo.delta_HEB for gptwo in gpTwoL if gptwo.pValue < p_cut]

    fig = plt.figure()
    bin_all = np.arange(-15,15,0.25)
    bin_pos = np.arange(0,15,0.25)
    bin_neg = np.arange(-15,0,0.25)
    delta_bias_sig_pos = [b for b in delta_bias_sig if b > 0]
    delta_bias_sig_neg = [b for b in delta_bias_sig if b < 0]
    plt.hist(delta_bias_all, bins=bin_all, color='#A9A9A9', label='Overall n={}'.format(len(delta_bias_all)))
    plt.hist(delta_bias_sig_pos, bins=bin_pos, color='#0000CD', label='Bias shift towards Nsyl, $n_0$={}'.format(len(delta_bias_sig_pos)))
    plt.hist(delta_bias_sig_neg, bins=bin_neg, color='#FF0000', label='Bias shift towards Ntom, $n_1$={}'.format(len(delta_bias_sig_neg)))
    plt.xlabel('Changes in Homeolog Expression Bias $\Delta HEB$')
    plt.title(f'$\Delta$HEB between {name1} and {name2}')
    plt.legend(loc='upper right', fontsize='small')
    plt.savefig(f'Homeologous Gene Expression Bias Changes between {name1} and {name2}.jpg', dpi=300, format='jpg', quality=95)