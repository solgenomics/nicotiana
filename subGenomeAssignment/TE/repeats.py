import bisect
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

def find_gt(a, x):
    #Find the index of leftmost value greater than x
    i = bisect.bisect_right(a, x)
    if i != len(a):
        return i
    raise ValueError

def find_lt(a, x):
    #Find the index of rightmost value less than x
    i = bisect.bisect_left(a, x)
    if i:
        return i-1
    raise ValueError



class TE(object):
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end

class chrom(object):
    def __init__(self, chrom, length):
        self.id = chrom
        self.size = length
        self.TELB = []
        self.TERB = []

    def addTE(self, TE):
        self.TELB.append(TE.start)
        self.TERB.append(TE.end)

    def query(self, qStart, qEnd, w):
        # Given a query region (specified by qStart and qEnd)
        # report the percentage of TEs w bps upstream and downstream of the given query region
        l,r = self.findFirst(qStart, w), self.findLast(qEnd, w)
        #print(f"l={l} and r={r}")
        accu = 0
        for i in range(l,r+1):
            t_start, t_end = self.TELB[i], self.TERB[i]
            lb = max(0, qStart-w, t_start)
            rb = min(self.size, qEnd+w, t_end)
            accu += (rb - lb)

        region = min(qEnd+w, self.size) - max(0, qStart - w)
        return accu/region

    def findFirst(self, qStart, w):
        return find_gt(self.TERB, max(0, qStart - w))

    def findLast(self, qEnd, w):
        return find_lt(self.TELB, min(self.size, qEnd + w))

   


class gene(object):
    def __init__(self, chrom, start, end, ID):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.ID = ID

class repeatCal(object):
    def __init__(self, w, cSize, TEFile, gffFile):
        self.w = w
        self.chromDict = {}
        self.geneDict = {}

        with open(cSize) as cs:
            line = cs.readline()
            while line:
                chrom, size = line.strip().split('\t')
                self.addChrom(chrom, int(size))
                line = cs.readline()

        # parse the TE.bed file
        with open(TEFile) as TE:
            line = TE.readline()
            while line:
                chrom, start, end = line.strip().split('\t')
                if chrom not in ['Nitab0S','Nitab0T','Nitab0U']:
                    self.addTE2Chrom(chrom, int(start), int(end))
                line = TE.readline()

        with open(gffFile) as gff:
            line = gff.readline()
            while line:
                if not line.startswith('#'):
                    seqid, source, feature, start, end, score, strand, phase, attributes = line.strip().split('\t')
                    if feature == 'gene' and seqid not in ['Nitab0S','Nitab0T','Nitab0U']:
                        gID = attributes.split(';')[0][3:]
                        g = gene(seqid, int(start), int(end), gID)
                        self.geneDict[gID] = self.TEDensity(g)
                line = gff.readline()



    def addChrom(self, chromID, size):
        self.chromDict[chromID] = chrom(chromID, size)

    def addTE2Chrom(self, chromID, start, end):
        if not chromID in self.chromDict:
            return
        self.chromDict[chromID].addTE(TE(chromID, start, end))

    def TEDensity(self, gene):
        # return the TE density around the given gene object
        if not gene.chrom in self.chromDict:
            return
        return self.chromDict[gene.chrom].query(gene.start, gene.end, self.w)

    def drawCDF(self, list1, list2, name1, name2):
        genDen_1 = [self.geneDict[g] for g in list1]
        genDen_2 = [self.geneDict[g] for g in list2]
        stat, p_value = stats.ks_2samp(genDen_1, genDen_2)
        print(stat)
        print(p_value)

        fig = plt.figure()
        bin = np.arange(0,1,0.001)
        plt.hist(genDen_1, bins=bin, color='#0000CD', label=f'{name1}',alpha=0.5, density=True, cumulative=True, histtype='step')
        plt.hist(genDen_2, bins=bin, color='#FF0000', label=f'{name2}',alpha=0.5, density=True, cumulative=True, histtype='step')
        plt.xlabel('TE Density')
        plt.title('TE Density Surrounding Genic Region')
        plt.legend(loc='upper right', fontsize='small') 
        plt.text(0.0,0.7,'K-S test={:.3f}(p={:.2E})'.format(stat, p_value))
        plt.savefig(f'TEDensity.hist.cdf.{name1} vs {name2}.jpg', dpi=300, format='jpg', quality=95)

    def drawCDFByChrom(self, chrom1, chrom2):
        genList1 = [g for g,d in self.geneDict.items() if chrom1 in g]
        genList2 = [g for g,d in self.geneDict.items() if chrom2 in g]
        self.drawCDF(genList1, genList2, chrom1, chrom2)
