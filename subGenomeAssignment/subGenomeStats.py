# calculate several basic stats in terms of subgenome
import sys
import bisect
import numpy as np
from functools import reduce

def calcComposition(file):
    Nsyl, Ntom, Unassigned, Ambiguous = 0,0,0,0
    with open(file) as f:
        line = f.readline()
        while line:
            _, start, end, label = line.strip().split('\t')
            span = int(end) - int(start)
            if (span <= 0):
                print(f"span smaller than 1 with start={start} and end={end} ")

            if label == 'Nsyl':
                Nsyl += span
            elif label == 'Ntom':
                Ntom += span
            elif label == 'unassigned':
                Unassigned += span
            else:
                Ambiguous += span
            line = f.readline()
    total = Nsyl + Ntom + Unassigned + Ambiguous
    print("Nsyl: {:.4f}".format(Nsyl/total))
    print("Ntom: {:.4f}".format(Ntom/total))
    print("Unassigned: {:4f}".format(Unassigned/total))
    print("Ambiguous: {:4f}".format(Ambiguous/total))


def find_gt(a, x):
    #Find the index of leftmost value greater than x
    i = bisect.bisect_right(a, x)
    if i != len(a):
        return i
    else:
        return len(a)-1

def find_lt(a, x):
    #Find the index of rightmost value less than x
    i = bisect.bisect_left(a, x)
    if i:
        return i-1
    else:
        return 0
    #raise ValueError



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
        return accu

    def findFirst(self, qStart, w):
        return find_gt(self.TERB, max(0, qStart - w))

    def findLast(self, qEnd, w):
        return find_lt(self.TELB, min(self.size, qEnd + w))

   


class interval(object):
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.TELength = 0
    
    def setTELength(self, L):
        self.TELength = L

class repeatCal(object):
    def __init__(self, w, cSize, TEFile, bedFile):
        self.w = w
        self.chromDict = {}
        self.blockList = []

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
                #if chrom not in ['Nitab0S','Nitab0T','Nitab0U']:
                self.addTE2Chrom(chrom, int(start), int(end))
                line = TE.readline()

        with open(bedFile) as bed:
            line = bed.readline()
            while line:
                chrom, start, end, label = line.strip().split('\t')
                #if chrom not in ['Nitab0S', 'Nitab0T' ,'Nitab0U']:
                stretch = interval(chrom, int(start), int(end))
                TELength = self.TEDensity(stretch)
                stretch.setTELength(TELength)
                self.blockList.append(stretch)
                line = bed.readline()

        total = 0
        TEtotal = 0
        for b in self.blockList:
            total += b.end - b.start
            TEtotal += b.TELength
        #total = reduce((lambda b1,b2: b1.end - b1.start + b2.end - b2.start), self.blockList)
        #TEtotal = reduce((lambda b1,b2: b1.TELength + b2.TELength), self.blockList)
        print(f'Total Length: {total}\nTotal TE Length: {TEtotal}')
        print("percent of TE:{:.4f}".format(TEtotal/total))

    def addChrom(self, chromID, size):
        self.chromDict[chromID] = chrom(chromID, size)

    def addTE2Chrom(self, chromID, start, end):
        if not chromID in self.chromDict:
            return
        self.chromDict[chromID].addTE(TE(chromID, start, end))

    def TEDensity(self, interval):
        return self.chromDict[interval.chrom].query(interval.start, interval.end, self.w)



# calculate 
#calcComposition(sys.argv[1])
obj = repeatCal(0, sys.argv[1], sys.argv[2], sys.argv[3])