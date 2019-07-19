# polish subgenome assignment for each individual chromosome

from optparse import OptionParser
import subprocess
import os
from collections import deque
import copy

frac = 0.01
parser = OptionParser()
parser.add_option("-g","--gap",action="store",type="string",dest="gap",
                  help="a bed file describing the gaps in the genome assembly. All gap regions will be automatically labelled as unassigned.")
parser.add_option("-f","--fileIn",action="store",type="string",dest="assign",
                  help="A subgenome assignment file in block form.")
parser.add_option("-c","--chrom",action="store",type="string",dest="chromName",
                  help="Chromosome Name")
(options,args) = parser.parse_args()
chrom = options.chromName

def collapse(segList):
    segQ = deque(segList)
    new_segQ = deque()

    if len(segQ) == 0: return
    prev = segQ.popleft()
    while len(segQ) > 2:
        curr = segQ.popleft()
        next = segQ.popleft()
        if prev['a'] != next['a']:
            new_segQ.append(prev)
            prev = curr
            segQ.appendleft(next)
        else:
            if (curr['e']-curr['s'])/(next['e']-prev['s']) < frac:
                new_segQ.append({'s':prev['s'],'e':next['e'],'a':prev['a']})
                prev = segQ.popleft()
            else:
                new_segQ.append(prev)
                prev = curr
                segQ.appendleft(next)

    while len(segQ) > 0:
        new_segQ.append(segQ.popleft())
    return new_segQ


def gapFill(d):
    # This current implemenation of gap filling is too buggy; let's try a simpler approach
    # use bedtools intersect -a <> -b <> -wa to store segments that has overlap with gaps
    # and also use bedtools intersect -a <> -b <> to get the position of each overlap
    # then we read the segment file, when encounter segments that's has overlap, do something;
    # otherwise do nothing.

    # the input deque d is the representation of subgenome assignment after collapsing
    # calculate the overlap between each assignment block and the gap region, using -wa
    new_d = deque()
    subprocess.call('bedtools intersect -a {}.polished.temp -b {} -wa > inter.gap.{}.wa.temp'.
                    format(chrom,options.gap,chrom), 
                    shell=True)
    subprocess.call('bedtools intersect -a {}.polished.temp -b {} > inter.gap.{}.temp'.
                    format(chrom,options.gap,chrom),
                    shell= True)
    dict = {}
    interFile = open('inter.gap.{}.wa.temp'.format(chrom),'r')
    interLine = interFile.readline()
    while interLine:
        content_inter = interLine.strip().split('\t')
        start = int(content_inter[1])
        if start in dict:
            dict[start] += 1
        else:
            dict[start] = 1
        interLine = interFile.readline()
    interFile.close()

    gapFile = open('inter.gap.{}.temp'.format(chrom),'r')
    gapLine = None
    while len(d) > 0:
        curr = d.popleft()
        if curr['s'] in dict and curr['a'] != 'unassigned':
            #print("case one")
            # need to split it up
            num = dict[curr['s']]
            # a non-unassigned block may have more than one overlaps with gaps
            # assume gap_start >= curr['s']
            curr_pos = curr['s']
            while num > 0:
                gapLine = gapFile.readline()
                content_gap = gapLine.strip().split('\t')
                gap_start,gap_end = int(content_gap[1]),int(content_gap[2])
                if gap_start > curr_pos:
                    new_d.append({'s':curr_pos,'e':gap_start,'a':curr['a']})
                    new_d.append({'s':gap_start,'e':gap_end,'a':'unassigned'})
                    curr_pos = gap_end
                num -= 1
            if curr_pos < curr['e']:
                new_d.append({'s':curr_pos,'e':curr['e'],'a':curr['a']})
            elif curr_pos > curr['e']:
                #if we are here, this indicates that a gap spans several blocks
                while curr_pos >= curr['e']:
                    curr = d.popleft()
                if curr['a'] != 'unassigned':
                    new_d.append({'s':curr_pos,'e':curr['e'],'a':curr['a']})
                else:
                    prev = new_d.pop()
                    new_d.append({'s':prev['s'],'e':curr['e'],'a':'unassigned'})
        elif curr['s'] in dict and curr['a'] == 'unassigned':
            #print("case2")
            num = dict[curr['s']]
            new_d.append(curr)
            # consume useless gap intersections
            while num > 0:
                gapLine = gapFile.readline()
                num -= 1
            # in case the last gap invades into a next block
            # note that the next block will be guaranteed to be non-unassigned by the way we define block
            content_gap = gapLine.strip().split('\t')
            gap_end = int(content_gap[2])
            if gap_end > curr['e']:
                next = d.popleft()
                new_d.append({'s':gap_end,'e':next['e'],'a':next['a']})
        else:
            #print("case3")
            new_d.append(curr)

    return new_d


def writeAssign(d):
    out = open('{}.polished.temp'.format(chrom),'w')
    print(len(d))
    while len(d) > 0:
        curr = d.popleft()
        out.write('{}\t{}\t{}\t{}\n'.format(chrom,curr['s'],curr['e'],curr['a']))
    out.close()


segment = [] # this will store a dictionary for each assignment
assignFile = open(options.assign,'r')
assignLine = assignFile.readline()
while assignLine:
    content = assignLine.strip().split('\t')
    start = int(content[1])
    end = int(content[2])
    assign = content[3]
    segment.append({'s':start, 'e':end, 'a':assign})
    assignLine = assignFile.readline()
assignFile.close()

d = collapse(segment)
d_copy = copy.deepcopy(d)
writeAssign(d)
d_new = gapFill(d_copy)
writeAssign(d_new)
