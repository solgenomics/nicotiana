# Draw histogram of size distribution of the block of interest.
# Supply two command line arguments: subgenome assignment txt file in blocks form
# and the block ID (eg. Nsyl, Ntom, ambiguous, unassigned)
# eg. python hist.py -f <path to block assignment file> -i ambiguous
import sys
import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser
import statistics as stats


parser = OptionParser()
parser.add_option("-f","--file",action="store",type="string",dest="subGenomeAssignFile",
                  help="path to the file produced by assignBase.py in blocks form")
parser.add_option("-i","--id",action="store",type="string",dest="ID",
                  help="the ID of the block whose size distribution you want to draw, eg. Nsyl, Ntom, unassigned, ambiguous")
(options,args) = parser.parse_args()

file = open(options.subGenomeAssignFile,'r')
line = file.readline()
id = options.ID
size = []
while line:
    content = line.strip().split('\t')
    if id == content[3]:
        size.append(int(content[2])-int(content[1]))
    line = file.readline()

#print("max size is {}".format(max(size)))
mean = stats.mean(size)
sd = stats.stdev(size)
size = [num for num in size if num >= mean - 3*sd and num <= mean + 3*sd]
plt.hist(size,bins='auto')
#plt.ylim(0,75000)
plt.savefig("Histogram of size of {} blocks".format(id))



