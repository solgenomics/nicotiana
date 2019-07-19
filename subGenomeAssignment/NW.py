# custom impelmentation of Nadaraya-Watson Kernel Smoother(test version)
# suppose we are only dealing with a single chromosome
# This is the initial naive implementation

import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser
import time

parser = OptionParser()
parser.add_option("-f","--file",action="store",type="string",dest="file",
                  help="path to two bedgraph files, delimited by comma")
parser.add_option("-w","--window",action="store",type="int",dest="winSize",
                  help="Window size for kernal regression")
parser.add_option("-n","--name",action="store",type="string",dest="name",
                  help="Names for the two subgenoms. Delimited by comma and provided in the order of the bedgraph file given.")
(options,args) = parser.parse_args()

name1 = options.name.strip().split(',')[0]
name2 = options.name.strip().split(',')[1]
file1 = open(options.file.strip().split(',')[0],'r')
file2 = open(options.file.strip().split(',')[1],'r')

width = 100
if options.winSize != None:
    width = options.winSize

line1 = file1.readline()
line2 = file2.readline()
count1 = []
count2 = []
pos = []

start_read = time.process_time()
print("Start Reading!")
while line1 and line2:
    content1 = line1.strip().split('\t')
    content2 = line2.strip().split('\t')
    count1.append(int(content1[2]))
    count2.append(int(content2[2]))
    pos.append(int(content1[1]))
    line1 = file1.readline()
    line2 = file2.readline()

file1.close()
file2.close()
print("Reading data takes {:.3f}".format(time.process_time()-start_read))

start_est = time.process_time()
array1 = np.array(count1)
array2 = np.array(count2)
adjusted1 = []
adjusted2 = []

print("Starting estimating!")

for i in list(range(array1.size)):
    print(i)
    diff = ((np.arange(array1.size))-i)/width
    K = 0.75*(np.ones(array1.size,)-diff*diff)
    K[diff <= -1] = 0
    K[diff >= 1] = 0
    adjusted1.append(array1.dot(K)/K.sum())
    adjusted2.append(array2.dot(K)/K.sum())



print("NW takes {:.3f}".format(time.process_time()-start_est))
line1, = plt.plot(pos,count1,color='red')
line2, = plt.plot(pos,count2,color='blue')
line3, = plt.plot(pos,adjusted1,color='black')
line4, = plt.plot(pos,adjusted2,color='magenta')
plt.xlabel('Position along the Chromosome')
plt.ylabel('Read Counts')
plt.title('Distribution of DNA Reads along the Chromosome')
plt.legend((line1,line2,line3,line4),(name1,name2,
    'NW estimate of {} with bandwidth {}'.format(name1,width),
    'NW estimate of {} with bandwidth {}'.format(name2,width)),loc='upper right',fontsize='x-large')
plt.show()
