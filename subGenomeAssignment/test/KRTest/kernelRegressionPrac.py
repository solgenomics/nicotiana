# python kernalRegressionPrac.py -f <bedgraph file1>,<bedgraph file2>
# -w <window size, default 100> -n <subgenome name1>,<subgenome name2>

import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser
from sklearn.neighbors import KernelDensity
import time

parser = OptionParser()
parser.add_option("-f","--file",action="store",type="string",dest="file",
                  help="path to two bedgraph files, delimited by comma")
parser.add_option("-w","--window",action="store",type="int",dest="winSize",
                  help="Window size for kernal regression")
parser.add_option("-n","--name",action="store",type="string",dest="name",
                  help="Names for the two subgenoms. Delimited by comma and provided in the order of the bedgraph file given.")
parser.add_option("-k","--kernel",action="store",type="string",dest="kernel_toUse",
                  help="The kernel to use. Valid kernels are [‘gaussian’|’tophat’|’epanechnikov’|’exponential’|’linear’|’cosine’] Default is ‘gaussian’")
(options,args) = parser.parse_args()

kernel_toUse = 'gaussian'
if options.kernel_toUse != None:
    kernel_toUse = options.kernel_toUse

name1 = options.name.strip().split(',')[0]
name2 = options.name.strip().split(',')[1]
file1 = open(options.file.strip().split(',')[0],'r')
file2 = open(options.file.strip().split(',')[1],'r')
line1 = file1.readline()
line2 = file2.readline()
count1 = []
count2 = []
list1 = []
list2 = []
pos = []

start_read = time.process_time()

while line1 and line2:
    content1 = line1.strip().split('\t')
    content2 = line2.strip().split('\t')
    count1.append(int(content1[2]))
    list1.extend([int(content1[1])]*int(content1[2]))
    count2.append(int(content2[2]))
    list2.extend([int(content2[1])]*int(content2[2]))
    pos.append(int(content1[1]))
    line1 = file1.readline()
    line2 = file2.readline()

file1.close()
file2.close()
print("Reading data takes {:.3f}".format(time.process_time()-start_read))

line1, = plt.plot(pos,count1,color='blue')
line2, = plt.plot(pos,count2,color='green')

start_fit = time.process_time()
kde1 = KernelDensity(kernel=kernel_toUse, bandwidth=100).fit(np.array(list1)[:,np.newaxis])
kde2 = KernelDensity(kernel=kernel_toUse, bandwidth=50).fit(np.array(list1)[:,np.newaxis])
kde3 = KernelDensity(kernel=kernel_toUse, bandwidth=10).fit(np.array(list1)[:,np.newaxis])
print("Fitting the kernel density model takes {:.3f}".format(time.process_time()-start_fit))
start_eval = time.process_time()
log_dens1 = kde1.score_samples(np.array(pos)[:,np.newaxis])
log_dens2 = kde2.score_samples(np.array(pos)[:,np.newaxis])
log_dens3 = kde3.score_samples(np.array(pos)[:,np.newaxis])
print("Evaluate at specified points takes {:.3f}".format(time.process_time()-start_eval))
line3, = plt.plot(np.array(pos)[:,np.newaxis],len(list1)*np.exp(log_dens1),color='red')
line4, = plt.plot(np.array(pos)[:,np.newaxis],len(list1)*np.exp(log_dens2),color='black')
line5, = plt.plot(np.array(pos)[:,np.newaxis],len(list1)*np.exp(log_dens3),color='magenta')
plt.xlabel('Position along the Chromosome')
plt.ylabel('Read Counts')
plt.title('Distribution of DNA Reads along the Chromosome')
plt.legend((line1,line2,line3,line4,line5),(name1,name2,
   '{} kernel density with bandwidth 100'.format(kernel_toUse),
   '{} kernel density with bandwidth 50'.format(kernel_toUse),
   '{} kernel density with bandwidth 10'.format(kernel_toUse)),loc='upper right',fontsize='x-large')
plt.show()
plt.savefig("Distribution of DNA Reads along the Chromosome with Kernel Density Estimate")


