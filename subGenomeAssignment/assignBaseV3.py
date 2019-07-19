# Assign each base to a subgenome A or B, or ambiguous, or unassigned
# python assignBase.py -f <bedgraph file1>,<bedgraph file2> 
# -n <subgenome name1>,<subgenome name2> -t <threshold,default is 2>
# -c <expected coverage 1>,<expected coverage 2> -k <kernel> -w <bandwidth>
from optparse import OptionParser
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
import numpy as np
import time

tol = 1e-4
parser = OptionParser()
parser.add_option("-f","--file",action="store",type="string",dest="filename",
                  help="path to two bedgraph files produced from bedtool genomecov -b, two file paths should be delimited by comma.")
parser.add_option("-t",action="store",type="float",dest="threshold",help="set the threshold when assigning subgenomes. A base is assigned to subgenome1\
                if the ratio of coverage from bed1 to bed2 is greater or equal to threshold. Default value is 2.")
parser.add_option("-n","--name",action="store",type="string",dest="subgenome_name",
                  help="Name for the two subgenomes. Delimited by comma. Must be in the same order as in the bedfiles.")
parser.add_option("-c","--coverage",action="store",type="string",dest="coverage",help="Coverage of the DNA sequencing that produce the input bedgraph file. Must be in the same order as in the bedfiles. Delimited by comma.")
parser.add_option("-k","--kernel",action="store",type="string",dest="kernel_to_use",
                  help="Kernel used in kernel density estimation")
parser.add_option("-w","--width",action="store",type="int",dest="width",help="bandwidth for KDE")
(options,args) = parser.parse_args()

threshold = 2
if options.threshold != None:
    threshold = options.threshold

kernel_to_use = "epanechnikov"
if options.kernel_to_use != None:
    kernel_to_use = options.kernel_to_use

width = 100
if options.width != None:
    width = options.width


bedfile_list = (options.filename).split(',')
subgenome_name_list = (options.subgenome_name).split(',')
cov = options.coverage.split(',')
correction_coef = float(cov[1])/float(cov[0])

bed1 = open(bedfile_list[0],'r')
subG1_name = subgenome_name_list[0]
bed2 = open(bedfile_list[1],'r')
subG2_name = subgenome_name_list[1]
outFile = open("out.KRassign.singlebase.txt",'w')
outFile2 = open("out.KRassign.block.txt",'w')

bed1_line = bed1.readline()
bed2_line = bed2.readline()

count_G1 = 0
count_G2 = 0
count_ambi = 0
count_un = 0

prev_chrom = ""
pos = []
list1 = []
list2 = []

def process_chrom(array1, array2, outfile_single, outfile_block, chrom):
    print("start processing a {}".format(chrom))
    global count_G1,count_G2,count_ambi,count_un
    start = 0
    prev_assign = ""
    for i in list(range(array1.size)):
        assign = ""
        l = max(0,i-width)
        r = min(1+i+width,len(pos))
        neigh = pos[l:r]
        neigh = (np.array(neigh)-i)/width
        K = 0.75*(np.ones(len(neigh),)-neigh*neigh)
        count1_adjusted = abs(np.array(list1[l:r]).dot(K)/K.sum())
        count2_adjusted = abs(np.array(list2[l:r]).dot(K)/K.sum())

        if count1_adjusted <= tol and count2_adjusted <= tol:
            count_un += 1
            outfile_single.write("{}\t{}\t{}\n".format(chrom, i+1, 'unassigned'))
            assign = 'unassigned'
        elif count1_adjusted <= tol and count2_adjusted > tol:
            count_G2 += 1
            outfile_single.write("{}\t{}\t{}\n".format(chrom, i+1, subG2_name))
            assign = subG2_name
        elif count1_adjusted > tol and count2_adjusted <= tol:
            count_G1 += 1
            outfile_single.write("{}\t{}\t{}\n".format(chrom, i+1, subG1_name))
            assign = subG1_name
        else:
            ratio = correction_coef*(count1_adjusted/count2_adjusted)
            if (ratio - threshold) >= tol:
                count_G1 += 1
                outfile_single.write("{}\t{}\t{}\n".format(chrom, i+1, subG1_name))
                assign = subG1_name
            elif (ratio - 1/threshold) <= tol:
                count_G2 += 1
                outfile_single.write("{}\t{}\t{}\n".format(chrom, i+1, subG2_name))
                assign = subG2_name
            else:
                count_ambi += 1
                outfile_single.write("{}\t{}\t{}\n".format(chrom, i+1, 'ambiguous'))
                assign = 'ambiguous'

        if start == 0: # deal with the initial case
                start = 1
                prev_assign = assign
            
        if assign != prev_assign:
                outfile_block.write("{}\t{}\t{}\t{}\n".format(chrom, start, i+1, prev_assign))
                start = i+1
                prev_assign = assign
    outfile_block.write("{}\t{}\t{}\t{}\n".format(chrom, start, array1.size+1, prev_assign))


while bed1_line and bed2_line: #although I checked both here, bed1 and bed2 are expected to have the
    # number of lines if they are correctly produced from bedtools genomecov -b.
    bed1_content = list(filter(None,bed1_line.strip().split('\t')))
    bed2_content = list(filter(None,bed2_line.strip().split('\t')))
    chrom = bed1_content[0] #chromosome name

    if chrom != prev_chrom:
        if prev_chrom != "":
            # you should consider using another thread to do the chromosome processing work
            # as process_chrom takes an awful amount of time
            process_chrom(np.array(list1), np.array(list2), outFile, outFile2, prev_chrom)
            pos = []
            list1 = []
            list2 = []
        prev_chrom = chrom

    base_pos = int(bed1_content[1]) #base position
    cov_G1 = int(bed1_content[2]) #coverage from bedfile1
    cov_G2 = int(bed2_content[2]) #coverage from bedfile2
    list1.append(cov_G1)
    list2.append(cov_G2)
    pos.append(base_pos)
    bed1_line = bed1.readline()
    bed2_line = bed2.readline()

process_chrom(np.array(list1),np.array(list2),outFile, outFile2, chrom)

print("Assignment Done")
total = count_G1 + count_G2 + count_ambi + count_un
print("A total of {} bases processed".format(total))
print("{}({:.4f}%) bases are assigned to {}".format(count_G1,100*count_G1/total,subG1_name))
print("{}({:.4f}%) bases are assigned to {}".format(count_G2,100*count_G2/total,subG2_name))
print("{}({:.4f}%) bases are ambiguous".format(count_ambi,100*count_ambi/total))
print("{}({:.4f}%) bases are unassigned".format(count_un,100*count_un/total))

bed1.close()
bed2.close()
outFile.close()
outFile2.close()


    