import os
from optparse import OptionParser
import subprocess

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
parser.add_option("--chrom",action="store",type="string",dest="chromName",
                  help="A file that contains a list of chromosome names. Must be consistent with the bedgraph file. One chromosome should occupy one line")
parser.add_option("-g","--gap",action="store",type="string",dest="gap",
                  help="a bed file describing the gaps in the genome assembly. All gap regions will be automatically labelled as unassigned.")

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

#first parse the chromosome name file
chromFile = open(options.chromName,'r')
chromLine = chromFile.readline()
chrom = []
while chromLine:
    chrom.append(chromLine.strip())
    chromLine = chromFile.readline()
chromFile.close()

#split the bedgraph file by chromName
bedfile_list = (options.filename).split(',')
bed1,bed2 = bedfile_list[0],bedfile_list[1]
subgenome_name_list = (options.subgenome_name).split(',')
sub1,sub2 = subgenome_name_list[0],subgenome_name_list[1]
processes = []
for chromosome in chrom:
    p1 = subprocess.Popen('grep \'{}\' {} > {}.{}.temp'.format(chromosome, bed1, sub1, chromosome),shell=True)
    p2 = subprocess.Popen('grep \'{}\' {} > {}.{}.temp'.format(chromosome, bed2, sub2, chromosome),shell=True)
    processes.append(p1)
    processes.append(p2)

# a blocking call to all subprocesses so that when we proceed, we are sure that all temp files are generated
for p in processes:
    p.wait()

processes.clear()

coverage = options.coverage.strip().split(',')
cov1 = coverage[0]
cov2 = coverage[1]
for chromosome in chrom:
    p1 = subprocess.Popen('python {}/assignBaseV3_parallel.py -f {}.{}.temp,{}.{}.temp -n {},{} -t {} -c {},{} -k {} -w {} --chrom {}'
                          .format(os.path.dirname(os.path.abspath(__file__)),
                                  sub1, chromosome, sub2, chromosome, sub1, 
                                  sub2, threshold, cov1, cov2, kernel_to_use, width, chromosome), shell=True)
    processes.append(p1)
    
for p in processes:
    p.wait()

print("Done Processing All Chromosomes")
print("Now Concatenating Results")
subprocess.call('cat *.w{}.singlebase.temp > singlebaseAssignment.txt'.format(width), shell=True)
subprocess.call('cat *.w{}.block.temp > blockAssignment.txt'.format(width), shell=True)
subprocess.call('cat *.log.temp > Log.beforePolish',shell=True)
subprocess.call('rm *.temp', shell=True)

