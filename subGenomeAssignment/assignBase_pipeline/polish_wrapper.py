# Given a subgenome assignment in blocks form, and a bed file for gaps(optional)
# produce a new subgenome assignment such that all gap regions are labelled as unassigned
# Also do some final polishing step for the subgeome assignment.


from optparse import OptionParser
import subprocess
import os
parser = OptionParser()
parser.add_option("--chrom",action="store",type="string",dest="chromName",
                  help="A file that contains a list of chromosome names. Must be consistent with the bedgraph file. One chromosome should occupy one line")
parser.add_option("-g","--gap",action="store",type="string",dest="gap",
                  help="a bed file describing the gaps in the genome assembly. All gap regions will be automatically labelled as unassigned.")
parser.add_option("-f","--fileIn",action="store",type="string",dest="assign",
                  help="A subgenome assignment file in block form.")
(options,args) = parser.parse_args()

#first parse the chromosome name file
chromFile = open(options.chromName,'r')
chromLine = chromFile.readline()
chrom = []
while chromLine:
    chrom.append(chromLine.strip())
    chromLine = chromFile.readline()
chromFile.close()

#subprocess.call('bedtools intersect -a {} -b {} > gap.intersect'.format(options.assign, options.gap))
#split assinment file and gap file by chromosome
#assignFile = open(options.assign,'r')
#gapFile = open('gap.intersect','r')
process = []
for c in chrom:
    p1 = subprocess.Popen('grep \'{}\' {} > assign.{}.temp'.format(c, options.assign, c), shell=True)
    p2 = subprocess.Popen('grep \'{}\' {} > gap.{}.temp'.format(c, options.gap, c), shell=True)
    process.append(p1)
    process.append(p2)

for p in process:
    p.wait()

# now start polishing each individual chromosome
for c in chrom:
    p = subprocess.Popen('python {}/polish.py -f assign.{}.temp -g gap.{}.temp -c {}'.format(
        os.path.dirname(os.path.abspath(__file__)),c,c,c), shell=True)
    process.append(p)

for p in process:
    p.wait()


subprocess.call('cat *.polished.temp > blockAssignment.polished', shell=True)
subprocess.call('rm *.temp', shell=True)

