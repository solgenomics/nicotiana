# Assign each base to a subgenome A or B, or ambiguous, or unassigned
# python assignBase.py -f <bedgraph file1>,<bedgraph file2> 
# -n <subgenome name1>,<subgenome name2> -t <threshold,default is 2>
# -c <sequencing coverage 1>,<sequencing coverage 2>
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f","--file",action="store",type="string",dest="filename",
                  help="path to two bedgraph files produced from bedtool genomecov -b, two file paths should be delimited by comma.")
parser.add_option("-t",action="store",type="float",dest="threshold",help="set the threshold when assigning subgenomes. A base is assigned to subgenome1\
                if the ratio of coverage from bed1 to bed2 is greater or equal to threshold. Default value is 2.")
parser.add_option("-n","--name",action="store",type="string",dest="subgenome_name",
                  help="Name for the two subgenomes. Delimited by comma. Must be in the same order as in the bedfiles.")
parser.add_option("-c","--coverage",action="store",type="string",dest="coverage",help="Coverage of the DNA sequencing that produce the input bedgraph file. Must be in the same order as in the bedfiles. Delimited by comma.")
(options,args) = parser.parse_args()

threshold = 2
if options.threshold != None:
    threshold = options[threshold]

bedfile_list = (options.filename).split(',')
subgenome_name_list = (options.subgenome_name).split(',')
cov = options.coverage.split(',')
correction_coef = float(cov[1])/float(cov[0])

bed1 = open(bedfile_list[0],'r')
subG1_name = subgenome_name_list[0]
bed2 = open(bedfile_list[1],'r')
subG2_name = subgenome_name_list[1]
outFile = open("out.assign.singlebase.txt",'w')

bed1_line = bed1.readline()
bed2_line = bed2.readline()

count_G1 = 0
count_G2 = 0
count_ambi = 0
count_un = 0

while bed1_line and bed2_line: #although I checked both here, bed1 and bed2 are expected to have the
    # number of lines if they are correctly produced from bedtools genomecov -b.
    bed1_content = list(filter(None,bed1_line.strip().split('\t')))
    bed2_content = list(filter(None,bed2_line.strip().split('\t')))
    #print(bed1_content)
    #print(bed2_content)
    chrom = bed1_content[0] #chromosome name
    base_pos = bed1_content[1] #base position
    cov_G1 = int(bed1_content[2]) #coverage from bedfile1
    cov_G2 = int(bed2_content[2]) #coverage from bedfile2 
    

    if cov_G1 == 0 and cov_G2 != 0:
        count_G2 += 1
        outFile.write("{}\t{}\t{}\n".format(chrom,base_pos,subG2_name))
    elif cov_G1 != 0 and cov_G2 == 0:
        count_G1 += 1
        outFile.write("{}\t{}\t{}\n".format(chrom,base_pos,subG1_name))
    elif cov_G1 == 0 and cov_G2 == 0:
        count_un += 1
        outFile.write("{}\t{}\t{}\n".format(chrom,base_pos,"unassigned"))
    else:
        ratio = correction_coef*(cov_G1/cov_G2)
        if ratio >= threshold:
            count_G1 += 1
            outFile.write("{}\t{}\t{}\n".format(chrom,base_pos,subG1_name))
        elif ratio <= (1/threshold):
            count_G2 += 1
            outFile.write("{}\t{}\t{}\n".format(chrom,base_pos,subG2_name))
        else:
            count_ambi += 1
            outFile.write("{}\t{}\t{}\n".format(chrom,base_pos,"ambiguous"))

    bed1_line = bed1.readline()
    bed2_line = bed2.readline()

print("Assignment Done")
total = count_G1 + count_G2 + count_ambi + count_un
print("A total of {} bases processed".format(total))
print("{}({:.2f}%) bases are assigned to {}".format(count_G1,100*count_G1/total,subG1_name))
print("{}({:.2f}%) bases are assigned to {}".format(count_G2,100*count_G2/total,subG2_name))
print("{}({:.2f}%) bases are ambiguous".format(count_ambi,100*count_ambi/total))
print("{}({:.2f}%) bases are unassigned".format(count_un,100*count_un/total))

bed1.close()
bed2.close()
outFile.close()