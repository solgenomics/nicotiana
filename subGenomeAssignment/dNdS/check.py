# This script tries to check how many HEB gene pairs have the same reciproacl best hits against tomato genes

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-p","--pair",action="store",type="string",dest="HEB",
                  help="Path to the file that specifies HEB gene pairs")
parser.add_option("-b","--blast",action="store",type="string",dest="blast",
                  help ="Path to the file that specifies recirpocal best hits based on blast result")
(options,args) = parser.parse_args()


pair = {}
with open(options.HEB) as HEB:
    line = HEB.readline()
    while line:
        g1, g2 = line.strip().split('\t')
        pair[g1] = g2
        line = HEB.readline()

rbh = {}
with open(options.blast) as blast:
    line = blast.readline()
    while line:
        gN,gT = line.strip().split('\t')
        rbh[gN] = gT
        line = blast.readline()

countS = 0
countN = 0
out = open('consensus.out','w')
for gS,gN in pair.items():
    if gS in rbh:
        countS += 1
    if gN in rbh:
        countN += 1


out.write(f'CountS:{countS}')
out.write(f'CountN:{countN}')
out.close()