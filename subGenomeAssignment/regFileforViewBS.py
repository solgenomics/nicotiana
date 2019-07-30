from optparse import OptionParser
import sys

parser = OptionParser()
parser.add_option("-p","--pair",action="store",type="string",dest="HEB",
                  help="Path to the file that specifies HEB gene pairs")
parser.add_option("-n","--name",action="store",type="string",dest="names",
                  help ="Names for the pair of genes provided. Eg, Nsyl, Ntom.")
(options,args) = parser.parse_args()

def featureDict(attribute):
    dict = {}
    attr_splitted = attribute.split(';')
    for p in attr_splitted:
        tag,value = p.split('=')
        dict[tag]=value
    return dict

class gene(object):
    def __init__(self, chrom, ID, start, end, strand):
        self.chrom = chrom
        self.ID = ID
        self.start = start
        self.end = end
        self.strand = strand

    def __str__(self):
        return f'{self.chrom}\t{self.start}\t{self.end}\t{self.ID}\t{self.strand}'

gffFile = sys.argv[1]
geneDic = {}
with open(gffFile) as gff:
    line = gff.readline()
    while line:
        if not line.startswith('#'):
            chr, _, feature, start, end, _, strand, _, attribute = line.strip().split('\t')
            if feature == 'gene':
                gID = featureDict(attribute)['ID']
                geneDic[gID] = gene(chr, gID, int(start),  int(end), strand)
        line = gff.readline()

n1, n2 = options.names.split(',')
with open(options.HEB) as pair:
    with open(f'{n1}.region','w') as out1:
        with open(f'{n2}.region','w') as out2:
            line = pair.readline()
            while line:
                g1, g2 = line.strip().split('\t')
                out1.write(f'{geneDic[g1]}\n')
                out2.write(f'{geneDic[g2]}\n')
                line = pair.readline()


