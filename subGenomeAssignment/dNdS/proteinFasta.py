### Prepare input file for blastn search

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-g","--gff",action="store",type="string",dest="gff",
                  help="path to the .gff file")
parser.add_option("-p","--protein",action="store",type="string",dest="proteinFasta",
                  help="Path to the protein .fasta file")
(options,args) = parser.parse_args()

def attr(attribute):
    dict = {}
    for attr in attribute.strip().split(';'):
        tag, val = attr.split('=')
        dict[tag] = val
    return dict

mRNA2gene = {}
with open(options.gff) as gff:
    line = gff.readline()
    while line:
        if not line.startswith('#'):
            seqid, source, feature, start, end, score, strand, phase, attributes = line.strip().split('\t')
            if feature == 'mRNA':
                d = attr(attributes)
                p = d['Parent'].find(':')
                mRNA, gene = d['Name'], d['Parent'][p+1:]
                mRNA2gene[mRNA] = gene
                # this should work even if there is no :, because in that case, p=-1, and plus one makes it 0, which is what we want
        line = gff.readline()

proteinSeq = {}
with open(options.proteinFasta) as f:
    line = f.readline()
    count = 0
    while line:
        if line.startswith('>'):
            count += 1
            print(f"Process {count}")
            gene = mRNA2gene[line.strip().split(' ')[0][1:]]
            proSeq = ""
            line = f.readline()
            while not line.startswith('>'):
                proSeq += line.strip()
                line = f.readline()
                if not line:
                    proteinSeq[gene] = proSeq
                    break
            proteinSeq[gene] = proSeq


with open('out.protein.fasta','w') as out:
    for g,p in proteinSeq.items():
        out.write(f'>{g}\n')
        if '*' in p:
            out.write(f'{p[:-1]}\n')
        else:
            out.write(f'{p}\n')