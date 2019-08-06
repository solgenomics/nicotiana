from reportlab.lib.units import cm
from Bio.Graphics import BasicChromosome
from optparse import OptionParser
from reportlab.lib import colors
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys

# python visualAssign.py -a <subgenome assignment at block level> -s <chrom size txt file>
# -n <subgenome name 1>,<subgenome name 2>
parser = OptionParser()
parser.add_option("-a",action="store",type="string",dest="assignmentFile",
                  help="Must be in block form. Produced by assignBase.py")
parser.add_option("-s","--size",action="store",type="string",dest="chromSize",
                  help="Chromosome size file")
parser.add_option("-n","--name",action="store",type="string",dest="subgenome_name",
                  help="Subgenome names. Delimited by comma.")
(options,args) = parser.parse_args()
subG1 = options.subgenome_name.split(',')[0]
subG2 = options.subgenome_name.split(',')[1]

# Get color RGB value from https://www.rapidtables.com/web/color/RGB_Color.html
def getFeature(start, end, id):
    if id == subG1:
        customColor = colors.Color(red=(0/255),green=(128.0/255),blue=(255.0/255))
        return (start,end,0,"",customColor,customColor)
    elif id == subG2:
        customColor = colors.Color(red=(255.0/255),green=(51.0/255),blue=(51.0/255))
        return (start,end,0,"",customColor,customColor)
    elif id == "ambiguous":
        customColor = colors.Color(red=(60.0/255),green=(179.0/255),blue=(113.0/255))
        return (start,end,0,"",customColor,customColor)
    elif id == "unassigned":
        customColor = colors.Color(red=(46.0/255),green=(139.0/255),blue=(87.0/255))
        return (start,end,0,"",customColor,customColor)
    else:
        print(f"Illegal block assignment.id={id}")
        sys.exit()

#parse the chromosome size file first
#store chromosome size in a dictionary
chromsize = {}
chromsizeFile = open(options.chromSize,'r')
chromsizeLine = chromsizeFile.readline()
while chromsizeLine:
    #format of chromsize file:<chromosome number>\t<chromsome size>\n
    content = chromsizeLine.strip().split('\t')
    chromsize[content[0]] = int(content[1])
    chromsizeLine = chromsizeFile.readline()

chromsizeFile.close()

maxlen = max(chromsize.values())
telomere_len = maxlen//100

chr_diagram = BasicChromosome.Organism()
chr_diagram.page_size = (29.7*cm, 21*cm) #A4 size

assignFile = open(options.assignmentFile,'r')
line = assignFile.readline()
for name,size in chromsize.items():
    cur_chromosome = BasicChromosome.Chromosome(name)
    cur_chromosome.scale_num = maxlen + 2*telomere_len

    #add an opening telomere
    start = BasicChromosome.TelomereSegment()
    start.scale = telomere_len
    cur_chromosome.add(start)

    # add blocks
    while line:
        content = line.strip().split('\t')
        if content[0] != name:
            break
        else:
            start = int(content[1])
            end = int(content[2])
            id = content[3]
            # maybe instead of adding bodies one at a time
            # you can produce a list of features here
            # must put feature tuple(s)in a list, otherwise you get an error
            #print(f"end={end}")
            #print(f"start={start}")
            body = BasicChromosome.AnnotatedChromosomeSegment(end-start,[getFeature(0,end-start,id)]) 
            body.scale = end-start
            cur_chromosome.add(body)

        line = assignFile.readline()


    #add a closing telomere
    end = BasicChromosome.TelomereSegment(inverted=True)
    end.scale = telomere_len
    cur_chromosome.add(end)

    #done with this chromosome
    chr_diagram.add(cur_chromosome)

chr_diagram.draw("Nicotiana.tabacum.subgenome.pdf","N.tabacum")







