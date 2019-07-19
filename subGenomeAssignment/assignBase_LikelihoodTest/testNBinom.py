

from likelihood.estNBinom import nbinomModel
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f","--file",action="store",type="string",dest="filename",
                  help="path to two bedgraph files produced from bedtool genomecov -b, two file paths should be delimited by comma.")
(options,args) = parser.parse_args()

files = options.filename.split(',')
NBinom = nbinomModel(files[0],files[1])


