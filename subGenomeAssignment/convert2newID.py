### Convert files made from old .gff gene annotation id to new id
from optparse import OptionParser
import sys
from os import walk


parser = OptionParser()
parser.add_option("--htseq",action="store",type="string",dest="htseq",
                  help="Path to directories that contain htseq count files.")
parser.add_option("--go",action="store",type="string",dest="go",
                  help ="Path to the gene2GO mapping  file.")
(options,args) = parser.parse_args()

map_file = sys.argv[1]
map = {} # a dictionary <old ID, new ID>
with open(map_file) as mapF:
    line = mapF.readline()
    while line:
        new_id,_,_,qry = line.strip().split('\t')
        old_id = qry[0:qry.find('|')]
        map[old_id] = new_id
        line = mapF.readline()

if options.htseq != None:
    # convert htseq files using new gene IDs
    for dir in options.htseq.split(','):
        # extract all files in that sub-directory
        f = []
        for (dirpath, dirnames, filenames) in walk(dir):
            f.extend(filenames)
            break

        for file in f:
            with open(f'{dir}/{file}') as htseq_in:
                with open(f'{file}.newAnnot','w') as htseq_out:
                    line = htseq_in.readline()
                    while line:
                        if line.startswith('__'):
                            htseq_out.write(line)
                        else:
                            gene, count = line.strip().split('\t')
                            if gene not in map:
                                print(f'{gene} does not have a new ID')
                            else:
                                htseq_out.write(f'{map[gene]}\t{count}\n')
                        line = htseq_in.readline()

if options.go != None:
    with open(options.go) as go_in:
        with open(f'{options.go}.newAnnot','w') as go_out:
            line = go_in.readline()
            while line:
                gene, go_term = line.strip().split('\t')
                go_out.write(f'{map[gene]}\t{go_term}\n')
                line = go_in.readline()


            