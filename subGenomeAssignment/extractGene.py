# Extract genes from a given gff file
# geneName(use the field "NAME" in the last column of a gff file) chromosome start_position end_position gene_length

from optparse import OptionParser
import subprocess
import geneOps

parser = OptionParser()
parser.add_option("-f","--file",action="store",type="string",dest="gff",
                  help="path to the GFF file")
parser.add_option("-i","--id",action="store",type="string",dest="ID",default="gene",
                  help="The feature to extract from the GFF")
parser.add_option("-a","--assign",action="store",type="string",dest="assign",
                  help="Path to the subgenome assignment(block form)")
(options,args) = parser.parse_args()

# Divide a given gene set into two sets according to which subgenome it belongs
# Criterion for deciding whether a gene belongs to S or T subgenome:
# Given a gene G, if at least half of its length belongs to a particular subgenome and is at least
# twice as long as the region belong to the other subgenome

# Do one chromosome at a time


#def compileGene(geneFile, chrom):
#    dict = {}
#    with open(geneFile) as geneF:
#        with open('gene.{}.bed.temp'.format(chrom),'w') as bed:
#            line = geneF.readline()
#            while line:
#                content = line.strip().split('\t')
#                gene_name = content[0]
#                start,end = int(content[2]),int(content[3])
#                #gene_length = content[4]
#                dict[gene_name] = (start, end)
#                bed.write('{}\t{}\t{}\t{}\n'.format(chrom, start, end, gene_name))
#                line = geneF.readline()
#    return dict


#def classify(geneFile, blockAssignment, chrom):
#    # split the geneFile by chromosome so that gene.bed is also split by chromosome
#    subprocess.call('grep \'{}\' {} > {}.{}.temp'.format(chrom, geneFile, geneFile, chrom), shell=True)
#    geneLen_dict = compileGene('{}.{}.temp'.format(geneFile, chrom), chrom) # the dictionary stores key value pair as <gene_name, gene_length>
#    tol_num_gene = len(geneLen_dict)
#    subprocess.call('grep \'{}\' {} | grep \'Nsyl\' > Nsyl.{}.temp'.
#                    format(chrom, blockAssignment, chrom), shell=True)
#    subprocess.call('grep \'{}\' {} | grep \'Ntom\' > Ntom.{}.temp'.
#                    format(chrom, blockAssignment, chrom), shell=True)
#    subprocess.call('bedtools intersect -a gene.{}.bed.temp -b Nsyl.{}.temp Ntom.{}.temp -wo > inter.{}.temp'.
#                    format(chrom, chrom, chrom, chrom), shell=True)
#    with open('gene.{}.Nsyl.bed.temp'.format(chrom),'w') as Nsyl:
#        with open('gene.{}.Ntom.bed.temp'.format(chrom),'w') as Ntom:
#            with open('{}.report.temp'.format(chrom),'w') as report:
#                # to be compatible with bedtools getfasta -name for downstream analysis
#                # for format for the bed file must be:
#                # <chromosome name> <start pos> <end pos>   <gene name>\n
#                interF = open('inter.{}.temp'.format(chrom),'r')
#                interL = interF.readline()
#                genComp_dict = {}
#                while interL:
#                    content = interL.strip().split('\t')
#                    gene = content[3]
#                    if not gene in genComp_dict:
#                        # in the tuple, the first integer is the number of base pairs that belong to Nsyl
#                        # the second integer is the number of base pairs that belong to Ntom
#                        genComp_dict[gene] = [0,0]
               
#                    if content[-2] == 'Nsyl':
#                        genComp_dict[gene][0] += int(content[-1])
#                    else:
#                        genComp_dict[gene][1] += int(content[-1])
                
#                    interL = interF.readline()

#                geneCount_Nsyl = 0
#                geneCount_Ntom = 0
#                for gene,comp in genComp_dict.items():
#                    num_Nsyl = comp[0]
#                    num_Ntom = comp[1]
#                    pos = geneLen_dict[gene]
#                    genLen = pos[1] - pos[0]
#                    if num_Nsyl/genLen > 0.5 and (num_Ntom == 0 or num_Nsyl/num_Ntom > 2):
#                        Nsyl.write('{}\t{}\t{}\t{}\n'.format(chrom, pos[0], pos[1], gene))
#                        geneCount_Nsyl += 1
#                    elif num_Ntom/genLen > 0.5 and (num_Nsyl == 0 or num_Ntom/num_Nsyl > 2):
#                        Ntom.write('{}\t{}\t{}\t{}\n'.format(chrom, pos[0], pos[1], gene))
#                        geneCount_Ntom += 1
#                report.write('{} has a total of {} annotated genes\n'.format(chrom, tol_num_gene))
#                report.write('{}({:.3f}%) are assigned as belong to S subgenome\n'.format(geneCount_Nsyl, 100*geneCount_Nsyl/tol_num_gene))
#                report.write('{}({:.3f}%) are assigned as belong to T subgenome\n'.format(geneCount_Ntom, 100*geneCount_Ntom/tol_num_gene))

with open(options.gff) as gff:
    gff_line = gff.readline()
    with open('gff.gene','w') as out:
        while gff_line:
            if not gff_line.startswith('#'):
                content = gff_line.strip().split('\t')
                if content[2] == 'gene':
                    chrom = content[0]
                    start,end = int(content[3]),int(content[4])
                    #print(content[8])
                    pos_1 = content[8].find('Name=') + 5
                    pos_2 = content[8][pos_1:].find(';')
                    name = content[8][pos_1:pos_2+pos_1] if pos_2 != -1 else content[8][pos_1:]
                    #print(name)
                    out.write('{}\t{}\t{}\t{}\t{}\n'.format(name,chrom,start,end,end-start))

            gff_line = gff.readline()





#TODO:
#Add a brief summary for each chromosome
#Finish parallel processing of each chromosome

chromList = ['Nitab01','Nitab02','Nitab03','Nitab04','Nitab05','Nitab06','Nitab07','Nitab08',
             'Nitab09','Nitab10','Nitab11','Nitab12','Nitab13','Nitab14','Nitab15','Nitab16',
             'Nitab17','Nitab18','Nitab19','Nitab20','Nitab21','Nitab22','Nitab23','Nitab24']

for c in chromList:
    geneOps.classify('gff.gene', options.assign, c)

subprocess.call('cat *.report.temp > final.report', shell=True)
subprocess.call('cat *.Nsyl.bed.temp > gene.Nsyl.bed', shell=True)
subprocess.call('cat *.Ntom.bed.temp > gene.Ntom.bed', shell=True)
subprocess.call('rm *.temp', shell=True)
