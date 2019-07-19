
import subprocess

def compileGene(geneFile, chrom):
    dict = {}
    with open(geneFile) as geneF:
        with open('gene.{}.bed.temp'.format(chrom),'w') as bed:
            line = geneF.readline()
            while line:
                content = line.strip().split('\t')
                gene_name = content[0]
                start,end = int(content[2]),int(content[3])
                #gene_length = content[4]
                dict[gene_name] = (start, end)
                bed.write('{}\t{}\t{}\t{}\n'.format(chrom, start, end, gene_name))
                line = geneF.readline()
    return dict


def classify(geneFile, blockAssignment, chrom):
    # split the geneFile by chromosome so that gene.bed is also split by chromosome
    subprocess.call('grep \'{}\' {} > {}.{}.temp'.format(chrom, geneFile, geneFile, chrom), shell=True)
    geneLen_dict = compileGene('{}.{}.temp'.format(geneFile, chrom), chrom) # the dictionary stores key value pair as <gene_name, gene_length>
    tol_num_gene = len(geneLen_dict)
    subprocess.call('grep \'{}\' {} | grep \'Nsyl\' > Nsyl.{}.temp'.
                    format(chrom, blockAssignment, chrom), shell=True)
    subprocess.call('grep \'{}\' {} | grep \'Ntom\' > Ntom.{}.temp'.
                    format(chrom, blockAssignment, chrom), shell=True)
    subprocess.call('bedtools intersect -a gene.{}.bed.temp -b Nsyl.{}.temp Ntom.{}.temp -wo > inter.{}.temp'.
                    format(chrom, chrom, chrom, chrom), shell=True)
    with open('gene.{}.Nsyl.bed.temp'.format(chrom),'w') as Nsyl:
        with open('gene.{}.Ntom.bed.temp'.format(chrom),'w') as Ntom:
            with open('{}.report.temp'.format(chrom),'w') as report:
                # to be compatible with bedtools getfasta -name for downstream analysis
                # for format for the bed file must be:
                # <chromosome name> <start pos> <end pos>   <gene name>\n
                interF = open('inter.{}.temp'.format(chrom),'r')
                interL = interF.readline()
                genComp_dict = {}
                while interL:
                    content = interL.strip().split('\t')
                    gene = content[3]
                    if not gene in genComp_dict:
                        # in the tuple, the first integer is the number of base pairs that belong to Nsyl
                        # the second integer is the number of base pairs that belong to Ntom
                        genComp_dict[gene] = [0,0]
               
                    if content[-2] == 'Nsyl':
                        genComp_dict[gene][0] += int(content[-1])
                    else:
                        genComp_dict[gene][1] += int(content[-1])
                
                    interL = interF.readline()

                geneCount_Nsyl = 0
                geneCount_Ntom = 0
                for gene,comp in genComp_dict.items():
                    num_Nsyl = comp[0]
                    num_Ntom = comp[1]
                    pos = geneLen_dict[gene]
                    genLen = pos[1] - pos[0]
                    if num_Nsyl/genLen > 0.5 and (num_Ntom == 0 or num_Nsyl/num_Ntom > 2): #version 1
                    #if num_Nsyl/genLen > 0.8 and num_Ntom/genLen < 0.1: #version 2
                        Nsyl.write('{}\t{}\t{}\t{}\n'.format(chrom, pos[0], pos[1], gene))
                        geneCount_Nsyl += 1
                    elif num_Ntom/genLen > 0.5 and (num_Nsyl == 0 or num_Ntom/num_Nsyl > 2): #version 1
                    #if num_Ntom/genLen > 0.8 and num_Nsyl/genLen < 0.1: #version 2
                        Ntom.write('{}\t{}\t{}\t{}\n'.format(chrom, pos[0], pos[1], gene))
                        geneCount_Ntom += 1
                report.write('{} has a total of {} annotated genes\n'.format(chrom, tol_num_gene))
                report.write('{}({:.3f}%) are assigned as belong to S subgenome\n'.format(geneCount_Nsyl, 100*geneCount_Nsyl/tol_num_gene))
                report.write('{}({:.3f}%) are assigned as belong to T subgenome\n'.format(geneCount_Ntom, 100*geneCount_Ntom/tol_num_gene))

