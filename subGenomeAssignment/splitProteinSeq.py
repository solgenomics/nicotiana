
# first build a dictionary <mRNA, gene that correspond to that mRNA> from gff
mRNA_dict = {}
with open('x4-nitab5-v1.all.maker.fixed.gff') as gff:
    line = gff.readline()
    while line:
        if not line.startswith('#'):
            gff_content = line.strip().split('\t')
            type = gff_content[2]
            if type == 'mRNA':
                addInfo = gff_content[8]
                start_gene = addInfo.find('Parent=')+7
                end_gene = addInfo[start_gene:].find(';')
                gene = addInfo[start_gene:][:end_gene] if end_gene != -1 else addInfo[start_gene:]

                start_mRNA = addInfo.find('Name=')+5
                end_mRNA = addInfo[start_mRNA:].find(';')
                mRNA = addInfo[start_mRNA:][:end_mRNA] if end_mRNA != -1 else addInfo[start_mRNA:]

                mRNA_dict[mRNA] = gene
        line = gff.readline()


#for mRNA,gene in mRNA_dict.items():
#    print('{}\t{}'.format(mRNA,gene))

# build a gene set for each subgenome
def build_gene_set(geneF, geneSet):
    with open(geneF) as f:
        line = f.readline()
        while line:
            if line.startswith('>'):
                geneSet.add(line.strip()[1:])
            line = f.readline()

geneNsyl = set()
geneNtom = set()
build_gene_set('gene.Nsyl.v1.fasta', geneNsyl)
build_gene_set('gene.Ntom.v1.fasta', geneNtom)

# Now build a protein sequence dictionary
protein_dict = {}
with open('nitab5-v1.all.maker.proteins.fasta') as profile:
    line = profile.readline()
    protein_name = ''
    protein_seq = ''
    first = True
    while line:
        if line.startswith('>'):
            if not first:
                gene_name = mRNA_dict[protein_name]
                protein_dict[gene_name] = protein_seq
            protein_name = line.strip().split(' ')[0][1:]
            #print(protein_name)
            protein_seq = ''
        else:
            protein_seq += line.strip()
        first = False
        line = profile.readline()
    protein_dict[mRNA_dict[protein_name]] = protein_seq

with open('Nsyl.protein.fasta','w') as NsylP:
    with open('Ntom.protein.fasta','w') as NtomP:
        for gene,seq in protein_dict.items():
            if gene in geneNsyl:
                NsylP.write('>{}\n'.format(gene))
                NsylP.write('{}\n'.format(seq))
            elif gene in geneNtom:
                NtomP.write('>{}\n'.format(gene))
                NtomP.write('{}\n'.format(seq))


