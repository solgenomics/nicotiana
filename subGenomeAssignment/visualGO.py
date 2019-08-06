import matplotlib.pyplot as plt
import numpy as np

def unpack(tupleL):
    des = [item[0] for item in tupleL]
    p = [item[1] for item in tupleL]
    return (des, p)

# Visualize GO
def visual(BP, MF, CC, tissue, bias):
    figure = plt.figure()
    des_BP, p_BP = unpack(BP)
    des_MF, p_MF = unpack(MF)
    des_CC, p_CC = unpack(CC)
    width = 0.25
    plt.title(f'GO terms: Genes biased towards {bias} in {tissue}')
    posCC = np.flip(np.arange(len(CC)), axis=0)
    posMF = np.flip(len(CC) + np.arange(len(MF)), axis=0)
    posBP = np.flip(len(CC) + len(MF) + np.arange(len(BP)), axis=0)
    plt.barh(posBP, -np.log10(p_BP), width, color='#FF8000')
    plt.barh(posMF, -np.log10(p_MF), width, color='#FF3333')
    plt.barh(posCC, -np.log10(p_CC), width, color='#FFB266')
    plt.yticks(np.arange(len(BP)+len(MF)+len(CC)), des_CC[::-1]+des_MF[::-1]+des_BP[::-1], fontsize='medium')
    plt.xlabel('-$\log_{10}(p\_value)$')
    plt.tight_layout()
    plt.show()
    #plt.savefig(f'go.{tissue}.{bias}.jpg', dpi=300, format='jpg', quality=95)




## Root, biased towards Nsyl
#GoList_BP = [('hexose metabolic process', 0.00049),('carbohydrate metabolic process', 0.00056),
#             ('fucose metabolic process', 0.00058),('monosaccharide metabolic process', 0.00067),
#             ('transmembrane transport', 0.0085),('chromatin-mediated maintenance of transcription', 0.01303)]
#GoList_MF = [('oxidoreductase activity (CH-OH donor group)', 0.00036),('transmembrane transporter activity', 0.00046),
#             ('oxidoreductase activity (NAD/NADP acceptor)', 0.00047),('transferase activity', 0.00097),
#             ('hydrolase activity', 0.00494),('inositol-1,4,5-trisphosphate 6-kinase activity', 0.00497)]
#GoList_CC = [('SAGA complex', 0.0014),('chloroplast thylakoid lumen', 0.0071),
#             ('plastid thylakoid lumen', 0.0071),('cytosolic ribosome', 0.0071)]
#visual(GoList_BP, GoList_MF, GoList_CC, 'Root', 'Nsyl')

## Root, biased towards Ntom
#GoList_BP = [('localization', 7.6e-6),('detoxification', 3.7e-5),
#             ('serine family amino acid metabolic process', 9.9e-5),('sulfur amino acid metabolic process', 0.00012),
#             ('cysteine metabolic process', 0.00064),('small molecule metabolic process', 0.00065)]
#GoList_MF = [('transmembrane transporter activity', 1e-7),('cofactor binding', 2.4e-6),
#             ('coenzyme binding', 6.3e-5),('antioxidant activity', 6.3e-5),
#             ('oxidoreductase activity', 7.4e-5),('small molecule binding', 0.00026)]
#GoList_CC = [('ribonuclease P complex', 0.0061),('H4/H2A histone acetyltransferase complex', 0.0121),
#             ('endoribonuclease complex', 0.0166),('endonuclease complex', 0.0166)]
#visual(GoList_BP, GoList_MF, GoList_CC, 'Root', 'Ntom')

#  Genes that show consistent bias towards Nsyl across all tissues
GoList_BP = [('organonitrogen compound \ncatabolic process', 0.00038),('peptide biosynthetic process', 0.00060),
             ('dicarboxylic acid metabolic process', 0.00064),('cellular amide metabolic process', 0.00078),
             ('peptide metabolic process', 0.00083),('protein catabolic process', 0.00088)]
GoList_MF = [('structural constituent of ribosome', 3.7e-6),('structural molecule activity', 3.2e-5),
             ('hydrolase activity', 0.0056),('betaine-aldehyde \ndehydrogenase activity', 0.0066),
             ('dTDP-4-dehydrorhamnose \n3,5-epimerase activity', 0.0098),('serine binding', 0.0098)]
GoList_CC = [('ribosome', 7.9e-5),('ribonucleoprotein complex', 0.00083),
             ('protein-containing complex', 0.00106),('non-membrane-bounded organelle', 0.00159)]
visual(GoList_BP, GoList_MF, GoList_CC, 'All Tissue', 'Nsyl')

# Genes that show consistent bias towards Ntom across all tissues

GoList_BP = [('establishment of localization', 7.9e-6),('organic substance transport', 4.8e-5),
             ('macromolecule localization', 0.00033),('vesicle-mediated transport', 0.00041),
             ('cellular modified amino acid \nmetabolic process', 0.00053),('methionine metabolic process', 0.00061)]
GoList_MF = [('symporter activity', 0.0018),('oxidoreductase activity', 0.0023),
             ('structural molecule activity', 0.0045),('S-methyl-5-thioribose-1-phosphate \nisomerase activity', 0.0062),
             ('unfolded protein binding', 0.0065), ('dihydroneopterin aldolase activity', 0.0093)]
GoList_CC = [('cytoplasm', 1.2e-5),('mitochondrial protein complex', 0.00234),
             ('intracellular', 0.00242),('mitochondrial membrane part', 0.00368)]
visual(GoList_BP, GoList_MF, GoList_CC, 'All Tissue', 'Ntom')