source("C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/GO/goModule.R")

fdr = 0.1
dir.htseq = "C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/DEseq2/htseq/NC95 vs NC95_nic(meJA)"
files.htseq = grep("htseq",list.files(dir.htseq),value=TRUE)
time = c()
genotype = c()

sampleTable <- data.frame(sampleName = files.htseq,
                          fileName = files.htseq,
                          time = time, strain = genotype)

library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = dir.htseq,
                                       design= ~ time + treatment)

