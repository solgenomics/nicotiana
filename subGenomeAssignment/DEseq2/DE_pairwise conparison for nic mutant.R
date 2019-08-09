
## DE analysis using DESeq2
## input data format: htseq-count output
source("C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/GO/goModule.R")

fdr = 0.05
dir.htseq = "C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/DEseq2/htseq/NC95 vs NC95_nic(meJA)"
files.htseq = grep("htseq",list.files(dir.htseq),value=TRUE)
condition = c("nc95_T0","nc95_T0",
              "nc95_T0","nc95_T0120","nc95_T0120",
              "nc95_T0120","nc95_T0480","nc95_T0480",
              "nc95_T1440","nc95_T1440","nc95_T1440",
              "nc95_nic_T0","nc95_nic_T0",
              "nc95_nic_T0120","nc95_nic_T0120","nc95_nic_T0480",
              "nc95_nic_T0480","nc95_nic_T01440","nc95_nic_T1440")
sampleTable <- data.frame(sampleName = files.htseq,
                          fileName = files.htseq,
                          condition=condition)                          
library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = dir.htseq,
                                       design= ~ condition)

## First let's compare nicT0 against wt_T0
ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "nc95_T0")
dds = DESeq(ddsHTSeq)
result = results(dds, contrast=c('condition','nc95_nic_T0','nc95_T0'))
sig.gene.t0.vs.t0 = rownames(result[which(result$padj < fdr),])
go.sig.t0.vs.t0 = go(sig.gene.t0.vs.t0)

resOrdered = result[order(result$padj),]
# write the results to output
write.csv(as.data.frame(resOrdered), file="nc95_nic_T0.vs.nc95_T0.csv")



ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "nc95_T0120")
dds = DESeq(ddsHTSeq)
result = results(dds, contrast=c('condition','nc95_nic_T0120','nc95_T0120'))
sig.gene.t0120.vs.t0120 = rownames(result[which(result$padj < fdr),])
go.sig.t0120.vs.t0120 = go(sig.gene.t0120.vs.t0120)

resOrdered = result[order(result$padj),]
# write the results to output
write.csv(as.data.frame(resOrdered), file="nc95_nic_T0120.vs.nc95_T0120.csv")




ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "nc95_T0480")
dds = DESeq(ddsHTSeq)
result = results(dds, contrast=c('condition','nc95_nic_T0480','nc95_T0480'))
sig.gene.t0480.vs.t0480 = rownames(result[which(result$padj < fdr),])
go.sig.t0480.vs.t0480 = go(sig.gene.t0480.vs.t0480)

resOrdered = result[order(result$padj),]
# write the results to output
write.csv(as.data.frame(resOrdered), file="nc95_nic_T0480.vs.nc95_T0480.csv")






ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "nc95_T1440")
dds = DESeq(ddsHTSeq)
result = results(dds, contrast=c('condition','nc95_nic_T1440','nc95_T1440'))
sig.gene.t1440.vs.t1440 = rownames(result[which(result$padj < fdr),])
go.sig.t1440.vs.t1440 = go(sig.gene.t1440.vs.t1440)

resOrdered = result[order(result$padj),]
# write the results to output
write.csv(as.data.frame(resOrdered), file="nc95_nic_T1440.vs.nc95_T1440.csv")
