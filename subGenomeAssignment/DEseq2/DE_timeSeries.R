source("C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/GO/goModule.R")

fdr = 0.05
dir.htseq = "C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/DEseq2/htseq/NC95 vs NC95_nic(meJA)"
files.htseq = grep("htseq",list.files(dir.htseq),value=TRUE)
time = c("T0","T0","T0","T2","T2","T2","T6","T6","T24","T24","T24",
         "T0","T0","T2","T2","T6","T6","T24","T24")
genotype = c("nc95","nc95","nc95","nc95","nc95","nc95","nc95","nc95","nc95","nc95","nc95",
             "nc95_nic","nc95_nic","nc95_nic","nc95_nic","nc95_nic","nc95_nic","nc95_nic","nc95_nic")

sampleTable = data.frame(sampleName = files.htseq,
                          fileName = files.htseq,
                          time = time, strain = genotype)

library("DESeq2")
ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = dir.htseq,
                                       design= ~ time + strain + strain:time)

# The LRT test here identifies genes that show different expression profile
# between the two genotypes (NC95&NC95_nic) over time (at one or multiple time points)
dds.LRT = DESeq(ddsHTSeq, test="LRT", reduced = ~ time + strain)
res.LRT = results(dds.LRT)
res.LRT.ordered = res.LRT[order(res.LRT$padj),]
write.csv(as.data.frame(res.LRT.ordered), file="NC95.vs.NC95_nic.LRT.csv")
sig.gene.NC95.vs.NC95_nic = rownames(res.LRT[which(res.LRT$padj < fdr),])
go.sig.NC95.vs.NC95_nic = go(sig.gene.NC95.vs.NC95_nic)

# Plot normalized counts of an individual gene over time
library(ggplot2)
data <- plotCounts(dds.LRT, which.min(res.LRT$padj),
                   intgroup=c("time","strain"), returnData=TRUE)
                              ggplot(data, aes(x=time, y=count, color=strain, group=strain)) +
                              geom_point() + stat_smooth(se=FALSE,method="loess") + scale_y_log10()


# log2 fold change using Wald Tests at individual time points
res.T2 = results(dds.LRT, name="timeT2.strainnc95_nic", test="Wald")
res.T6 = results(dds.LRT, name="timeT6.strainnc95_nic", test="Wald")
res.T24 = results(dds.LRT, name="timeT24.strainnc95_nic", test="Wald")

                              
betas = coef(dds.LRT)
library("pheatmap")
topGenes = head(order(res.LRT$padj),200)
#exclude columns that store betas for "Intercept" and "strain_nc95_nic_vs_nc95"
#only columns for time points and the interaction between time:strain are retained
# use colnames(betas[topGenes,]) to check this
mat = betas[topGenes, -c(1,5)]
thr = 3
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),cluster_col=FALSE)
                              
                              
                              
                              
                                                   
