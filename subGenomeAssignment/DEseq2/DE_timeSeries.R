source("C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/GO/goModule.R")

fdr = 0.05
dir.htseq = "C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/DEseq2/htseq/NC95 vs NC95_nic(meJA)"
files.htseq = grep("htseq",list.files(dir.htseq),value=TRUE)
time = c("T0","T0","T0","T0120","T0120","T0120","T0480","T0480","T1440","T1440","T1440",
         "T0","T0","T0120","T0120","T0480","T0480","T1440","T1440")
genotype = c("nc95","nc95","nc95","nc95","nc95","nc95","nc95","nc95","nc95","nc95","nc95",
             "nc95_nic","nc95_nic","nc95_nic","nc95_nic","nc95_nic","nc95_nic","nc95_nic","nc95_nic")

sampleTable = data.frame(sampleName = files.htseq,
                          fileName = files.htseq,
                          time = time, strain = genotype)

library("DESeq2")
ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = dir.htseq,
                                       design= ~ time + strain + strain:time)

# Do some exploratory analysis first
rld = rlog(ddsHTSeq)
# calculate Euclidean distance between each sample
sampleDist = dist(t(assay(rld)))
# draw heatmap of sample distance
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDist )
rownames(sampleDistMatrix) <- paste( genotype, time, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDist,
         clustering_distance_cols=sampleDist,
         col=colors)





# The LRT test here identifies genes that show different expression profile
# between the two genotypes (NC95&NC95_nic) over time (at one or multiple time points)
dds.LRT = DESeq(ddsHTSeq, test="LRT", reduced = ~ time + strain)
res.LRT = results(dds.LRT)
res.LRT.ordered = res.LRT[order(res.LRT$padj),]
write.csv(as.data.frame(res.LRT.ordered), file="NC95.vs.NC95_nic.LRT.csv")
sig.gene.NC95.vs.NC95_nic = rownames(res.LRT[which(res.LRT$padj < fdr),])
go.sig.NC95.vs.NC95_nic = go(sig.gene.NC95.vs.NC95_nic)

# Plot normalized counts of an individual gene over time
num.sig = sum(res.LRT$padj <= fdr,na.rm=TRUE)
library(ggplot2)
# center the title
# since ggplot 2.0, the title is left-aligned by default
theme_update(plot.title = element_text(hjust = 0.5))
for (i in 1:num.sig){
  data <- plotCounts(dds.LRT, rownames(res.LRT.ordered[i,]),
                   intgroup=c("time","strain"), returnData=TRUE)
  ggplot(data, aes(x=time, y=count, color=strain, group=strain)) +
  geom_point() + stat_smooth(se=FALSE,method="loess") + scale_y_log10() + ggtitle(rownames(res.LRT.ordered[i,]))
  ggsave(paste(rownames(res.LRT.ordered[i,]),'.png',sep=""),
  path='C:/Users/10453/Desktop/2019 Summer Research/alkaloid RNA-seq/plotCountsNC95vsNC95_nic',
  width=6.78,height=3.12)
}


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
                              
                              
                              
                              
                                                   
