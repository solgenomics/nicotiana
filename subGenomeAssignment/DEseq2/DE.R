## DE analysis using DESeq2
## input data format: htseq-count output
go = function(gList){
  ### GO enrichment analysis
  setwd("C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/GO")
  library(topGO)
  geneID2GO <- readMappings(file = "Ntab.gene2go.newAnnot")
  # the gene universe is given by the following line
  # A total of 34056 genes get one or more GO terms assigned to it
  geneNames = names(geneID2GO)
  gene.of.interest = gList
  gene.of.interest = as.integer(geneNames %in% gene.of.interest)
  gene.of.interest = factor(gene.of.interest)
  names(gene.of.interest) = geneNames # make it intno a named vector
  GOdata = new("topGOdata", ontology="BP", allGenes=gene.of.interest,
               annot=annFUN.gene2GO, gene2GO=geneID2GO)
  test.stat = new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher = getSigGroups(GOdata, test.stat)
  pvals = score(resultFisher)
  allRes = GenTable(GOdata, classic = resultFisher,
                    orderBy = "classic", ranksOf = "classic", topNodes = 20)
  return (allRes)
}

fdr = 0.1
dir.htseq = "C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/DEseq2/htseq"
files.htseq = grep("htseq",list.files(dir.htseq),value=TRUE)
condition = c("Control_T0","Control_T0","Control_T0",
              "Control_T2","Control_T2","Control_T2",
              "meJA_T2","meJA_T2","meJA_T2",
              "Control_T6","Control_T6","Control_T6",
              "meJA_T6","meJA_T6",
              "Control_T24","Control_T24","Control_T24",
              "meJA_T24","meJA_T24","meJA_T24")
sampleTable <- data.frame(sampleName = files.htseq,
                          fileName = files.htseq,
                          condition = condition)
library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = dir.htseq,
                                       design= ~ condition)

## First let's compare meJA_T6 against Control_T6
## The reason for such a comparison:
## We observed that T6 has the highest number of biased gene paris; and the number
## of biased genes rises from T0 to T6 and drops from T6 to T24, so T6 might represent
## the climax of the effect of Jasmonate treatment.
ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "Control_T6")
dds = DESeq(ddsHTSeq)
result = results(dds)
sig.gene.t6.vs.t6 = rownames(result[which(result$padj < fdr),])
go.sig.t6.vs.t6 = go(sig.gene.t6.vs.t6)

resOrdered = result[order(result$padj),]
# write the results to output
write.csv(as.data.frame(resOrdered), file="Control_T6.vs.meJA_T6.csv")





## Now compare meJA_T2 against Control_T2
## This comparison seems to yield most DE genes
ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "Control_T2")
dds = DESeq(ddsHTSeq)
result = results(dds, contrast=c('condition','meJA_T2','Control_T2'))
sig.gene.t2.vs.t2 = rownames(result[which(result$padj < fdr),])
go.sig.t2.vs.t2 = go(sig.gene.t2.vs.t2)

resOrdered = result[order(result$padj),]
# write the results to output
write.csv(as.data.frame(resOrdered), file="Control_T2.vs.meJA_T2.csv")






##Now compare meJA_T24 against Control_T24
ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "Control_T24")
dds = DESeq(ddsHTSeq)
result = results(dds, contrast=c('condition','meJA_T24','Control_T24'))
sig.gene.t24.vs.t24 = rownames(result[which(result$padj < fdr),])
go.sig.t24.vs.t24 = go(sig.gene.t6.vs.t6)

resOrdered = result[order(result$padj),]
# write the results to output
write.csv(as.data.frame(resOrdered), file="Control_T24.vs.meJA_T24.csv")



