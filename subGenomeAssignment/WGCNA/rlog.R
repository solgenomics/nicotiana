library(DESeq2)
rawCounts = read.csv("./Round5/rawCounts.csv", header=TRUE)
rawMatrix = as.matrix(rawCounts[,-1])
rownames(rawMatrix)=rawCounts$X
tMatrix = rlog(rawMatrix, fitType="parametric", blind=TRUE)
rownames(tMatrix) = rownames(rawMatrix)
write.csv(tMatrix, "rlog.counts.csv")
