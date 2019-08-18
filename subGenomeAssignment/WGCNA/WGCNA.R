library(WGCNA)
setwd("C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/WGCNA/round3")
options(stringsAsFactors=FALSE)
NtabData = read.csv("rlog.counts.csv")
# each row is the gene expression data (rlog transformed) of all genes in each condition
NtabExpr = as.data.frame(t(NtabData[,-1]))
names(NtabExpr) = NtabData$X
rownames(NtabExpr) = colnames(NtabData[,-1])

# first cluster samples to see whether there are any outliers
sampleTree = hclust(dist(NtabExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#it appears that there is no obvious outliers. So nothing is done.
# choose appropriate soft-thresholding powers
# Choose a set of soft-thresholding powers
powers = 1:20
# Call the network topology analysis function
sft = pickSoftThreshold(NtabExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
pdf(file = "Plots/soft-thresholding.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
abline(h=0.85,col="black")
abline(h=0.80,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# from the plot, it appears that 12 is a good candidate for soft-thresholding power

bwnet = blockwiseModules(NtabExpr, power = 12,maxBlockSize=14000,
                       TOMType = "unsigned", minModuleSize = 20,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "Ntab-blockwise",networkType="signed",
                       verbose = 3)

# plotting
# open a graphics window
sizeGrWindow(12,12)
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 2
plotDendroAndColors(bwnet$dendrograms[[2]], bwModuleColors[bwnet$blockGenes[[2]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 2",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 3
plotDendroAndColors(bwnet$dendrograms[[3]], bwModuleColors[bwnet$blockGenes[[3]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 3",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

bwModuleLabels = bwnet$colors
bwModuleColors = labels2colors(bwnet$colors)
bwMEs = bwnet$MEs
geneTree1 = bwnet$dendrograms[[1]]
geneTree2 = bwnet$dendrograms[[2]]
geneTree3 = bwnet$dendrograms[[3]]

save(NtabExpr, bwnet, bwModuleLabels, bwModuleColors, bwMEs, 
     geneTree1, geneTree2, geneTree3, file="Ntab.network.RData")


# interfacing WGCNA with other tools
