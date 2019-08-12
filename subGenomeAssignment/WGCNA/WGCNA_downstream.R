setwd("C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/WGCNA")
options(stringsAsFactors=FALSE)
library(WGCNA)
lnames = load(file="Ntab.network.RData")

nGenes = ncol(NtabExpr)
nSamples = nrow(NtabExpr)

# explore the similarities between module eigengenes
# Recalculate module eigengenes
MEs = moduleEigengenes(NtabExpr, bwModuleColors)$eigengenes
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MEs, "Eigengene Dendeogram", marDendro = c(0,4,1,2), 
                      plotHeatmaps = FALSE)
par(cex = 0.9)
plotEigengeneNetworks(MEs, "Eigengene heatmap", marHeatmap=c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle=90)

