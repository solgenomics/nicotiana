
### GO enrichment analysis
source("C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/GO/goModule.R")
### Files that provide list of genes should reside in the following directory.
### Otherwise please change to the appropriate working directory.
setwd("C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/GO")
go.res = goFromFile("Intersection of genes towards Nsyl.txt", 1)
#go.res = goFromCSV("HEB.meJA_T6.towards.Ntom.csv", 2, 0.05)
#write.csv(go.res,'meJA.T6.Root.Ntom.BP.csv')
