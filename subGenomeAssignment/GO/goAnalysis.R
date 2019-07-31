
### GO enrichment analysis
source("C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/GO/goModule.R")
### Files that provide list of genes should reside in the following directory.
### Otherwise please change to the appropriate working directory.
setwd("C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/GO")
#go.res = goFromFile("HEB.Shoot_ZT12.towards.Ntom.txt", 2)
go.res = goFromCSV("HEB.Root_ZT12.towards.Ntom.csv", 2, 0.05)
write.csv(go.res,'Root.Ntom.BP.csv')
