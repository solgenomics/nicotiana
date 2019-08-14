setwd("C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/WGCNA")
options(stringsAsFactors=FALSE)
library(WGCNA)
lnames = load(file="./round1/Ntab.network.2.RData")

Glist = function(Infile){
  glist = list()
  con = file(Infile,'r')
  while(TRUE){
    line = readLines(con,n=1)
    if(length(line)==0){
      break
    }
    glist = c(glist, trimws(line))
  }
  return (glist)
}

#first get a list of genes in each subgenome
Nsyl.list = Glist("gene.Nsyl.v1.ID")
Ntom.list = Glist("gene.Ntom.v1.ID")
color.set = unique(bwModuleColors)
colnames = c("total","S subgenome","T subgenome", "ambiguous")
data = matrix(NA,nrow=length(color.set),ncol=length(colnames))
# get the number of nodes in each module

for (color.index in 0:(length(color.set)-1)){
  count = sum(bwnet$colors == color.index)
  gene.subset = colnames(NtabExpr[,bwnet$colors==color.index])
  Nsyl.count = sum(Nsyl.list %in% gene.subset)
  Ntom.count = sum(Ntom.list %in% gene.subset)
  other.count = count - Nsyl.count - Ntom.count
  data[color.index+1,] = c(count,Nsyl.count,Ntom.count,other.count)
}

df = as.data.frame(data, row.names=1:length(color.set))
colnames(df) = c("total","S subgenome","T subgenome","ambiguous")
write.csv(df,"WGCNA_summaryBySubgenome.2.csv")

source("C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/GO/goModule.R")
for (color in color.set){
  gene.subset = colnames(NtabExpr[,bwModuleColors == color])
  go.result = go(gene.subset, FALSE)
  write.csv(go.result, paste('C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/WGCNA/round1/GO.network.',color,'.csv',sep=""))
}