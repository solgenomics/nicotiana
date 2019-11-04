setwd("C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/WGCNA")
options(stringsAsFactors=FALSE)
library(WGCNA)
lnames = load(file="./round5/Ntab.network.RData")
# create a mapping from geneID to its TF family
library(hash)
h = hash()
TFs = read.csv('../../TFs/Ntab.hmmer.table.per.seq.summary.csv', header=TRUE)
for(TFfamily in TFs$TF.family.name){
  TFdata <- TFs[which(TFs$TF.family.name == TFfamily),
                names(TFs) %in% c("TF.family.name","number.of.TFs","list.of.TFs")]
  
  h[[TFfamily]] = strsplit(strsplit(TFdata$list.of.TFs,',')[[1]],' ')
}


color.set = unique(bwModuleColors)
data = matrix(NA,nrow=length(color.set),ncol=length(h))
#colnames(data) = TFs$TF.family.name
for (color.index in 0:(length(color.set)-1)){
  gene.subset = colnames(NtabExpr[,bwnet$colors==color.index])
  count = c()
  for(TFfamily in TFs$TF.family.name){
    #list.of.TFs = strsplit(h[[TFfamily]],',')
    count = c(count, sum(h[[TFfamily]] %in% gene.subset))
  }
  data[color.index+1,] = count
}

df = as.data.frame(data, row.names=1:length(color.set))
colnames(df) = TFs$TF.family.name
write.csv(df,"combined.module.TF.summary.csv")


#compare the root and leaf modules



