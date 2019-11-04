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
data_raw = matrix(NA,nrow=length(color.set),ncol=length(h))
data_test = matrix(NA, nrow=length(color.set), ncol=length(h))

total.gene.count = length(NtabExpr)
total.TFs.count = sum(TFs$number.of.TFs)
for (color.index in 0:(length(color.set)-1)){
  gene.subset = colnames(NtabExpr[,bwnet$colors==color.index])
  count.list = c()
  p.value.list = c()
  for(TFfamily in TFs$TF.family.name){
    overlap = sum(h[[TFfamily]] %in% gene.subset)
    count.list = c(count.list, overlap)
    data = matrix(data=c(overlap, length(gene.subset)-overlap, length(h[[TFfamily]])-overlap, 
                          total.gene.count-length(gene.subset)-length(h[[TFfamily]])+overlap), 
                   nrow=2)
    p.value = fisher.test(data)$p.value
    p.value.list = c(p.value.list, p.value)
  }
  data_raw[color.index+1,] = count.list
  data_test[color.index+1,] = p.value.list
}

df_raw = as.data.frame(data_raw, row.names=1:length(color.set))
df_test = as.data.frame(data_test, row.names=1:length(color.set))
colnames(df_raw) = TFs$TF.family.name
colnames(df_test) = TFs$TF.family.name
#write.csv(df_raw,"combined.module.TF.summary.csv")
write.csv(df_test, "combined.module.TF.fisher.test.csv")


#compare the root and leaf modules
#goal: test whether there is significant overlap between root and leaf co-expression modules
#avoid variable clobbering!
load(file="./round3/Ntab.network.RData")
NtabExpr.root = NtabExpr
bwnet.root = bwnet
load(file="./round4/Ntab.network.RData")
NtabExpr.leaf = NtabExpr
bwnet.leaf = bwnet

total.gene.count.root = length(NtabExpr.root)
total.gene.count.leaf = length(NtabExpr.leaf)
color.set.root = unique(bwnet.root$colors)
color.set.leaf = unique(bwnet.leaf$colors)
data.p.value = matrix(NA, nrow=length(color.set.root), ncol=length(color.set.leaf))
for (root.index in 0:(length(color.set.root)-1)){
  gene.subset.root = colnames(NtabExpr.root[,bwnet.root$colors==root.index])
  p.value.list = c()
  for(leaf.index in 0:(length(color.set.leaf)-1)){
    gene.subset.leaf = colnames(NtabExpr.leaf[,bwnet.leaf$colors==leaf.index])
    overlap = sum(gene.subset.root %in% gene.subset.leaf)
    data = matrix(data = c(overlap, length(gene.subset.leaf)-overlap,length(gene.subset.root)-overlap,
                                     total.gene.count.leaf-length(gene.subset.leaf)-length(gene.subset.root)+overlap), nrow=2)
    p.value.list = c(p.value.list, fisher.test(data)$p.value)
  }
  data.p.value[root.index+1, ]=p.value.list
}
df = as.data.frame(data.p.value, row.names=1:length(color.set.root), col.names=1:length(color.set.leaf))
write.csv(df, 'module.overlap.fisher.test.root.vs.leaf.csv')
