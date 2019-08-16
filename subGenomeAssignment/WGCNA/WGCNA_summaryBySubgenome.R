setwd("C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/WGCNA")
options(stringsAsFactors=FALSE)
library(WGCNA)
lnames = load(file="./round2/Ntab.network.RData")

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
  close(con)
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
write.csv(df,"WGCNA_summaryBySubgenome.csv")


# GO term anaylsis for each module
source("C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/GO/goModule.R")
for (color in color.set){
  gene.subset = colnames(NtabExpr[,bwModuleColors == color])
  go.result = go(gene.subset, FALSE)
  write.csv(go.result, paste('C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/WGCNA/round1/GO.network.',color,'.csv',sep=""))
}

# analyze genes in terms of homeologous pairs
# task: how many homeologous genes belong to the same module?
# if not, how further apart they are?

# first load a file specifying homeologous gene pairs and store them in a dictionary
library(hash)
con = file('non.sig.txt','r')
h = hash()
while(TRUE){
  line = readLines(con,n = 1)
  if(length(line) == 0){
    break
  }
  pair = strsplit(trimws(line),"\t")[[1]]
  #print(pair)
  h[[pair[[1]]]] = pair[[2]]
}
close(con)
# query WGCNA object to see whether two genes in a pair belong to the same module
total.pair = 0
sameModule = 0
similarModule = 0
divergentModule = 0
eigenGene.matrix = t(as.matrix(bwMEs))
maxDist = max(dist(eigenGene.matrix, method="euclidean"))

for (Nsyl.g in keys(h)){
  Ntom.g = h[[Nsyl.g]]
  Nsyl.label = bwModuleLabels[Nsyl.g]
  Ntom.label = bwModuleLabels[Ntom.g]
  if (is.na(Nsyl.label) | is.na(Ntom.label)){
    next
  }
  total.pair = total.pair + 1
  if (Nsyl.label != Ntom.label){
    # check how divergent these two modules are
    Nsyl.ME = bwMEs[,paste('ME',Nsyl.label,sep="")]
    Ntom.ME = bwMEs[,paste('ME',Ntom.label,sep="")]
    if (sqrt(sum((Nsyl.ME-Ntom.ME)^2)) <= 0.5*maxDist){
      similarModule = similarModule + 1
    }else{
      divergentModule = divergentModule + 1
    }
  }
  else{
    if (Nsyl.label == 0){
      divergentModule = divergentModule + 1
    }
    else{
      sameModule = sameModule + 1
    }
  }
}

library(ggplot2)
df = data.frame(
  category=c(sprintf("same module(%s)", sameModule),
             sprintf("similar module(%s)", similarModule),
             sprintf("divergent module(%s)", divergentModule)),
  num = c(sameModule, similarModule, divergentModule)
)
bp<- ggplot(df, aes(x="", y=num, fill=category))+
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
# Apply blank theme
library(scales)
pie + scale_fill_grey() +  blank_theme +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(label = percent(num/total.pair)),size=5, 
            position = position_stack(vjust = 0.5))
