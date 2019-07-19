### Read a file that contains a list of genes of interest, 
### and return a R list object of these genes of interest
geneList = function(Infile, col){
  gene.list = list()
  con = file(Infile, "r")
  while(TRUE){
    line = readLines(con, n=1)
    if (length(line) == 0){
      break
    }
    gene.list = c(gene.list, strsplit(trimws(line),"\t")[[1]][col])
    
  }
  close(con)
  return (gene.list)
}

### GO enrichment analysis
setwd("C:/Users/10453/source/repos/subgenomeAssignment/GO")
library(topGO)
geneID2GO <- readMappings(file = "Ntab.gene2go")
# the gene universe is given by the following line
# A total of 34056 genes get one or more GO terms assigned to it
geneNames = names(geneID2GO)
gene.of.interest = geneList("HEB.Shoot_ZT12.towards.Ntom.txt", 2)
gene.of.interest = as.integer(geneNames %in% gene.of.interest)
gene.of.interest = factor(gene.of.interest)
names(gene.of.interest) = geneNames # make it intno a named vector
GOdata = new("topGOdata", ontology="MF", allGenes=gene.of.interest,
             annot=annFUN.gene2GO, gene2GO=geneID2GO)
test.stat = new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher = getSigGroups(GOdata, test.stat)
pvals = score(resultFisher)
allRes = GenTable(GOdata, classic = resultFisher,
                  orderBy = "classic", ranksOf = "classic", topNodes = 20)
