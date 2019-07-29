## Some functions used in GO analysis

### Read a file that contains a list of gene pairs of interest, 
### and return a R list object of these genes of interest
### This function is specifically written for the output of HEB test.
goFromCSV = function(csv, ncol, fdr){
  csv.file = read.csv(file=csv, sep=",", header=TRUE)
  #print(colnames(csv.file))
  glist = csv.file[which(csv.file$q_value <= fdr),][,ncol]
  return (go(glist))
}

goFromFile = function(Infile, col){
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
  return (go(gene.list))
}

### Given a list of genes, return the result of GO Enrichment Analysis.
go = function(gList){
  #print(gList)
  ### GO enrichment analysis
  setwd("C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/GO")
  library(topGO)
  geneID2GO <- readMappings(file = "gene2go.ahrd")
  # the gene universe is given by the following line
  # A total of 34056 genes get one or more GO terms assigned to it
  geneNames = names(geneID2GO)
  gene.of.interest = gList
  gene.of.interest = as.integer(geneNames %in% gene.of.interest)
  gene.of.interest = factor(gene.of.interest)
  names(gene.of.interest) = geneNames # make it intno a named vector
  GOdata = new("topGOdata", ontology="BP", allGenes=gene.of.interest,
               annot=annFUN.gene2GO, gene2GO=geneID2GO)
  test.stat = new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher = getSigGroups(GOdata, test.stat)
  pvals = score(resultFisher)
  allRes = GenTable(GOdata, classic = resultFisher,
                    orderBy = "classic", ranksOf = "classic", topNodes = 20)
  return (allRes)
}