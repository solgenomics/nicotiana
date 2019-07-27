## DE analysis using DESeq2
## input data format: htseq-count output
source("C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/GO/goModule.R")

fdr = 0.05
dir.htseq = "C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/DEseq2/htseq/Control vs meJA(NC95)"
files.htseq = grep("htseq",list.files(dir.htseq),value=TRUE)
time = c("T0","T0","T0","T2","T2","T2","T2","T2","T2",
         "T6","T6","T6","T6","T6",
         "T24","T24","T24","T24","T24","T24")
treated = c("N","N","N","N","N","N","Y","Y","Y",
            "N","N","N","Y","Y","N","N","N","Y","Y","Y")
sampleTable <- data.frame(sampleName = files.htseq,
                          fileName = files.htseq,
                          time=time, treat=treated)
library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = dir.htseq,
                                       design= ~ time+treat+time:treat)
