library(karyoploteR)

header.lines = readLines('Ntab.v5.chrom.sizes.13-24')
ll = header.lines[grepl(header.lines, pattern = 'Nitab')]
gg = data.frame(do.call(rbind, strsplit(ll, split = '\t')))
gg[,2] = as.numeric(as.character(gg[,2]))
gg[,3] = as.numeric(as.character(gg[,3]))
Ntab.genome = toGRanges(gg)
kp <- plotKaryotype(genome=Ntab.genome,ideogram.plotter=NULL,main='N.tabacum Gene Density')

library(rtracklayer)
gff.nitab = import('Nitab.13-24.gff')
genes = gff.nitab[gff.nitab$type == 'gene']
kpPlotDensity(kp, genes, data.panel='ideogram',window.size=0.5e6,col="#3388FF",border="#3388FF")
