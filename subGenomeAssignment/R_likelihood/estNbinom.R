S.chrom = 'Nitab15'
T.chrom = 'Nitab08'
S.filename = '../Nsyl.merged.sorted.coverage'
T.filename = '../Ntom.merged.sorted.coverage'

system(sprintf('grep \'%s\' %s > SinS.bed', S.chrom, S.filename), wait=TRUE)
#system(sprintf('grep \'%s\' %s > SinT.bed', S.chrom, T.filename), wait=TRUE)
#system(sprintf('grep \'%s\' %s > TinT.bed', T.chrom, T.filename), wait=TRUE)
#system(sprintf('grep \'%s\' %s > TinS.bed', T.chrom, S.filename), wait=TRUE)

SinS.exp = 58.6
TinT.exp = 78.2
SinT.exp = 0
TinS.exp = 0

logProbSum=fn(p,mean,X){
  return (sum(dnbinom(X,prob=p,mu=mean,log=TRUE)))
}

file.SinS = read.table('SinS.bed', header=FALSE, sep='\t')
counts.SinS = file.SinS[,3]
counts.SinS.sample = sample(counts.SinS,size=100000,replace=FALSE)
optim(par=0.5,fn=logProbSum,mean=SinS.exp,X=counts.SinS.sample, method='CG')

system('rm *.bed', wait=TRUE)