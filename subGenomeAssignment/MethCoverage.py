# Given a bismark genome-wide cytosine report, return a cdf of cysotine coverage

import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser
import time

parser = OptionParser()
parser.add_option("-r","--report",action="store",type="string",dest="report",
                  help="Path to bismark genome-wide cytosine report.")
(options,args) = parser.parse_args()

CG = []
CHH = []
CHG = []
i = 0
with open(options.report) as r:
    line = r.readline()
    while line:
        print(i)
        i += 1
        chrom, pos, strand, meth, un_meth, epi_context, tri_context = line.strip().split('\t')
        cov = int(meth)+int(un_meth)
        if epi_context == 'CG':
            CG.append(cov)
        elif epi_context == 'CHH':
            CHH.append(cov)
        else:
            CHG.append(cov)
        line = r.readline()


def empirical_cdf(covList):
    covList = np.array(covList)
    bins = list(range(np.max(covList)+1))
    return [sum((covList <= bin).tolist())/len(covList) for bin in bins]
    #for bin in bins:
    #    cdf = sum((covList <= bin).tolist())/len(covList)

# start plotting cdf
CG_cdf = empirical_cdf(CG)
print('CG cdf calculation done')
CHH_cdf = empirical_cdf(CHH)
print('CHG cdf calculation done')
CHG_cdf = empirical_cdf(CHG)
print('CHG cdf calculation done')
MAX = max(len(CG_cdf), len(CHH_cdf), len(CHG_cdf))
bins = np.arange(0, 1 + MAX)
figure = plt.figure()
plt.plot(bins, CG_cdf + [1]*(MAX-len(CG_cdf)),'r', label='CG')
plt.plot(bins, CHH_cdf + [1]*(MAX-len(CHH_cdf)),'m', label='CHH')
plt.plot(bins, CHG_cdf + [1]*(MAX-len(CHG_cdf)),'b', label='CHG')
plt.legend(loc='lower right', fontsize='medium')
plt.xlabel('coverage')
plt.ylabel('cdf')
plt.savefig('MethCoverage.jpg',dpi=300, format='jpg', quality=95)
