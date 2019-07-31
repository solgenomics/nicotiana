# Given a bismark genome-wide cytosine report, return a cdf of cysotine coverage

import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser
import time

parser = OptionParser()
parser.add_option("-r","--report",action="store",type="string",dest="report",
                  help="Path to bismark genome-wide cytosine report.")
(options,args) = parser.parse_args()

CG = {}
CHH = {}
CHG = {}
i = 0

def add(dict, cov):
    if cov in dict:
        dict[cov] += 1
    else:
        dict[cov] = 1

with open(options.report) as r:
    line = r.readline()
    while line:
        print(i)
        i += 1
        chrom, pos, strand, meth, un_meth, epi_context, tri_context = line.strip().split('\t')
        cov = int(meth)+int(un_meth)
        if epi_context == 'CG':
            add(CG, cov)
        elif epi_context == 'CHH':
            add(CHH, cov)
        else:
            add(CHG, cov)
        line = r.readline()

def emp_cdf_fromDict(dict):
    covList = [cov for cov, count in dict.items()]
    TOTAL = sum([count for cov, count in dict.items()])
    bins = list(range(1+np.max(covList)))
    cdf = []
    for bin in bins:
        percent = sum([count for cov, count in dict.items() if cov <= bin])/TOTAL
        cdf.append(percent)
        if percent > 0.8:
            break
    return cdf

def empirical_cdf(covList):
    covList = np.array(covList)
    bins = list(range(np.max(covList)+1))
    return [sum((covList <= bin).tolist())/len(covList) for bin in bins]
    #for bin in bins:
    #    cdf = sum((covList <= bin).tolist())/len(covList)

# start plotting cdf
CG_cdf = emp_cdf_fromDict(CG)
print('CG cdf calculation done')
CHH_cdf = emp_cdf_fromDict(CHH)
print('CHG cdf calculation done')
CHG_cdf = emp_cdf_fromDict(CHG)
print('CHG cdf calculation done')
MAX = max(len(CG_cdf), len(CHH_cdf), len(CHG_cdf))
bins = np.arange(0, MAX)
figure = plt.figure()
plt.xlim(0, MAX)
plt.plot(bins, CG_cdf + [1]*(MAX-len(CG_cdf)),'r', label='CG')
plt.plot(bins, CHH_cdf + [1]*(MAX-len(CHH_cdf)),'m', label='CHH')
plt.plot(bins, CHG_cdf + [1]*(MAX-len(CHG_cdf)),'b', label='CHG')
plt.legend(loc='lower right', fontsize='medium')
plt.xlabel('coverage')

plt.ylabel('cdf')
plt.savefig('MethCoverage.jpg',dpi=300, format='jpg', quality=95)
