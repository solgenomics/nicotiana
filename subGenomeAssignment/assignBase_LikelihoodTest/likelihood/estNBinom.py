import subprocess
import numpy as np
from scipy.stats import nbinom
from scipy.optimize import minimize

S_chrom = 'Nitab15'
T_chrom = 'Nitab08'




class nbinomModel(object):
    
    def __init__(self, SbedFile, TbedFile):
        subprocess.call('grep \'{}\' {} > SinS.bed'.format(S_chrom, SbedFile), shell=True)
        subprocess.call('grep \'{}\' {} > TinT.bed'.format(T_chrom, TbedFile), shell=True)
        subprocess.call('grep \'{}\' {} > SinT.bed'.format(S_chrom, TbedFile), shell=True)
        subprocess.call('grep \'{}\' {} > TinS.bed'.format(T_chrom, SbedFile), shell=True)
        SinS = open('SinS.bed','r')
        TinT = open('TinT.bed','r')
        SinT = open('SinT.bed','r')
        TinS = open('TinS.bed','r')
        SinSline = SinS.readline()
        TinTline = TinT.readline()
        SinTline = SinT.readline()
        TinSline = TinS.readline()
        S_countsInS = []
        T_countsInT = []
        S_countsInT = []
        T_countsInS = []
        
        while SinSline:
            content = SinSline.strip().split('\t')
            S_countsInS.append(int(content[2]))
            SinSline = SinS.readline()

        while TinTline:
            content = TinTline.strip().split('\t')
            T_countsInT.append(int(content[2]))
            TinTline = TinT.readline()

        while SinTline:
            content = SinTline.strip().split('\t')
            S_countsInT.append(int(content[2]))
            SinTline = SinT.readline()

        while TinSline:
            content = TinSline.strip().split('\t')
            T_countsInS.append(int(content[2]))
            TinSline = TinS.readline()

        SinS.close()
        TinT.close()
        SinT.close()
        TinS.close()

        
        
    def logProbSum(countArray, n, p):
        array = np.array(countArray)
        return sum(np.log(nbinom(array, n, p)))



        







