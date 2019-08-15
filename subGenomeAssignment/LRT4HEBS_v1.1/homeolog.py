import matlab.engine
import math
from scipy.stats import chi2
from scipy import mean
import numpy as np
import matplotlib.pyplot as plt
import itertools
from operator import attrgetter
import sys
import gffOp
import csv

class genePair(object):
    def __init__(self, homeolog1, homeolog2, count1, count2, len1, len2, isTestable):
        self.h1 = homeolog1
        self.h2 = homeolog2
        self.count1 = count1
        self.count2 = count2
        self.exonLen1 = len1
        self.exonLen2 = len2
        self.isTestable = isTestable
        self.tpm1 = 0
        self.tpm2 = 0
        self.testStats = 0
        self.pValue = 0
        self.qValue = 0
        self.HEB = 0
        self.isSig = False
        

    def test(self, L1, L0):
        # the test stats calculated here follows chi2 distribution by Wilk's Theorem
        self.testStats = 2*(L1-L0) # L1 and L0 is already log transformed before
        self.pValue = 1 - chi2.cdf(self.testStats, 1) #degree of freedom = 1
        if (math.isnan(self.testStats)):
            print("{} and {}".format(self.h1, self.h2))
            print(self.count1)
            print(self.count2)

        #print(self.pValue)

    #def setHEB(self, b):
    #    self.HEB = b

    def setQvalue(self, q, alpha):
        self.qValue = q
        if q <= alpha:
            self.isSignificant()

    def isSignificant(self):
        self.isSig = True

    def setTPM(self, tpm1, tpm2):
        self.tpm1 = tpm1
        self.tpm2 = tpm2
        self.HEB = math.log2(tpm1/tpm2)
    
    def __eq__(self, other):
        if isinstance(other, genePair):
            return self.h1 == other.h1 and self.h2 == other.h2
        return False

    def __hash__(self):
        return hash((self.h1, self.h2))


        

class experiment(object):
    def __init__(self, aggParam, numReads, genePairDict, gc_dict, exonDict, name):
        self.eng = matlab.engine.start_matlab()
        self.name = name
        print(f"Start Processing Experiment {self.name}")
        self.aggParam = aggParam
        self.numReads = numReads
        self.gPairList = []

        # calculate the scaling factor for TPM for each replicate
        exonList = []
        countList = []
        for g,count in gc_dict.items():
            exonList.append(exonDict[g])
            countList.append(count)
        T = np.sum((np.array(countList)).T/np.array(exonList), axis=1)[:,np.newaxis] # make it into a column vector
        TPM = np.exp(np.log(np.array(countList).T+1)+np.log(1e6)-np.log(exonList)-np.log(T))
        
        for g1,g2 in genePairDict.items():
            count1, count2 = gc_dict[g1], gc_dict[g2]
            isTestable = len([c for c in count1 if c > 0])/len(count1) >= 0.5 and len([c for c in count2 if c > 0])/len(count2) >= 0.5
            #if len([c for c in count1 if c > 0])/len(count1) >= 0.5 and len([c for c in count2 if c > 0])/len(count2) >= 0.5:
            #isTestable = mean(count1) >= 10 and mean(count2) >= 10
            exonLen1, exonLen2 = exonDict[g1], exonDict[g2]
            gp = genePair(g1, g2, count1, count2, exonLen1, exonLen2, isTestable)
            self.__addHomeologousGene(gp, T)

    def __addHomeologousGene(self, gPair, T):
        self.gPairList.append(gPair)
        if gPair.isTestable:
            L1,L0,v_H1,y_H1,v_H0,stats = self.eng.LRT_NB_HEB_v8(matlab.int64(gPair.count1), matlab.int64(gPair.count2), 
                                                       gPair.exonLen1, gPair.exonLen2,
                                                       self.aggParam, self.aggParam,
                                                       self.numReads, 0, nargout=6)
        
            # Does the following calculation make sense for paired-end data?
            # Can we do better with RPKM?
            # Is this still comparable across samples?
            tpm1 = np.mean(np.exp(np.log(np.array(gPair.count1)+1)+np.log(1e6)-np.log(gPair.exonLen1)-np.log(T)))
            tpm2 = np.mean(np.exp(np.log(np.array(gPair.count2)+1)+np.log(1e6)-np.log(gPair.exonLen2)-np.log(T)))
            #rpkm1 = np.mean(np.array(gPair.count1)/((gPair.exonLen1/1000)*np.array(self.numReads)))
            #rpkm2 = np.mean(np.array(gPair.count2)/((gPair.exonLen2/1000)*np.array(self.numReads)))
            #gPair.setHEB(math.log(rpkm1/rpkm2,2)) # use log base 2
            #gPair.test(L1, L0)
            gPair.setTPM(tpm1, tpm2)
            #gPair.setHEB(math.log(tpm1/tpm2))
            gPair.test(L1, L0)

    def fdr(self, alpha):
        # performs Benjamini-Hochberg multiple testing correction
        testList = [g for g in self.gPairList if g.isTestable]
        testList.sort(key=lambda gp:gp.pValue)
        p_sorted = [g.pValue for g in testList]
        bias_all = [g.HEB for g in testList] # store HEB for all homeologous gene pairs, used later for plotting
        bias_sig = [] # store HEB for all homeologous gene paris that show statistically significant biased expression, used for plotting later

        numGene = len(testList)
        m = np.arange(1,numGene+1,1)
        q_vals = numGene*np.array(p_sorted)/m
        p_cut = max(list(itertools.compress(p_sorted, (q_vals < alpha).tolist())))
        bias_sig = [g.HEB for g in testList if g.pValue <= p_cut]

        # write test result to csv
        with open(f'HEB.{self.name}.towards.Nsyl.csv','w',newline='') as Nsyl_csv:
            with open(f'HEB.{self.name}.towards.Ntom.csv','w',newline='') as Ntom_csv:
                Nsyl_writer = csv.writer(Nsyl_csv, delimiter=',')
                Ntom_writer = csv.writer(Ntom_csv, delimiter=',')
                Nsyl_writer.writerow(['Nsyl','Ntom','HEB','p_value','q_value'])
                Ntom_writer.writerow(['Nsyl','Ntom','HEB','p_value','q_value'])
                index = 0
                for gp in testList:
                    gp.setQvalue(q_vals[index], alpha) #because testList is sorted by qValue, so this is fine
                    index += 1
                    if gp.HEB > 0:
                        Nsyl_writer.writerow([gp.h1, gp.h2, gp.HEB, gp.pValue, gp.qValue])
                    elif gp.HEB < 0:
                        Ntom_writer.writerow([gp.h1, gp.h2, gp.HEB, gp.pValue, gp.qValue])


        #with open(f'HEB.{self.name}.towards Nsyl.txt','w') as out1:
        #    with open(f'HEB.{self.name}.towards Ntom.txt','w') as out2:
        #        i = 0
        #        p = p_sorted[i]
        #        while p <= p_cut:
        #            if testList[i].HEB > 0:
        #                out1.write('{}\t{}\n'.format(testList[i].h1, testList[i].h2))
        #            else:
        #                out2.write('{}\t{}\n'.format(testList[i].h1, testList[i].h2))

        #            bias_sig.append(testList[i].HEB)
        #            testList[i].isSignificant()
        #            i += 1
        #            p = p_sorted[i]
            
        self.__plot(bias_all, bias_sig)
        gffOp.drawScatterHEB(self.gPairList)

    def __plot(self, bias_all, bias_sig):
        fig = plt.figure()
        bin_all = np.arange(-15,15,0.25)
        bin_pos = np.arange(0,15,0.25)
        bin_neg = np.arange(-15,0,0.25)
        bias_sig_pos = [b for b in bias_sig if b > 0]
        bias_sig_neg = [b for b in bias_sig if b < 0]
        plt.hist(bias_all, bins=bin_all, color='#A9A9A9', label='Overall n={}'.format(len(bias_all)))
        plt.hist(bias_sig_pos, bins=bin_pos, color='#0000CD', label='Towards Nsyl, $n_0$={}'.format(len(bias_sig_pos)))
        plt.hist(bias_sig_neg, bins=bin_neg, color='#FF0000', label='Towards Ntom, $n_1$={}'.format(len(bias_sig_neg)))
        plt.xlabel('Homeolog Expression Bias $\log_2 \frac{r^{Nsyl}}{r^{Ntom}}$')
        plt.title(f'Homeologous Gene Expression Bias in N.tabacum {self.name}')
        plt.legend(loc='upper right', fontsize='small')
        plt.savefig(f'Homeologous Gene Expression Bias in N.tabacum {self.name}.jpg', dpi=300, format='jpg', quality=95)

       


class gpTwo():
    # gp1 and gp2 must be the same homeologous gene pair coming from two different experiment conditions
    def __init__(self, gp1, gp2, L1, L0):
        if gp1.h1 != gp2.h1 or gp1.h2 != gp2.h2:
            raise ValueError("Must give the same gene pair under different experiments")
        self.gp1 = gp1
        self.gp2 = gp2
        self.HEB1 = gp1.HEB
        self.HEB2 = gp2.HEB
        self.delta_HEB = gp2.HEB-gp1.HEB
        self.testStats = 2*(L1-L0)
        self.pValue =  1 - chi2.cdf(self.testStats, 1)
        #self.isSig = False


# Take in two experiment objects, and see how many gene pairs show significant changes in HEB.
def compExpr(expr1, expr2):
    # First we get a list of genes that are both testable in expr1 and expr2
    # expr1.gPairList and expr2.gPairList should have gene pairs in the same order (and they would have the same number of gene pairs)
    if len(expr1.gPairList) != len(expr2.gPairList):
        raise ValueError("Invalid comparison between experiments")
    
    print(f'Starting Comparing Change in Expression Bias between {expr1.name} and {expr2.name}')
    gpl = [(expr1.gPairList[i], expr2.gPairList[i]) 
           for i in range(0,len(expr1.gPairList)) 
           if expr1.gPairList[i].isTestable and expr2.gPairList[i].isTestable]
    gpTwoL = []
    eng = matlab.engine.start_matlab()
    for gp in gpl:
        gp_in_expr1, gp_in_expr2 = gp
        L1, L0, v1_H1, v2_H1, y1_H1, y2_H1, v1_H0, v2_H0, y_H0, stats = eng.LRT_NB_HEBS_v8(
            matlab.int64(gp_in_expr1.count1), matlab.int64(gp_in_expr2.count1), 
            matlab.int64(gp_in_expr1.count2), matlab.int64(gp_in_expr2.count2), 
            gp_in_expr1.exonLen1, gp_in_expr2.exonLen2,
            expr1.aggParam, expr2.aggParam, expr1.numReads, expr2.numReads, nargout=10)
        gpTwoL.append(gpTwo(gp_in_expr1, gp_in_expr2, L1, L0))

    print('Start Benjamini-Hochberg Procedure...')
    gpTwoL.sort(key=attrgetter('pValue'))
    numGp = len(gpTwoL)
    m = np.arange(1,numGp+1,1)
    p_sorted = [gptwo.pValue for gptwo in gpTwoL]
    q_vals = numGp*np.array(p_sorted)/m

    valid_p = list(itertools.compress(p_sorted, (q_vals < 0.05).tolist()))
    if len(valid_p) == 0:
        print('No gene has statistically significant change in bias expression at FDR=0.05')
        print('A plot of p-values are provided for reference.')
        gffOp.draw_scatter(gpTwoL, -1, expr1.name, expr2.name)
        sys.exit()

    p_cut = max(valid_p)
    with open(f'Change in HEB towards Nsyl from {expr1.name} to {expr2.name}.txt','w') as out1:
        with open(f'Change in HEB towards Ntom from {expr1.name} to {expr2.name}.txt','w') as out2:
            for gptwo in gpTwoL:
                if gptwo.pValue <= p_cut: 
                    if gptwo.delta_HEB > 0:
                        out1.write(f'{gptwo.gp1.h1}\t{gptwo.gp1.h2}\n')
                    else:
                        out2.write(f'{gptwo.gp1.h1}\t{gptwo.gp1.h2}\n')
                else:
                    break


    print('Calculation Done. Start Plotting...')
    gffOp.draw_hist_changeHEB(gpTwoL, p_cut, expr1.name, expr2.name)
    gffOp.draw_scatter(gpTwoL, p_cut, expr1.name, expr2.name)
    




    




