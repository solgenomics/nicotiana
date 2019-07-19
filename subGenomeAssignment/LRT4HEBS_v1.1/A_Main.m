close all; clc; tic

% This script will perform the likelihood ratio test for changes in 
% homeolog expression bias reported in our paper, "A likelihood ratio test
% for changes in homeolog expression bias", Smith et al, 2017.
%
% In addition to this script, contained in this directory should be:
%
%   1) LRT_NB_HEBS_v8.m: The function that performs the LRT for
%                        changes in HEB
%
%   2) LRT_NH_HEB_v8.m:  The function that performs the LRT for HEB
%
%   3) get_alf.m:        A function for converting the LRT test statistic
%                        to a significance level 
%
%   4) get_W.m:          A function for converting a significance level to
%                        to a test statistic
%
%   5) get_R.m:          A function for determining aggregation parameters
%                        as described in our paper
%
%   6) FitVarByMean.m:   A function used by get_R.m for the curve fitting
%                        procedure
%  
%   7) proc_ShowResults_HEBS.m: A script to plot the results.
%
%   8) W_Mimulus_LeafPetal_3895.mat: The data that was used in the test.
%
%   9) fdr_bh.m:        A function to perform the Benjamini-Hochberg
%                       correction.  
%                       Source: https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr-bh
%                       date: April 6 2017
%
%
%
% Author: 
%   Ronald D. Smith
%   Graduate Student, Applied Science
%   The College of William & Mary
%   rdsmith@email.wm.edu
%   April 6, 2017

%% Load data
% This data contains the count of mapped reads for 5 replicates from the 
% leaf and flower of Mimulus luteus, as well as the total exon length of 
% each gene, and mean RPKM values.  
% 
% The table F contains the flower data; the table L contains the leaf data.  
% These tables are required to compute the total sequencing depths, and to 
% determine the aggregation parameters.
% 
% The tables F_Pars and L_Pars have been filtered to homeologous gene 
% pairs, and they are sorted in the same order. HEB was computed ahead of
% time and is the last column on each table. These tables contain the
% homeolgous pairs that will be tested.
load('W_Mimulus_LeafPetal_3895.mat')

%% Set desired significance level
sig=0.05;

%% Get r's and sequencing depths
[rF, dF] = get_R(F{:,1:5});
[rL, dL] = get_R(L{:,1:5});

%% Initialize output
N = height(F_Pars);
L1 = nan(N,1);
L0 = nan(N,1);

% Commence testing.  Show a waitbar so we know the computer didn't freeze.
h = waitbar(0, 'Doing LRT...');
for idx=1:N
    % Update the waitbar
    waitbar(idx/N, h)
    
    % Get the data for the current homeolog pair from the data tables
    %
    % condition 1 = Flower, condition 2 = Leaf.  
    % a = M. Gutttaus-like homeolog, b = other homeolog
    
    a1 = F_Pars.GutLikeReads(idx, :);
    b1 = F_Pars.OtherReads(idx, :);
    
    a2 = L_Pars.GutLikeReads(idx, :);
    b2 = L_Pars.OtherReads(idx, :);
    
    % Get the total exon length of each gene, in kilobases
    Ka = L_Pars.GutLength(idx)/1e3;
    Kb = L_Pars.OthLength(idx)/1e3;
    
    if any(a1) && any(a2) && any(b1) && any(b2)
        [L1(idx), L0(idx)] = ...
            LRT_NB_HEBS_v8(a1, a2, b1, b2, Ka, Kb, rF, rL, dF, dL);
    end
end

close(h)

%% Compute the test statistics
% Results are in the same order as the input data
W = 2*(L1-L0);
HEBS = L_Pars.HEB - F_Pars.HEB;

%% Make the Benjamini-Hochberg correction for multiple testing
% Exclude un-testable pairs
idx = isfinite(W);
% Convert W's to p-values
p = get_alf(W(idx),1);
% W values greater than Wadj should ensure FDR < [sig]
[~, pcrit]=fdr_bh(p, sig);
Wadj = get_W(pcrit,1);

%% Make the results table 
% Will be sorted in the same order as original input tables
HEBS_Results = table;
HEBS_Results.L0 = L0;
HEBS_Results.L1 = L1;
HEBS_Results.W = W;
HEBS_Results.HEBS = HEBS;

% Show histogram and scatter plot of results
proc_ShowResults_HEBS

% Print total runtime
toc

% Clean up
clearvars a1 a2 b1 b2 D_F D_L h HEBS idx Ka Kb L0 L1 N p pcrit rF rL sig