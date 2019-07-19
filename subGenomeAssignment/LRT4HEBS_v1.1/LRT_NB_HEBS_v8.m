% LRT_NB_HEBS_v8() - Performs the likelihood ratio test for changes in
%                    homeolog expression bias (Smith et al, 2017). Given
%                    two genes in two conditions, the total sequencing
%                    depths of the RNA-seq experiments, and the length of
%                    the coding region of the two genes, computes the
%                    log-likelihoods of the null hypothesis, that the bias
%                    is the same in both conditions, and the alternative
%                    hypothesis that the bias is different.  These values
%                    can be used to compute the test statistic W = 2(L1-L0) 
%                    which, by Wilks' theorem, can be compared to a
%                    chi-squared distribution to determine whether there is
%                    sufficient evidence to reject the null at some
%                    significance level.
%
% Usage:
% 
% >> [L1, L0, v1_H1, v2_H1, y1_H1, y2_H1, v1_H0, v2_H0, y_H0, stats] = ...
%       LRT_NB_HEBS_v8(a1_data, a2_data, b1_data, b2_data, ...
%                      Ka, Kb, R1, R2, D1, D2)
%
% Required inputs:
%       a1_data - a row vector of length N containing the number of 
%                 mapped reads for gene A in condition 1
%
%       a2_data - a row vector of length N containing the number of 
%                 mapped reads for gene A in condition 2
%
%       b1_data - a row vector of length N containing the number of 
%                 mapped reads for gene B in condition 1
%
%       b2_data - a row vector of length N containing the number of 
%                 mapped reads for gene B in condition 2
%
%       Ka - the length of the coding region of gene A, in kilobases
%
%       Kb - the length of the coding region of gene B, in kilobases
%
%       R1 - a row vector of length N containing the aggregation parameters
%            for each replicate in condition 1
%
%       R2 - a row vector of length N containing the aggregation parameters
%            for each replicate in condition 2
%
%       D1 - a row vector of length N containing the number of mapped
%            reads, in millions, for each replicate in condition 1
%
%       D2 - a row vector of length N containing the number of mapped
%            reads, in millions, for each replicate in condition 2
%
% Outputs:
%
%       L1 - the log-likelihood of the alternative hypothesis, that bias is
%            different between the two conditions, excluding constant terms
%            that cancel in the calculation of the test statistic 
%            W = 2(L1-L0)
%
%       L0 - the log-likelihood of the null hypothesis, that bias is not
%            different between the two conditions, excluding constant terms
%            that cancel in the calculation of the test statistic 
%            W = 2(L1-L0)
%
%   The rest of the outputs are included for debugging:
%
%       v1_H1 - the value obtained during maximum likelihood estimation for
%               the parameter v1 under the alternative hypothesis
%
%       v2_H1 - the value obtained during maximum likelihood estimation for
%               the parameter v2 under the altnernative hypothesis
%
%       y1_H1 - the value obtained during maximum likelihood estimation for
%               the parameter y1 under the alternative hypothesis
%
%       y2_H1 - the value obtained during maximum likelihood estimation for
%               the parameter y2 under the alternative hypothesis
%       
%       v1_H0 - the value obtained during maximum likelihood estimation for
%               the parameter v1 under the null hypothesis
%       
%       v2_H0 - the value obtained during maximum likelihood estimation for
%               the parameter v2 under the null hypothesis
%
%       y_H0  - the value obtained during maximum likelihood estimation for
%               the parameter y under the null hypothesis
%
%       stats - for debugging. Includes flags returned by
%               fsolve and the value of the objective functions
%               evaluated at their most likely parameters.
%
%
% Author: 
%   Ronald D. Smith
%   Graduate Student, Applied Science
%   The College of William & Mary
%   rdsmith@email.wm.edu
%   April 6, 2017

function [L1, L0, v1_H1, v2_H1, y1_H1, y2_H1, v1_H0, v2_H0, y_H0, stats] = ...
            LRT_NB_HEBS_v8(a1_data, a2_data, b1_data, b2_data, Ka, Kb, R1, R2, D1, D2)
       
    opts_f = optimoptions('fsolve', 'Display', 'off', ...
                          'FunctionTolerance', 1e-12, ...
                          'MaxIterations', 5000, ...
                          'MaxFunctionEvaluations', 40000, ...
                          'FiniteDifferenceType', 'central');
    global a1 b1 a2 b2 ka kb r1 r2 d1 d2 n sa1 sa2 sb1 sb2 sa sb s1 s2
    a1 = double(a1_data);
    b1 = double(b1_data);
    a2 = double(a2_data);
    b2 = double(b2_data);
    
    ka = double(Ka);
    kb = double(Kb);
    
    r1 = R1;
    r2 = R2;
    
    d1 = D1;
    d2 = D2;
    
    n = length(d1);
    sa1 = sum(a1);
    sa2 = sum(a2);
    sb1 = sum(b1);
    sb2 = sum(b2);
    sa = sa1 + sa2;
    sb = sb1 + sb2;
    s1 = sa1 + sb1;
    s2 = sa2 + sb2;
    
    % H1
    [~, ~, v1_H1, y1_H1, ~, stats_1] = LRT_NB_HEB_v8(a1, b1, ka, kb, r1, r1, d1, 1);
    [~, ~, v2_H1, y2_H1, ~, stats_2] = LRT_NB_HEB_v8(a2, b2, ka, kb, r2, r2, d2, 1);
    
    flag_H1 = [stats_1(1) stats_2(1)];
    fval_H1 = [stats_1(2:3) stats_2(2:3)];
    
    % H0
    [params, fval_H0, flag_H0] = fsolve(@objfun, [v1_H1 v2_H1 0], opts_f);
    
    v1_H0 = params(1);
    v2_H0 = params(2);
    y_H0 = params(3);
    
    stats = [flag_H0 fval_H0 flag_H1 fval_H1];
    
    L1 = likely(v1_H1, v2_H1, y1_H1, y2_H1);
    L0 = likely(v1_H0, v2_H0, y_H0, y_H0); 
end

function F = objfun(params)
    global a1 b1 a2 b2 ka kb r1 r2 d1 d2 sa1 sa2 sb1 sb2
    % This function only needs to handle null hypothesis.
    %
    % The alternative hypothesis is two independent copies of the
    % alternative hypothesis for HEB.
    
    v1 = params(1);
    v2 = params(2);
    y = params(3);
    
    A1 = sa1 + sum(-(a1+r1).*(ka*d1)./(ka*d1+r1*exp(-v1)*exp(+y)) );
    A2 = sa2 + sum(-(a2+r2).*(ka*d2)./(ka*d2+r2*exp(-v2)*exp(+y)) );
    B1 = sb1 + sum(-(b1+r1).*(kb*d1)./(kb*d1+r1*exp(-v1)*exp(-y)) );
    B2 = sb2 + sum(-(b2+r2).*(kb*d2)./(kb*d2+r2*exp(-v2)*exp(-y)) );
    
    F(1) = A1 + B1;
    F(2) = A2 + B2;
    F(3) = -A1 - A2 + B1 + B2;
end

function L = likely(v1,v2,y1,y2)
    global a1 b1 a2 b2 ka kb r1 r2 d1 d2
    % Terms that cancel when computing W = 2*(L1-L0) are not computed.
    L_matters4test = sum( a1.*log(exp(v1-y1)*ka*d1) ) - sum( (a1+r1).*log(exp(v1-y1)*ka*d1+r1) ) + ...
                     sum( b1.*log(exp(v1+y1)*kb*d1) ) - sum( (b1+r1).*log(exp(v1+y1)*kb*d1+r1) ) + ...
                     sum( a2.*log(exp(v2-y2)*ka*d2) ) - sum( (a2+r2).*log(exp(v2-y2)*ka*d2+r2) ) + ...
                     sum( b2.*log(exp(v2+y2)*kb*d2) ) - sum( (b2+r2).*log(exp(v2+y2)*kb*d2+r2) );
    L = L_matters4test;
end