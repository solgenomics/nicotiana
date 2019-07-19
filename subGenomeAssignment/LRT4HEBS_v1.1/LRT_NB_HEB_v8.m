% LRT_NB_HEB_v8() - Performs the likelihood ratio test for homeolog
%                   expression bias (Smith et al, 2017). Given two genes in
%                   a single condition, the total sequencing depth of the
%                   RNA-experiments, and the length of the coding region of
%                   the two genes, computes the log-likelihood of the null
%                   hypothesis, that the expression levels are the same,
%                   and the alternative hypothesis, that the expression
%                   levels are different. These values can be used to 
%                   compute the test statistic W = 2(L1-L0) which, by 
%                   Wilks' theorem, can be compared to a chi-squared 
%                   distribution to determine whether there is sufficient 
%                   evidence to reject the null at some significance level.
%
% Usage:
%   >> [L1, L0, v_H1, y_H1, v_H0, stats] = ...
%         LRT_NB_HEB_v8(a_data, b_data, Ka, Kb, Ra, Rb, D, skip_null)
%
%
% Required inputs:
%
%       a_data - a row vector of length N containing the number of mapped
%                reads for gene A
%
%       b_data - a row vector of length N containing the number of mapped
%                reads for gene B
%
%       Ka - the length of the coding region of gene A
%
%       Kb - the length of the coding region of gene B
%
%       Ra, Rb - should be equal for this implementation.  A row vector
%                of length N containing the aggregation parameters for 
%                each replicate
%
%       D - a row vector of length N containing the number of mapped reads
%           for each replicate
%       
%       skip_null - should be 0 if testing HEB.  When testing HEBS, this
%                   procedure is called but does not require any 
%                   calculation of the null hypothesis
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
%       v_H1 - the value obtained during maximum likelihood estimation for
%              the parameter v under the alternative hypothesis
%
%       y_H1 - the value obtained during maximum likelihood estimation for
%              the parameter y under the alternative hypothesis
%       
%       v_H0 - the value obtained during maximum likelihood estimation for
%              the parameter v under the null hypothesis
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

function [L1, L0, v_H1, y_H1, v_H0, stats] = ...
        LRT_NB_HEB_v8(a_data, b_data, Ka, Kb, Ra, Rb, D, skip_null)
    
    if nargin < 8
        skip_null = 0;
    end

    opts = optimoptions('fsolve', 'Display', 'off', ...
                        'FunctionTolerance', 1e-12, ...
                        'MaxIterations', 1000, ...
                        'MaxFunctionEvaluations', 6000, ...
                        'FiniteDifferenceType', 'central');
    
    global a b ka kb ra rb d n sa sb hyp s
    a = double(a_data);
    b = double(b_data);
    ka = double(Ka);
    kb = double(Kb);
    ra = Ra;
    rb = Rb;
    d = D;
    n = length(d); % Add error handling for unequal vector lengths.
    sa = sum(a);
    sb = sum(b);
    s = sa+sb;
    init_v = log(mean([mean(a/ka./d) mean(b/kb./d)]));
    init_y = 0;
    
    % H1
    hyp = 1;
    [params, fval_H1, flag_H1] = fsolve(@objfun, [init_v init_y], opts);
    
    v_H1 = params(1);
    y_H1 = params(2);
    
    if ~skip_null
        % H0
        hyp = 0;
        [params, fval_H0, flag_H0] = fsolve(@objfun, init_v, opts);
        %if TWICE
        %    [params, fval_H0, flag_H0] = fsolve(@objfun, params, opts);
        %end
        v_H0 = params(1);
        stats = [flag_H0 fval_H0 flag_H1 fval_H1];
        
        L1 = likely(v_H1, y_H1);
        L0 = likely(v_H0, 0);
    else
        v_H0 = nan;
        stats = [flag_H1 fval_H1];
        
        L1 = nan;
        L0 = nan;
    end
end

function F = objfun(params)
    global a b ka kb ra rb d sa sb hyp
    
    v = params(1);
    
    if hyp==1
        y = params(2);
    else
        y = 0;
    end
    
    A = sa - sum( (a+ra).*(ka*d)./(ka*d+ra*exp(-v)*exp(y)) );
    B = sb - sum( (b+rb).*(kb*d)./(kb*d+rb*exp(-v)*exp(-y)) );
    
    if hyp == 1
        F(2) = -A+B;
    end
    F(1) = A+B;
end

function L = likely(v,y)
    global a b ka kb ra rb d
    
    L_matters4test = sum( a.*log(exp(v-y)*ka*d) ) - sum( (a+ra).*log(exp(v-y)*ka*d+ra) ) + ...
                     sum( b.*log(exp(v+y)*kb*d) ) - sum( (b+rb).*log(exp(v+y)*kb*d+rb) );
    
    L = L_matters4test;
    
end