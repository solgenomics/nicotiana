% get_alf() - Function to find the significance level [alf] that would 
%             produce a critical value [W_CRIT] from a chi-squared 
%             distribution with [deg] degrees of freedom, for a 
%             one-tailed test. This procedure is used to convert W values 
%             from a likelihood ratio test into p-values for the purposes 
%             of multiple-testing correction.
%
% Usage: 
%   >> alf = get_Alf(W_CRIT, deg)
%
% Input:
%
%       W_CRIT - The value of W we wish to find [alf] for
%       
%       deg - The degrees of freedom between the null and alternative
%             hypotheses.  If omitted, deg = 1.
%
% Output:
%
%       alf - The significance level [alf] that yields the critical value
%             W_CRIT
%
% Author: 
%   Ronald D. Smith
%   Graduate Student, Applied Science
%   The College of William & Mary
%   rdsmith@email.wm.edu
%   April 6, 2017

function alf = get_alf(W_CRIT, deg)
    if nargin < 2
        deg = 1;
    end
    alf = 1 - chi2cdf((W_CRIT),deg);
end