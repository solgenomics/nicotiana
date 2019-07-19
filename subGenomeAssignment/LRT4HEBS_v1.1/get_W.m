% get_W() - Function to find a critical value [W_CRIT] corresponding to a
%           significance level [alf] for a chi-squared distribution with 
%           [deg] degrees of freedom.  Assumes a one-tailed test.
%
%           This procedure is used to find the critical value for a 
%           likelihood ratio test of composite hypotheses at a desired 
%           significance level.
%
% Usage:
%   >> W_CRIT = get_W(alf, deg, prec)
%
% Input:
%
%       alf - The desired significance level
%
%       deg - The degrees of freedom between the null and alternative
%       hypotheses.
%   
%       prec - The desired decimal precision.  If omitted, 10 decimal
%              places are computed.
%
% Author: 
%   Ronald D. Smith
%   Graduate Student, Applied Science
%   The College of William & Mary
%   rdsmith@email.wm.edu
%   April 6, 2017

function W_CRIT = get_W(alf, deg, prec)
    if nargin < 3
        prec = 10;
    end
    w=0:exp(-prec):100;
    c = chi2cdf(w,deg);
    a = find(c>=(1-alf),1,'first');
    W_CRIT = w(a);
end