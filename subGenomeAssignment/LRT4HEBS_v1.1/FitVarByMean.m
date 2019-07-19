% Function to find the value b that best fits the data to:
%   variance = mean + b * mean^2
%
% Used to determine aggregation parametrs r = 1/b for a data set under the
% assumption of a negative binomial distribution.
%
% Input:
%
%       means - a list of observed means
%       vars  - a list of observed variances
%
% Output:
%
%       b - the value of b that provided the best fit
%       fval - the error at each point
%       r2 - the coefficient of determination
%
%  Note: Data points with mean=variance=0 should be filtered from the data
%  prior to calling this function.
%
% Author: 
%   Ronald D. Smith
%   Graduate Student, Applied Science
%   The College of William & Mary
%   rdsmith@email.wm.edu
%   April 6, 2017

function [b, fval, r2] = FitVarByMean(means, vars)

global x y 
    opts = optimset('Display', 'off');
    x = means;
    y = vars;
    
    b = 1;
    b = lsqnonlin(@objfun,b,0,inf,opts);
    fval = feval(@objfun, b);
    
    r2 = 1 - sum( fval.^2 ) / sum( (y-mean(y)).^2 );
end

function F = objfun(params)
    global x y
    b = params(1);
    F = x + b*x.^2 - y;
end