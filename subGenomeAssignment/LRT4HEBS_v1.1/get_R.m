% get_R() - Procedure to determine aggregation parameters (r_i) by 
%           creating N auxillary data sets and rescaling the i-th 
%           replicate by D_j/D_i where D_j is the total depth of the 
%           reference data set j, for j = 1,...,N.
%
%           The values r_i, for i=1,...,N, are then determined by 
%           minimizing the sum of squares error from the rescaled data 
%           set to the curve: 
%               
%                   variance = mean + 1/r_i * mean^2
% 
%           for each rescaled data set i.
%
% Usage:
%   >> [r, D] = get_R(data, plot)
%
% Input:
%       data - An array whose rows are genes and columns are total number
%              of mapped reads for each replicate.  Data should include the
%              data for all genes, even if some (e.g., genes who do not have
%              homeologs) will not be tested.
%
%       plot - optional. 1 = Show plots of the data and fitted curve.
%
% Output:
%       r - an Nx1 vector of aggregation parameters, 
%           where N is the number of replicates.
%       
%       D - The total sequencing depth of each replicate.  This is returned
%           for convenience since it needs to be computed along the way.
%
% Author: 
%   Ronald D. Smith
%   Graduate Student, Applied Science
%   The College of William & Mary
%   rdsmith@email.wm.edu
%   April 6, 2017

function [r, D] = get_R(data, plot)
	
    % If [plot] was omitted, do not show plots
    if nargin < 2
        plot = 0;
    end
    
    % Compute the total sequencing depth, in millions.  Rescaling by
    % millions is important for later steps, so that all the numbers will 
    % be approximately the same order of magnitude.
    D = sum(data)/1e6;
    % Determine the number of replicates
    n = length(D);
    % Remove observations that were all zeroes
    data = data(any(data,2), :);
    % Initialize the (inverse of) the output
    phi = nan(1,n);
    
    for i = 1:n
        % Initialize the rescaled data set
        temp = nan(size(data));
        % Rescale the counts as if they had the depth of the i-th replicate
        for j = 1:n
            temp(:, j) = data(:, j)*D(i)/D(j);
        end
        
        % Compute the mean and the variance of the rescaled data set
        m = mean(temp,2); % take the mean by column-wise (i.e, calculate the mean for each re-scaled replicate)
        v = var(temp, [],2); % take variance column-wise
        
        % Fit the curve
        [phi(i)] = FitVarByMean(m,v);
        
        % Optionally show a plot of the data and fitted curve
        if plot
            figure
            loglog(m,v,'.k'); hold on
            r = refline(1,0);
            r.LineWidth = 2;
            loglog(m, m+phi(i)*m.^2, '.r');
            xlabel('Mean of normalized counts')
            ylabel('Variance of normalized counts')
            legend('Data', 'Poisson', 'Best Fit NB', 'location', 'northwest')
        end
    end
% The r's are the inverses of the phi's obtained from the curve fit
r = 1./phi;
end