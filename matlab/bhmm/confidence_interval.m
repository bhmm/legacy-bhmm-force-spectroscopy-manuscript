function [x_low, x_high] = confidence_interval(x, P_low, P_high) 
% CONFIDENCE_INTERVAL
%
% Compute a confidence interval from a sampled dataset.
%
% [x_low, x_high] = confidence_interval(x, P_low, P_high) 

% sort the dataset
xsorted = sort(x);

% get length
N = length(x);

% find appropriate fractions
n_low = floor(N * P_low) + 1;
n_high = floor(N * P_high) + 1;

x_low = xsorted(n_low);
x_high = xsorted(n_high);

return
