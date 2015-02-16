function [low, high] = combine_confidence_intervals(xlow, xhigh, interval)
% Combine multiple confidence intervals using Gaussian posterior assumption.
%
% [low, high] = combine_confidence_intervals(xlow, xhigh, interval)
%
% ARGUMENTS
%   xlow (1D array) - low end of symmetric confidence intervals about the mean to be combined
%   xhigh (1D array) - high end of symmetric confidence intervals about the mean to be combined
%   interval - fraction of probability encompassed by symmetric confidence intervals about the mean (e.g. 0.95 for 95% CI)
%
% RETURNS
%   low - low end of combined confidence interval
%   high - high end of combined confidence interval

% Determine how many sigma the given confidence interval is equal to.
nsigma = sqrt(2) * erfinv(interval); % 'interval' represents 'nsigma' std devs of gaussian

% Compute means and standard deviations.
mu_n = (xlow + xhigh) / 2.0;
sigma_n = (xhigh - xlow) / 2.0 / nsigma;

% Combine means and standard deviations.
sigma = sum(sigma_n.^(-2))^(-1/2);
mu = sum(sigma_n.^(-2) .* mu_n) / sum(sigma_n.^(-2));

% DEBUG
[mu, sigma]

% Compute new confidence interval.
low = mu - nsigma * sigma;
high = mu + nsigma * sigma;

return

