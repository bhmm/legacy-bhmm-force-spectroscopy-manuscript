function [low, high] = empirical_confidence_interval(sample, interval)
% Compute specified symmetric confidence interval for empirical sample.
%
% ARGUMENTS
%  sample - empirical sample
%  interval - desired confidence interval (0 < interval < 1) 
%    e.g. 0.68 for 68% confidence interval, 0.95 for 95% confidence interval

% Sort sample in increasing order.
sample = sort(sample);

% Determine sample size.
N = length(sample);

% Compute low and high indices.
low_index = round((N-1) * (0.5 - interval/2)) + 1;
high_index = round((N-1) * (0.5 + interval/2)) + 1;

% Compute low and high.
low = sample(low_index);
high = sample(high_index);

return
