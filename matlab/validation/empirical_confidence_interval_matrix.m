function [low, high] = empirical_confidence_interval_matrix(sample, interval)
% Compute specified symmetric confidence interval for empirical sample.
%
% ARGUMENTS
%  sample - empirical sample (could be a matrix)
%  interval - desired confidence interval (0 < interval < 1) 
%    e.g. 0.68 for 68% confidence interval, 0.95 for 95% confidence interval

shape = size(sample);
ndim = length(shape) - 1;
nelements = prod(shape);
N = shape(1);
low = zeros(shape(2:end));
high = zeros(shape(2:end));

for i = 1:nelements
  % Collect elements.
  if ndim == 1
    x_n = sample(n,:);
  elseif ndim == 2
    x_n = sample(n,:,:);
  else
    error('not supported');
  end
  
  % Sort sample in increasing order.
  sorted = sort(x_n);

  % Determine sample size.
  N = length(sorted);
  
  % Compute low and high indices.
  low_index = round((N-1) * (0.5 - interval/2)) + 1;
  high_index = round((N-1) * (0.5 + interval/2)) + 1;
  
  % Compute low and high.
  low(i) = sorted(low_index);
  high(i) = sorted(high_index);
end
  
return
