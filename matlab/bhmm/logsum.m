function log_sum = logsum(log_a_n)
% Compute log(sum(exp(log_a_n))) in a numerically stable manner.
%
% ARGUMENTS
%  log_a_n (array of N elements) - logarithm of elements to be summed
%
% RETURNS
%  log_sum - log(sum(exp(log_a_n)))

% Determine maximum argument.
max_arg = max(log_a_n);

% Compute sum of shifted arguments.
log_sum = max_arg + log(sum(exp(log_a_n - max_arg)));

return
