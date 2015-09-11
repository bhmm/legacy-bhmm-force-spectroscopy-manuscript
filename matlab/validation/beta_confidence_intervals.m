function [Plow, Phigh] = beta_confidence_intervals(ci_X, ntrials, ci)
% Compute confidence interval of beta distributions.
%
% ARGUMENTS
%   ci_X (vector) - computed confidence interval estimate from 'ntrials' experiments
%   ntrials (int) - number of trials
%   ci (float) - confidence interval to report for (e.g. 0.95 for 95% confidence interval)

% Compute low and high confidence interval for symmetric CI about mean.
ci_low = 0.5 - ci/2;
ci_high = 0.5 + ci/2;

% Compute for every element of ci_X.
Plow = ci_X * 0.0;
Phigh = ci_X * 0.0;
for i = 1:size(ci_X,1)
  for j = 1:size(ci_X,2)
    X = betainv([ci_low, ci_high], ci_X(i,j) * ntrials, (1-ci_X(i,j)) * ntrials);
    Plow(i,j) = X(1);
    Phigh(i,j) = X(2);
  end
end

return 
