function log_likelihood = compute_gmm_log_likelihood(mu, sigma, logPi, data)
% Compute the log-likelihood of the GMM given the data.
%
% log_likelihood = compute_gmm_log_likelihood(model, data)
%
% ARGUMENTS
%   model (structure) - model parameters, containing state trajectories % DIFFERENT
%   data (cell array of 1D arrays) - observed trajectories of some real-valued signal
%
% RETURNS
%   log_likelihood (double) - log-likelihood of the GMM 'model' given the data 'data', or P(O | \Theta)

% Initialize log likelihood.
log_likelihood = 0.0;

% Compute log likelihood summed over all trajectories.
for trajectory_index = 1:length(data)
  % Extract trajectory of observations.
  o_t = data{trajectory_index}; % o_t(t) is the observation at time t \in 1...T.

  % Accumulate log-likelihood.
  log_likelihood = log_likelihood + bhmm_helper.computeGmmLogLikelihood(o_t, mu, sigma, logPi);
end
  
return

