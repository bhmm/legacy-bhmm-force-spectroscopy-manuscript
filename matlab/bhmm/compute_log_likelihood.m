function log_likelihood = compute_log_likelihood(model, data, options, use_state_assignments)
% Compute the log-likelihood of the HMM given the data.
%
% log_likelihood = compute_log_likelihood(model, data, options, use_state_assignments)
%
% ARGUMENTS
%   model (structure) - model parameters, containing state trajectories
%   data (cell array of 1D arrays) - observed trajectories of some real-valued signal
%  options - optional options structure
%    if options.reversible is true, then reversible transition matrix sampling will be used 
%    if options.diagonally_dominant is true, then only diagonally-dominant transition matrices will be accepted
%    if options.verbosity sets level of verbose output (0 = none, 1 = summary report, 2 = output each iteration)
%    if options.equilibrium is true, trajectories are initially in equilibrium
%  use_state_assignments - if true, will use state assignments; otherwise, will marginalize them out
%
% RETURNS
%   log_likelihood (double) - log-likelihood of the HMM 'model' given the data 'data', or P(O | \Theta)
%
% NOTE
%   We compute the probability of a trajectory conditioned on the first observation, so that the trajectory may not necessarily be in equilibrium.

% Calculate the stationary probability distribution from dominant eigenvector of transition matrix.
model.logTij = log(model.Tij); % elementwise logarithm - can we use a more numerically stable approach later?

% Initialize log likelihood.
log_likelihood = 0.0;

% Compute log likelihood summed over all trajectories.
for trajectory_index = 1:length(data)
  % Extract trajectory of observations.
  o_t = data{trajectory_index}; % o_t(t) is the observation at time t \in 1...T.

  % Aggregate model parameters into vectors.
  mu = zeros(model.nstates,1);
  sigma = zeros(model.nstates,1);
  for i = 1:model.nstates
    mu(i) = model.states{i}.mu;
    sigma(i) = model.states{i}.sigma;
  end  

  % Accumulate log-likelihood.
  if use_state_assignments
    % Extract hidden state trajectory.
    s_t = model.state_trajectories{trajectory_index};  
  
    log_likelihood = log_likelihood + bhmm_helper.computeTrajectoryLogLikelihood(s_t, o_t, mu, sigma, model.logPi, model.logTij, options.equilibrium);
  else
    log_likelihood = log_likelihood + bhmm_helper.computeModelLogLikelihood(o_t, mu, sigma, model.logPi, model.logTij, options.equilibrium);
  end
end
  
return

