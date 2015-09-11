function model = update_transition_matrix(data, model, options)
% Update transition matrix given state trajectories.
%
% model = update_transition_matrix(data, model, options)
%
% ARGUMENTS
%   data (cell array of 1D arrays) - observed trajectories of some real-valued signal
%   model (structure) - model parameters, containing state trajectories
%   options (structure) - options structure
%
% RETURNS
%   model (structure) - model with updated transition matrix
%
% TODO
%   Allow more options to choose the prior for transition matrix sampling.

% Construct transition count matrix from state trajectories.
Nij = zeros(model.nstates, model.nstates, 'int32'); % Nij(i,j) is number of observed transitions to j given system initially in i
for trajectory_index = 1:length(data)
  % Accumulate transition counts from this trajectory.
  Nij = Nij + count_transitions(model.state_trajectories{trajectory_index}, model.nstates);
end

% Update transition matrix using specified method.
switch (options.updateMethod)  
 case 'bayesian'
  if (options.verbosity >= 4)
    disp('Bayesian transition matrix update');
  end
  % Update transition matrix with Bayesian sampling.
  model.Tij = transition_matrix_sample(model.Tij, Nij, options);  
  % Compute stationary probability.
  model.Pi = stationary_probability(model.Tij);
 case 'maximum-likelihood'
  % Determine maximum-likelihood estimator of transition matrix.
  % model.Tij = transition_matrix_mle(Nij, options.tau);  
  error('Not currently implemented');
 otherwise
  error(sprintf('options.updateMethod = "%s" is unknown', options.updateMethod));
end     

return

