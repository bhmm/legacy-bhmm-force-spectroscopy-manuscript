function model = generate_initial_model(data, nstates, options)
% Generate a (poor) initial model for hidden Markov model (HMM).
%
% model = generate_initial_model(data, nstates, options)
%
% ARGUMENTS
%   data (cell array of 1D arrays) - observed trajectories of some real-valued signal
%   nstates - number of states
%   options - options data structure
%
% RETURNS
%   model (structure) - initial model parameters
%
% NOTES
%   A (poor) guess at the initial model parameters is generated.
%   This guess must be subsequently refined by iterations of either the maximum-likeihood or Bayesian HMM algorithms.
%
% TODO
%   Implement initial guess via EM algorithm
%   Guess number of states by multiple-normal density fitting?

if (options.verbosity >= 1)
  disp('******************************************************************************');
  disp('Bayesian HMM initial model generation');
end
  
% Create initial model from Gaussian mixture model.
model = em_gaussian_mixture(data, nstates, options);

% Generate initial transition matrix.
% Compute expected transition counts using only fractional state assignemnts.
Nij = zeros(nstates, nstates); % Nij(i,j) is fractional expected transition counts from i to j
% Construct transition count matrix from state trajectories.
mu = zeros(nstates, 1);
sigma = zeros(nstates, 1);
Pi = zeros(nstates, 1);
for i = 1:nstates
  mu(i) = model.states{i}.mu;
  sigma(i) = model.states{i}.sigma;
  Pi(i) = model.Pi(i);
end
for trajectory_index = 1:length(data)
  % Extract trajectory of observables.
  o_t = data{trajectory_index};
  T = length(o_t);

  % Compute fractional assignment of samples to states
  p_ti = bhmm_helper.computeStateProbabilities(o_t, mu, sigma, Pi);

  % Accumulate transition counts from this trajectory.
  for t = 1:(T-1)
    Nij = Nij + p_ti(t,:)' * p_ti(t+1,:);
  end
end
% Update transition matrix estimate.
model.Tij = transition_matrix_mle(Nij, options);
clear mu sigma Pi Nij;

% Compute stationary probability.
model.Pi = stationary_probability(model.Tij);
model.logPi = log(model.Pi);

% Assign samples based on probability alone.
%disp('Generating initial guess of state trajectory...');
%mloptions = options;
%%mloptions.updateMethod = 'maximum-likelihood';
%model = update_state_trajectories(data, model, mloptions);

if (options.verbosity >= 1)
  disp('******************************************************************************');
end

return
