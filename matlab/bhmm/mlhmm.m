function mlmodel = mlhmm(data, nstates, options, initial_model)
% Construct maximum-likelihod hidden Markov model (HMM), useful for initializing a Bayesian HMM.
%
% mlmodel = mlhmm(data, nstates, options, [initial_model])
%
% ARGUMENTS
%   data (cell array of 1D arrays) - observed trajectories of some real-valued signal (e.g. trap extension, FRET efficiency)
%   nstates (int) - number of states to fit
%   options - options structure to control model generation
%
% OPTIONAL ARGUMENTS
%   initial_model - model to initialized MLHMM procedure
%
% RETURNS
%   mlmodel (structure) - maximum likelihood model after fit to data
%
% NOTES
%   Trajectories provided in 'data' cell array may be of different lengths. 
%   The expectation-maximization (EM) algorithm is used for updates; this is a variation on Baum-Welch.
%   This may not produce the global maximum-likelihood model.
%
% TODO
%   Add ability for user to pass in optional 'options' structure.

% Import fast Java helper code.
import bhmm_helper.*;

% DEBUG
RELATIVE_TOLERANCE = options.convergenceTolerance;

if (nargin == 3)
  % Generate a valid initial model (which may be a poor guess).
  model = generate_initial_model(data, nstates, options);
else
  % Use specified initial model.
  model = initial_model;
end

% Determine number of trajectories provided in 'data'.
ntrajectories = length(data);

if (options.verbosity >= 1)
  disp('Fitting maximum-likelihood HMM with EM algorithm');
end

% Run cycles of maximum-likelihood updating.
for iteration = 1:options.maximumIterations

  if (options.verbosity >= 2)
    disp('******************************************************************************');
    disp(sprintf('EM iteration %d', iteration));
  end

  % Store old model.
  [old_Tij, old_mu, old_sigma] = extract_model_parameters(model);

  % Update emission probability functions associated with each state.
  for i = 1:model.nstates
    state = model.states{i};
    % Define normal density for emission probability.
    model.states{i}.emission_probability = @(o) (2*pi)^(-1/2) * state.sigma^(-1) * exp(-(1/2)*((o - state.mu) / state.sigma)^2); 
    model.states{i}.log_emission_probability = @(o) - (1/2)*log(2*pi) - log(state.sigma) - (1/2)*((o - state.mu) / state.sigma)^2;
  end

  % Calculate the stationary probability distribution from dominant eigenvector of transition matrix.
  model.logTij = log(model.Tij); % elementwise logarithm - can we use a more numerically stable approach later?

  % Compute expected transition and state-occupation quantities.  
  Nij = zeros(model.nstates, model.nstates);
  o_n = []; % o_n(n) is observation n from concatenated trajectories
  w_ni = []; % w_ni(n,i) is probability observation n came from state i
  for trajectory_index = 1:length(data)
    % Extract trajectory of observations.
    o_t = data{trajectory_index}; % o_t(t) is the observation at time t \in 1...T.

    % Compute log probabilities of transitions (log_xi_tij) and observations (log_gamma_ti) with Baum-Welch.
    [log_xi_tij, log_gamma_ti] = baum_welch(o_t, model, options);    

    % Compute expected transition counts.
    Nij = Nij + squeeze(sum(exp(log_xi_tij),1));
    
    % Compute observations and their weights for each state.
    o_n = [o_n; o_t'];
    w_ni = [w_ni; exp(log_gamma_ti)];
  end
  
  % Determine maximum likelihood transition matrix using fractional transition counts.
  model.Tij = transition_matrix_mle(Nij, options);
  model.Pi = stationary_probability(model.Tij);
  model.logPi = log(model.Pi);
  
  % Update state observation probabilities.
  for i = 1:model.nstates
    state = model.states{i};

    % Extract weights for this state.
    w_n = w_ni(:,i);
    
    % Determine weighted sample statistics.
    state.mu = sum(w_n .* o_n) / sum(w_n);
    state.sigma = sqrt(sum(w_n .* (o_n - state.mu).^2) / sum(w_n));
   
    model.states{i} = state;
  end
  
  if (options.verbosity >= 2)
    show_model(model);
  end
    
  % Check convergence criteria by computing relative change.
  [Tij, mu, sigma] = extract_model_parameters(model);  
  if (norm(Tij - old_Tij) / norm(Tij) < RELATIVE_TOLERANCE) && (norm(mu - old_mu) / norm(mu) < RELATIVE_TOLERANCE) && (norm(sigma - old_sigma) / norm(sigma) < RELATIVE_TOLERANCE)
    if (options.verbosity >= 1)
      disp(sprintf('Relative convergence tolerance of %e achieved.', RELATIVE_TOLERANCE));
    end
    break;
  end
end
clear old_model;

% Determine maximum-likelihood state trajectories.
model.state_trajectories = {};
for trajectory_index = 1:length(data)
  o_t = data{trajectory_index}; % o_t(t) is the observation at time t \in 1...T.
  [s_t, log_trajectory_probability] = viterbi(o_t, model);
  model.state_trajectories{trajectory_index} = s_t;
end


if (options.verbosity >= 1)
  disp('******************************************************************************');
  disp('Final converged results:');
  %log_likelihood
  show_model(model);
  disp('******************************************************************************');
end
  
% Return maximum-likelihood model.
mlmodel = model;

return


