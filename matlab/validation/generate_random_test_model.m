function [model, data] = generate_random_test_model(observation_lengths)
% Generate a test model and synthetic data for a random number of states and parameters.
%
%  [model, data] = generate_test_model()
%
% ARGUMENTS
%  observation_lengths - vector of number of observations per trace to generate
%
% RETURNS
%  model - the model from which the data is simulated
%  data  - the synthetic test data generated from the model

% Set up synthetic test system.
disp('******************************************************************************');
disp('Creating synthetic test system...');

% PARAMETERS
min_states = 2; 
max_states = 6;

min_mu = -10;
max_mu = +10;
min_sigma = 0.1;
max_sigma = 3.0;

% DEBUG Seed random number generator deterministically.
%random_seed = 1;
%rand('seed', random_seed);

% Define a model.
model = struct();

% Define number of states.
model.nstates = randi([min_states max_states]);

% Define state emission probabilities.
mu_i = sort(min_mu + (max_mu - min_mu)*rand(1, model.nstates)); % random mu values in ascending order
%sigma_i = min_sigma + (max_sigma - min_sigma)*rand(1, model.nstates); % uniform
sigma_i = exp(log(min_sigma) + (log(max_sigma) - log(min_sigma))*rand(1, model.nstates)); % uniform over log distribution
for i = 1:model.nstates
  model.states{i} = struct('mu', mu_i(i), 'sigma', sigma_i(i));
end

% Store model parameters.
model.tau = 1.0;

% Generate a transition matrix.
Nij = ones(model.nstates,model.nstates);
options = bhmm_default_options();
options.reversibe = 1;
%options.diagonally_dominant = 1;
options.diagonally_dominant = 0;
options.equilibrium = 1;
options.verbosity = 5;
options.prior_transition_pseudocounts = ones(model.nstates, model.nstates);
generate_new_Tij = 1;
while generate_new_Tij
  model.Tij = 0.9 * eye(model.nstates) + 0.1/model.nstates * ones(model.nstates, model.nstates);
  for iteration = 1:10
    model.Tij = transition_matrix_sample_reversible(model.Tij, Nij, options);
  end
  disp('model.Tij = ');
  model.Tij
  % Compute stationary distribution.
  model.Pi = stationary_probability(model.Tij);
  % Check transition matrix is OK.
  generate_new_Tij = 0;
  if any(model.Pi < 0) 
    generate_new_Tij = 1;
  end
  if any(any(model.Tij <= 0))
    generate_new_Tij = 1;
  end
end

% Generate synthetic data.
ntrajectories = length(observation_lengths);
data = cell(ntrajectories,1);
model.state_trajectories = cell(ntrajectories,1);
for trajectory_index = 1:ntrajectories
  % Generate state trajectory.
  T = observation_lengths(trajectory_index);
  s_t = zeros(1,T,'int32'); % s_t(t) is state visited at time t
  s_t(1) = draw(model.Pi);
  for t = 2:T
    s_t(t) = draw(model.Tij(s_t(t-1),:));
  end
  model.state_trajectories{trajectory_index} = s_t; 
  
  % Generate observation trajectory.
  T = observation_lengths(trajectory_index);
  o_t = zeros(1,T,'single'); % o_t(t) is observation at time t
  for t = 1:T
    state = model.states{s_t(t)};
    o_t(t) = randn() * state.sigma + state.mu;
  end
  
  % Pack data into cellarray.
  data{trajectory_index} = o_t;
end

show_model(model);

% DEBUG
%clf;
%plot(o_t, '.');
%pause;

disp('******************************************************************************');

return
