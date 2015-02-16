function [model, data] = generate_test_model(T)
% Generate a test model and synthetic data.
%
%  [model, data] = generate_test_model()
%
% RETURNS
%  model - the model from which the data is simulated
%  data  - the synthetic test data generated from the model
%
% NOTE
%  The random number seed is set deterministically.

% Set up synthetic test system.
disp('******************************************************************************');
disp('Creating synthetic test system...');

% DEBUG Seed random number generator deterministically.
rand('seed', 1);

% Define a model.
model = struct();

% Define number of states.
model.nstates = 3;

% Define state emission probabilities.
model.states{1} = struct('mu', 3.0, 'sigma', 1.0); % unfolded
model.states{2} = struct('mu', 4.7, 'sigma', 0.3); % intermediate
model.states{3} = struct('mu', 5.6, 'sigma', 0.2); % folded

% Define a transition matrix that satisfies detailed balance and comes from a true rate matrix.
model.Tij = [    0.980 0.019 0.001;
		 0.05  0.90  0.05;
		 0.001 0.009 0.990];
        
model.Pi = stationary_probability(model.Tij);

% Ensure that Tij satisfies detailed balance.
for i = 1:model.nstates
  for j = 1:(i-1)
    model.Tij(i,j) = model.Pi(j) * model.Tij(j,i) / model.Pi(i);
  end
end
model.Pi = stationary_probability(model.Tij);

% Store model parameters.
model.tau = 0.001;

% Generate synthetic data.
%T = 50000; % trajectory length

% Generate state trajectory.
s_t = zeros(1,T,'int32'); % s_t(t) is state visited at time t
s_t(1) = draw(model.Pi);
for t = 2:T
  s_t(t) = draw(model.Tij(s_t(t-1),:));
end
model.state_trajectories{1} = s_t; 

% Generate observation trajectory.
o_t = zeros(1,T,'single'); % o_t(t) is observation at time t
for t = 1:T
  state = model.states{s_t(t)};
  o_t(t) = randn() * state.sigma + state.mu;
end

% Pack data into cellarray.
data{1} = o_t;

show_model(model);

disp('******************************************************************************');

return
