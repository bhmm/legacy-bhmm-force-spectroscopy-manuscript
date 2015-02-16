function [model] = four_state_model()
% Generate a four-state test model.
%
%  [model] = four_state_model()
%
% RETURNS
%  model - a four-state test model

% Define a model.
model = struct();

% Define number of states.
model.nstates = 4;

% Define state emission probabilities.
model.states{1} = struct('mu', 0.1, 'sigma', 0.1);
model.states{2} = struct('mu', 0.3, 'sigma', 0.1);
model.states{3} = struct('mu', 0.6, 'sigma', 0.1);
model.states{4} = struct('mu', 0.9, 'sigma', 0.1);

% Define a transition matrix that satisfies detailed balance and comes from a true rate matrix.
model.Tij = [0.98 0.01  0.005 0.005; 
	     0.01 0.975 0.01  0.05; 
	     0.05 0.01  0.975 0.01; 
	     0.05 0.05  0.01  0.98];
model.Pi = stationary_probability(model.Tij);

% Store model parameters.
model.tau = 1.0;

return
