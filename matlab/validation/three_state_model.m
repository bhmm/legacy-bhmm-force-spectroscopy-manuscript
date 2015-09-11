function [model] = three_state_model(sigma)
% Generate a three-state test model.
%
%  model = three_state_model(sigma)
%
% ARGUMENTS
%  sigma - the std dev for the state widths (states are spaced 1 unit apart)
%
% RETURNS
%  model - a three-state test model

% Define a model.
model = struct();

% Define number of states.
model.nstates = 3;

% Define state emission probabilities.
model.states{1} = struct('mu', -1, 'sigma', sigma);
model.states{2} = struct('mu',  0, 'sigma', sigma);
model.states{3} = struct('mu', +1, 'sigma', sigma);

% Define row-stochastic rate matrix that satisfies detailed balance, and compute transition matrix from this.
Kij = [-0.10  0.10  0.00;
        0.10 -0.15  0.05;
        0.00  0.05 -0.05];
model.Tij = expm(Kij);

% Comptue stationary probability.
model.Pi = stationary_probability(model.Tij);

% Store model parameters.
model.tau = 1.0;

return
