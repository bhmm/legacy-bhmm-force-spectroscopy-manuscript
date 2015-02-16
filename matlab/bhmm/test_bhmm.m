% Test Bayesian HMM with synthetic data.

% Clear workspace.
clear;

% PARAMETERS
T = 10000; % observation length

% Generate test model and data.
[true_model, data] = generate_test_model(T);
true_model.Pi

% Set parameters.
options = bhmm_default_options(); % get default options structure
options.maximumIterations = 10; % maximum number of allowed iterations
options.convergenceTolerance = 1.0e-4; % relative convergence tolerance in likelihood
options.tau = 1;

% Options for transition matrix estimation.
options.reversible = true; % infer reversible transition matrices
options.diagonally_dominant = false; % don't enforce diagonally-dominant transition matrices
options.verbosity = 2; % set verbosity level
options.equilibrium = false; % trajectory data is not initially drawn from equilibrium

% Determine maximum-likelihood HMM.
disp('Determining maximum-likelihood model...');
mlmodel = mlhmm(data, true_model.nstates, options);

% Sample from Bayesian HMM starting from maximum-likelihood model.
disp('Sampling from Bayesian posterior...');
nsamples = 100;
models = bhmm(data, mlmodel, nsamples, options);

% Plot
plot_state_assignments(data, models(end), options);


