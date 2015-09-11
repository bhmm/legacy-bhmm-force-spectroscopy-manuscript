% Test prediction of how uncertainty will go decrease with additional data collection.

% Clear workspace.
clear;

% Add BHMM paths.
addpath ../../matlab/bhmm
javaaddpath ../../matlab/bhmm

% Test parameters.
observation_length = 1000;
max_traces = 7;

% Generate true model and initial observation data for first trace.
[true_model, data] = generate_test_model(observation_length);
nstates = true_model.nstates;

cis = 0.05:0.05:0.95; % confidence intervals to evaluate
ncis = length(cis);

nsynthetic = 50; % number of synthetic experiments to attempt for prediction
nburnin = 50; % number of samples to discard to equilibration
nmodels = 50; % number of models to sample with Bayesian HMM (after burn-in)

% Set parameters.
options = bhmm_default_options(); % get default options structure
options.maximumIterations = 1000; % maximum number of allowed EM iterations
options.convergenceTolerance = 1.0e-4; % relative convergence tolerance in EM iterations
options.tau = 0.001; % time interval between observations (s)
options.time_units = 's';

% Options for transition matrix estimation.
options.reversible = true; % infer reversible transition matrices
options.diagonally_dominant = false; % don't enforce diagonally-dominant transition matrices
options.verbosity = 1; % set verbosity level
options.equilibrium = true; % trajectory data is not initially drawn from equilibrium

% Show the true model.
disp('true model:');
show_model(true_model);

% Compute true uncertainty from realizations from true model.
replicates = cell(max_traces,nsynthetic); % replicates(trace_index,synthetic_index) is a struct containing observed data, mlhmm, and bhmm
for ntraces = 1:max_traces
  % Sample a certain number of true experiments.
  for synthetic_index = 1:nsynthetic
    [ntraces, synthetic_index]

    % Generate synthetic data for true model.
    data = generate_synthetic_data(true_model, observation_length * ones(ntraces,1));

    % Build MLHMM.
    mlhmm_model = mlhmm(data, true_model.nstates, options);

    % Generate Bayesian HMM and store predicted models.
    bhmm_models = bhmm(data, mlhmm_model, nburnin+nmodels, options);
    
    % Discard models to burnin.
    bhmm_models = bhmm_models(nburnin+1:end);

    % Remove hidden state sequences.
    clear pruned_models;
    for i = 1:nmodels
      model = bhmm_models(i);
      model = rmfield(model, 'state_trajectories');
      pruned_models(i) = model;
    end

    % Store BHMM models.
    replicate = struct();
    replicate.data = data;
    %replicate.mlhmm_model = mlhmm_model;
    replicate.bhmm_models = pruned_models;
    replicates{ntraces,synthetic_index} = replicate;
  end
end

% Save data.
save true_uncertainties.mat;

