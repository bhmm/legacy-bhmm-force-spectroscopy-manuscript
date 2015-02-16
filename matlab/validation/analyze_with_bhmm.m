% Analyze traces with Bayesian HMM to compute average forces and lifetimes for different trap extensions.

clear;

addpath ../../matlab/bhmm
javaaddpath ../../matlab/bhmm

fibername = 'fiber3';
fiber_data = sprintf('calibrated/%s.mat', fibername); % file to load, containing 'trace' cellarray with multiple traces to be analyzed.
subsample = 50; % subsampling fraction -- subsample to 1 kHz
sampletime = (1/50000) * subsample; % time between samples, in seconds
nstates = 3; % number of states -- fixed beforehand
nsamples = 1000; % number of models to sample with Bayesian HMM
nburnin = 50; % number of samples to discard to equilibration

% Load data.
load(fiber_data);
ntraces = length(trace);

% Allocate storage for estimates
estimates = cell(ntraces,1); % estimates{trace_index} is a struct containing estimates of state properties for each trace

for trace_index = 1:ntraces  
  % Extract data.
  o_t = trace{trace_index}';

  % Subsample data.
  o_t = o_t(1:subsample:end);

  % Pack just this trace into a cell array.
  data = {};
  data{1} = o_t;

  % Determine maximum-likelihood HMM.
  %disp('Determining maximum-likelihood model...');
  %mlmodel = mlhmm(data, model.nstates);
  
  options = bhmm_default_options();
  options.maximumIterations = 10; % maximum number of allowed iterations
  options.convergenceTolerance = 1.0e-4; % relative convergence tolerance in likelihood

  % Options for transition matrix estimation.
  options.reversible = true; % infer reversible transition matrices
  options.diagonally_dominant = false; % don't enforce diagonally-dominant transition matrices
  options.verbosity = 2; % set verbosity level
  options.equilibrium = false; % trajectory data is not initially drawn from equilibrium
  options.maximumIterations = 2; % maximum number of EM iterations for maximum-likelihood HMM (WARNING: THESE ITERATIONS ARE SLOW)

  options.tau = sampletime; % seconds between observations

  % Generate an initial (poor) model from data.
  %initial_model = generate_initial_model(data, nstates, options);
  
  % Determine maximum-likelihood HMM.
  disp('Determining maximum-likelihood model...');
  mlmodel = mlhmm(data, nstates, options);
  
  % Initialize Bayesian HMM.
  disp('Sampling from Bayesian posterior...');
  models = bhmm(data, mlmodel, nsamples, options);

  %% Analayze Bayesian models.
  
  % Plot state assignments.
  %plot_fractional_state_assignments(data, models(nburning+1:end), options);
  plot_state_assignments(data, models(end), options);
  filename = sprintf('plots/%s-trace%d-stateassignments.png', fibername, trace_index);
  print('-dpng', filename);
  %command = sprintf('epstopdf %s', filename);
  %system(command);
  
  % Plot states.
  options.sameaxis = false;
  estimates{trace_index} = analyze_states(models(nburnin+1:end), options);
  filename = sprintf('plots/%s-trace%d-states.eps', fibername, trace_index);
  print('-depsc', filename);
  command = sprintf('epstopdf %s', filename);
  system(command);
end

% Save all estimates of state mean forces and lifetimes.
filename = sprintf('analysis/%s-estimates.mat', fibername);
save(filename, 'estimates');

