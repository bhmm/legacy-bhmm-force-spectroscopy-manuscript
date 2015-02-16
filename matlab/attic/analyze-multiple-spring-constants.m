% Analyze data from multiple spring constants.

% Clear all
clear all;
clear java;

% Add paths to BHMM package.
addpath /Users/jchodera/Documents/latex/doc/papers/2009/bhmm/matlab/bhmm
javaaddpath /Users/jchodera/Documents/latex/doc/papers/2009/bhmm/matlab/bhmm

% PARAMETERS
matfilename = '10khz/fibers_l_10kHz.mat'; % file to analyze
nsplit = 10; % number of traces to split into.

% Set BHMM options.
nstates = 2;
nsamples = 100; % number of samples to generate

options = bhmm_default_options();
options.maximumIterations = 3; % only need 30 iterations of EM for mlhmm
options.reversible = true; % infer reversible transition matrices
options.diagonally_dominant = false; % don't enforce diagonally-dominant transition matrices
options.verbosity = 2; % set verbosity level
options.equilibrium = false; % trajectory data is not initially drawn from equilibrium
options.tau = 1.0 / (10.0e3); % seconds between observations

% Get list of fibers in each file.
disp(sprintf('querying fibers in "%s"...', matfilename));
fiber_names = who('-file', matfilename);
nfibers = length(fiber_names);

for fiber_index = 1:nfibers
  fiber_name = fiber_names{fiber_index};
  disp(sprintf('processing fiber "%s"...', fiber_name));

  % Load fiber data.
  contents = load(matfilename, fiber_name);
  fiber = getfield(contents, fiber_name);
  clear contents;

  % Get number of traces.
  ntraces = length(fiber);
  for trace_index = 1:ntraces
    disp(sprintf('trace %d / %d', trace_index, ntraces));

    % Extract trace.
    o_t = fiber{trace_index};
    T = length(o_t);

    % Break trace into smaller parts that are less sensitive to drift.
    data = {};
    for i = 1:nsplit
      indices = ((i-1)*floor(T/nsplit)+1):(i*floor(T/nsplit));
      data{i} = o_t(indices)';
    end

    % Generate an initial model from data.
    initial_model = generate_initial_model(data, nstates, options);
    
    % Fit a maximum-likelihood HMM.
    %disp('Determining maximum-likelihood model...');
    %mlmodel = mlhmm(data, nstates, options);
    
    % Initialize Bayesian HMM.
    disp('Sampling from Bayesian posterior...');    
    models = bhmm(data, initial_model, nsamples, options);

    % Plot
    nburnin = 50; % number of samples to discard to equilibration
    plot_state_assignments(data, models(end), options);
    options.sameaxis = false;
    analyze_states(models(nburnin+1:end, options));

    pause;
  end
end




