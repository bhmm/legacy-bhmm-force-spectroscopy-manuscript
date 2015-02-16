% Generate synthetic data to test the confidence intervals genereated from BHMM sampling.

% Clear workspace.
clear;

% Add BHMM paths.
addpath ../../matlab/bhmm
javaaddpath ../../matlab/bhmm

% Test parameters.
observation_lengths = [10000]; % length of observation trajectories
sigma = 0.2;
true_model = three_state_model(sigma); % test model to use.
nstates = true_model.nstates;
nmodels = 1000; % number of models to sample with Bayesian HMM
nburnin = 50; % number of samples to discard to equilibration

% Set parameters.
options = bhmm_default_options(); % get default options structure
options.maximumIterations = 10; % maximum number of allowed EM iterations
options.convergenceTolerance = 1.0e-4; % relative convergence tolerance in EM iterations
options.tau = 1; % time interval between observations

% Options for transition matrix estimation.
options.reversible = true; % infer reversible transition matrices
options.diagonally_dominant = false; % don't enforce diagonally-dominant transition matrices
options.verbosity = 3; % set verbosity level
options.equilibrium = true; % trajectory data is not initially drawn from equilibrium

% Options for plots.
options = rmfield(options, 'observable_units'); % no units for observable axis
options = rmfield(options, 'time_units'); % no units for time axis
options.observable_name = 'obs';

% Turn on state sampling.
options.enforce_state_ordering = 0; % don't enforce state ordering for now
options.sample_states = 1; % sample over number of states

% Show the true model.
show_model(true_model);

%% Generate initial data and model to illustrate.

% Generate synthetic data from the true model.
data = generate_synthetic_data(true_model, observation_lengths);

% Determine MLHMM.
disp('Determining MLHMM...');
mlmodel = mlhmm(data, true_model.nstates, options);

% DEBUG
models = [mlmodel];
for i = 1:100
  model = models(end);
  models = bhmm(data, model, 1, options);
  plot_state_assignments(data, models(end), options);
  pause
end

% Initialize Bayesian HMM.
disp('Sampling from Bayesian posterior...');
models = bhmm(data, mlmodel, nburnin+nmodels, options);

% Discard initial samples to burn-in.
models = models(nburnin+1:end)

% Plot MLHMM state assignments.
%plot_fractional_state_assignments(data, models(nburnin+1:end), options);
plot_state_assignments(data, mlmodel, options);
%filename = 'stateassignments.eps';
%print('-deps', filename);
%command = sprintf('epstopdf %s', filename);
%system(command);
% Export figure.
%exportfig(gcf, 'model-stateassignments.eps', 'color', 'cmyk', 'width', 7.5, 'height', 1.5);
%system('epstopdf model-stateassignments.eps');
filename = 'model-variablestates-mlhmm-stateassignments.eps';
% Generate separate text and figure
resolution = 450; % dpi for figures
scale = 1.4; % figure scaling
addpath exportfig
addpath epscombine
%exportfig(gcf, 'zbuffer', 'color', 'cmyk', 'width', scale*7.5, 'height', scale*1.5, 'separatetext', 1, 'renderer', 'zbuffer', 'resolution', resolution);
exportfig(gcf, 'zbuffer', 'color', 'cmyk', 'width', scale*7.5, 'height', scale*1.5, 'separatetext', 1, 'renderer', 'opengl', 'resolution', resolution);
exportfig(gcf, 'painter', 'color', 'cmyk', 'width', scale*7.5, 'height', scale*1.5, 'separatetext', 1, 'renderer', 'painters');
epscombine('painter_t.eps','zbuffer.eps',filename);
system(sprintf('epstopdf %s', filename));

% Plot BHMM state assignments.
%plot_fractional_state_assignments(data, models(nburnin+1:end), options);
plot_state_assignments(data, models(end), options);
%filename = 'stateassignments.eps';
%print('-deps', filename);
%command = sprintf('epstopdf %s', filename);
%system(command);
% Export figure.
%exportfig(gcf, 'model-stateassignments.eps', 'color', 'cmyk', 'width', 7.5, 'height', 1.5);
%system('epstopdf model-stateassignments.eps');
filename = 'model-variablestates-stateassignments.eps';
% Generate separate text and figure
resolution = 450; % dpi for figures
scale = 1.4; % figure scaling
addpath exportfig
addpath epscombine
%exportfig(gcf, 'zbuffer', 'color', 'cmyk', 'width', scale*7.5, 'height', scale*1.5, 'separatetext', 1, 'renderer', 'zbuffer', 'resolution', resolution);
exportfig(gcf, 'zbuffer', 'color', 'cmyk', 'width', scale*7.5, 'height', scale*1.5, 'separatetext', 1, 'renderer', 'opengl', 'resolution', resolution);
exportfig(gcf, 'painter', 'color', 'cmyk', 'width', scale*7.5, 'height', scale*1.5, 'separatetext', 1, 'renderer', 'painters');
epscombine('painter_t.eps','zbuffer.eps',filename);
system(sprintf('epstopdf %s', filename));

% Plot states.
options.sameaxis = true;
estimates = analyze_states(models, options);
%filename = 'states.eps';
%print('-depsc', filename);
%command = sprintf('epstopdf %s', filename);
%system(command);
exportfig(gcf, 'model-variablestates-state-parameters.eps', 'color', 'cmyk', 'width', 3.5, 'height', 3.5);
system('epstopdf model-variablestates-state-parameters.eps');

