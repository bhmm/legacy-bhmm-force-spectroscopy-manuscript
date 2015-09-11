% Generate synthetic data to test the confidence intervals genereated from BHMM sampling.

% Clear workspace.
clear;

% Add BHMM paths.
addpath ../../matlab/bhmm
javaaddpath ../../matlab/bhmm

% Test parameters.
L_n = [1e3 1e4 1e5];
NL = length(L_n);

% Generate model and synthetic data for maximum length.
[true_model, all_data] = generate_test_model(max(L_n));
nstates = true_model.nstates;

cis = 0.05:0.05:0.95; % confidence intervals to evaluate
ncis = length(cis);
nmodels = 1000; % number of models to sample with Bayesian HMM
nburnin = 50; % number of samples to discard to equilibration

% Set parameters.
options = bhmm_default_options(); % get default options structure
options.maximumIterations = 10; % maximum number of allowed EM iterations
options.convergenceTolerance = 1.0e-4; % relative convergence tolerance in EM iterations
options.tau = 0.001; % time interval between observations (s)
options.time_units = 's';

% Options for transition matrix estimation.
options.reversible = true; % infer reversible transition matrices
options.diagonally_dominant = false; % don't enforce diagonally-dominant transition matrices
options.verbosity = 3; % set verbosity level
options.equilibrium = true; % trajectory data is not initially drawn from equilibrium

% Options for plots.
%options = rmfield(options, 'observable_units'); % no units for observable axis
%options = rmfield(options, 'time_units'); % no units for time axis
%options.observable_name = 'obs';

% Show the true model.
disp('true model:');
show_model(true_model);

% Determine MLHMM.
disp('Determining MLHMM...');
observation_length = 1e4;
data = cell(1,1);
o_t = all_data{1};
data{1} = o_t(1:observation_length);
mlmodel = mlhmm(data, true_model.nstates, options);

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
filename = 'model-mlhmm-stateassignments.eps';
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

% Storage for models.
model_store = cell(NL, 1);

for length_index = 1:NL
  % Get observable trace length to use.
  L = L_n(length_index);

  % Extract subset of data.
  data = cell(1,1);
  o_t = all_data{1};
  data{1} = o_t(1:L);

  % Initialize Bayesian HMM.
  disp('Sampling from Bayesian posterior...');
  models = bhmm(data, mlmodel, nburnin+nmodels, options);
  
  % Discard initial samples to burn-in.
  models = models(nburnin+1:end)
  
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
  filename = sprintf('model-stateassignments-%d.eps', L);
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

  % Plot state characterization.
  %options.sameaxis = true;
  %estimates = analyze_states(models, options);
  %exportfig(gcf, 'model-state-parameters.eps', 'color', 'cmyk', 'width', 3.5, 'height', 3.5);
  %system('epstopdf model-state-parameters.eps');

  % save workspace
  %save test-example-data.mat

  % Store models.
  model_store{length_index} = models;
end

% DEBUG: STOP.
%return;

% Save data.
save test-example-data.mat;

