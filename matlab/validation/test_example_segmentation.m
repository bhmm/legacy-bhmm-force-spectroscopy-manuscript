% Generate synthetic data to test the estimates produced by simple segmentation of the observable.

% Clear workspace.
clear;

% Add BHMM paths.
addpath ../../matlab/bhmm
javaaddpath ../../matlab/bhmm

% Test parameters.
L_n = [1e3 1e4 1e5 1e6];
NL = length(L_n);

generate_plots = true;

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
options.convergenceTolerance = 0.01; % relative convergence tolerance in EM iterations
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

% Generate segmentation models.

% Storage for models.
model_store = cell(NL, 1);

for length_index = 1:NL
  % Get observable trace length to use.
  L = L_n(length_index);

  % Extract subset of data.
  data = cell(1,1);
  o_t = all_data{1};
  data{1} = o_t(1:L);

  % Determine segmentation estimate
  segmentation_model = generate_segmentation_model(data, true_model.nstates, options);

  if generate_plots
    % Plot state assignments.
    plot_state_assignments(data, segmentation_model, options);
    filename = sprintf('model-segmentation-stateassignments-%d.eps', L);
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
  end

  % Store models.
  model_store{length_index} = segmentation_model;
end

% DEBUG: STOP.
%return;

% Save data.
save test-example-data-segmentation.mat;

test_table_segmentation