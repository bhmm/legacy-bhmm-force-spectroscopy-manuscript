% Generate synthetic data to test the confidence intervals genereated from BHMM sampling.

% Clear workspace.
clear;

% Add BHMM paths.
addpath ../../matlab/bhmm
javaaddpath ../../matlab/bhmm

% Test parameters.
ntrials = 50; % number of independent realizations of the experiment to generate
cis = 0.05:0.05:0.95; % confidence intervals to evaluate
ncis = length(cis);
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

% Show the true model.
show_model(true_model);

%% Generate initial data and model to illustrate.

% Generate synthetic data from the true model.
data = generate_synthetic_data(true_model, observation_lengths);

% Determine MLHMM.
disp('Determining MLHMM...');
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
filename = 'model-stateassignments.eps';
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
exportfig(gcf, 'model-state-parameters.eps', 'color', 'cmyk', 'width', 3.5, 'height', 3.5);
system('epstopdf model-state-parameters.eps');

% DEBUG: STOP.
return;

%% Run many uncorrelated experiments to test confidence intervals.

% Allocate data.
ci_Pi = zeros(ncis,nstates);
ci_Tij = zeros(ncis,nstates,nstates);
ci_mu_i = zeros(ncis,nstates);
ci_sigma_i = zeros(ncis,nstates);

for trial = 1:ntrials
  % Generate synthetic data from the true model.
  data = generate_synthetic_data(true_model, observation_lengths);

  % Determine MLHMM.
  disp('Determining MLHMM...');
  mlmodel = mlhmm(data, true_model.nstates, options);
  
  % Initialize Bayesian HMM.
  disp('Sampling from Bayesian posterior...');
  models = bhmm(data, mlmodel, nburnin+nmodels, options);

  % Discard initial samples to burn-in.
  models = models(nburnin+1:end)

  % Extract model parameters.
  Pi_n = zeros(nstates,nmodels);
  Tij_n = zeros(nstates,nstates,nmodels);
  mu_i_n = zeros(nstates,nmodels);
  sigma_i_n = zeros(nstates,nmodels);
  for n = 1:nmodels
    Pi_n(:,n) = models(n).Pi;
    Tij_n(:,:,n) = models(n).Tij;
    for i = 1:nstates
      mu_i_n(i,n) = models(n).states{i}.mu;
      sigma_i_n(i,n) = models(n).states{i}.sigma;
    end
  end

  % Accumulate statistics for confidence intervals.
  for ci_index = 1:ncis
    % Get confidence interval to check.
    ci = cis(ci_index);
    ci
    
    for i = 1:nstates
      % Pi
      [low, high] = empirical_confidence_interval(Pi_n(i,:), ci);
      if (low <= true_model.Pi(i)) && (true_model.Pi(i) <= high)
	ci_Pi(ci_index,i) = ci_Pi(ci_index,i) + 1;
      end

      % Tij
      for j = 1:nstates
	[low, high] = empirical_confidence_interval(Tij_n(i,j,:), ci);
	if (low <= true_model.Tij(i,j)) && (true_model.Tij(i,j) <= high)
	  ci_Tij(ci_index,i,j) = ci_Tij(ci_index,i,j) + 1;
	end
      end

      % mu
      [low, high] = empirical_confidence_interval(mu_i_n(i,:), ci);
      if (low <= true_model.states{i}.mu) && (true_model.states{i}.mu <= high)
	ci_mu_i(ci_index,i) = ci_mu_i(ci_index,i) + 1;
      end

      % sigma
      [low, high] = empirical_confidence_interval(sigma_i_n(i,:), ci);
      if (low <= true_model.states{i}.sigma) && (true_model.states{i}.sigma <= high)
	ci_sigma_i(ci_index,i) = ci_sigma_i(ci_index,i) + 1;
      end      

    end    
    
  end
end

ci_Pi = ci_Pi / ntrials;
ci_Tij = ci_Tij / ntrials;
ci_mu_i = ci_mu_i / ntrials;
ci_sigma_i = ci_sigma_i / ntrials;

% Save all data.
save three_state_model.mat


