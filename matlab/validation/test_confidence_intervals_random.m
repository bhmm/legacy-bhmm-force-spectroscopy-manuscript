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
options.equilibrium = true; % trajectory data is initially drawn from equilibrium

% Options for plots.
options = rmfield(options, 'observable_units'); % no units for observable axis
options = rmfield(options, 'time_units'); % no units for time axis
options.observable_name = 'obs';

%% Run many uncorrelated experiments to test confidence intervals.

% Allocate data.
ci_Pi = zeros(ncis,1);
ci_Pi_denom = zeros(ncis,1);
ci_Tij = zeros(ncis,1);
ci_Tij_denom = zeros(ncis,1);
ci_mu_i = zeros(ncis,1);
ci_mu_i_denom = zeros(ncis,1);
ci_sigma_i = zeros(ncis,1);
ci_sigma_i_denom = zeros(ncis,1);

for trial = 1:ntrials
  % Generate synthetic model and observation data.
  [true_model, data] = generate_random_test_model(observation_lengths);

  % Get number of states.
  nstates = true_model.nstates;

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
	ci_Pi(ci_index) = ci_Pi(ci_index) + 1;
      end
      ci_Pi_denom(ci_index) = ci_Pi_denom(ci_index) + 1;

      % Tij
      for j = 1:nstates
	[low, high] = empirical_confidence_interval(Tij_n(i,j,:), ci);
	if (low <= true_model.Tij(i,j)) && (true_model.Tij(i,j) <= high)
	  ci_Tij(ci_index) = ci_Tij(ci_index) + 1;
	end
	ci_Tij_denom(ci_index) = ci_Tij_denom(ci_index) + 1;
      end

      % mu
      [low, high] = empirical_confidence_interval(mu_i_n(i,:), ci);
      if (low <= true_model.states{i}.mu) && (true_model.states{i}.mu <= high)
	ci_mu_i(ci_index) = ci_mu_i(ci_index) + 1;
      end
      ci_mu_i_denom(ci_index) = ci_mu_i_denom(ci_index) + 1;

      % sigma
      [low, high] = empirical_confidence_interval(sigma_i_n(i,:), ci);
      if (low <= true_model.states{i}.sigma) && (true_model.states{i}.sigma <= high)
	ci_sigma_i(ci_index) = ci_sigma_i(ci_index) + 1;
      end      
      ci_sigma_i_denom(ci_index) = ci_sigma_i_denom(ci_index) + 1;
    end    
    
  end
end

ci_Pi = ci_Pi ./ ci_Pi_denom;
ci_Tij = ci_Tij ./ ci_Tij_denom;
ci_mu_i = ci_mu_i ./ ci_mu_i_denom;
ci_sigma_i = ci_sigma_i ./ ci_sigma_i_denom;

% Save all data.
save random_model_validation.mat


