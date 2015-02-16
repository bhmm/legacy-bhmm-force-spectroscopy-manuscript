% Test prediction of uncertainty.

generate_bhmm = true;
if generate_bhmm
  clear;

  % Add BHMM paths.
  addpath ../../matlab/bhmm
  javaaddpath ../../matlab/bhmm

  load true_uncertainties.mat;

  % Generate true model and initial observation data for first trace.
  [true_model, data] = generate_test_model(observation_length);
  nstates = true_model.nstates;

  % Build MLHMM.
  mlhmm_model = mlhmm(data, true_model.nstates, options);

  % Generate Bayesian HMM and store predicted models.
  nburnin = 100;
  nmodels = 10000;
  bhmm_models = bhmm(data, mlhmm_model, nburnin+nmodels, options);    
  
  % Discard initial models to burn-in.
  bhmm_models = bhmm_models(nburnin+1:end);

  % Remove hidden state sequences.
  clear pruned_models;
  for i = 1:nmodels
    model = bhmm_models(i);
    model = rmfield(model, 'state_trajectories');
    pruned_models(i) = model;
  end
  bhmm_models = pruned_models;
  clear pruned_models;
  
  save uncertainty_prediction_models data bhmm_models;
else
  load uncertainty_prediction_models;
end

% Compute original log probabilities of sampled models.
log_sampled_probabilities = zeros(nmodels,1);
for model_index = 1:nmodels
  log_sampled_probabilities(model_index) = compute_log_posterior(bhmm_models(model_index), data, options, false);
end    

% Define property to examine.
%property = @(model) model.Tij(1,2) / model.Tij(1,3);
property = @(model) model.Tij(1,2);
%property = @(model) model.Pi(2);
%property = @(model) 1.0 / (1 - model.Tij(2,2)); % lifetime of intermediate state

% Compute property for all models.
property_n = zeros(nmodels,1); % storage for property
for n = 1:nmodels
  property_n(n) = property(bhmm_models(n));
end

% Predict uncertainty with BHMM modeling.
predicted_uncertainty = zeros(max_traces,1);
predicted_uncertainty_error = zeros(max_traces,1);
for ntraces = 1:max_traces
  % Sample a certain number of synethetic experiments for prediction.
  mlhmm_estimates = zeros(nsynthetic,1); 
  bhmm_predicted_estimates = zeros(nsynthetic,1); 
  bhmm_predicted_stddevs = zeros(nsynthetic,1); 
  for synthetic_index = 1:nsynthetic
    [ntraces, synthetic_index]

    % Select random model to simulate from
    model_index = randi(nmodels,1);
    simulation_model = bhmm_models(model_index);

    % Generate synthetic data for unobserved traces from model, keeping original data.
    synthetic_data = generate_synthetic_data(simulation_model, observation_length * ones(ntraces,1));

    % Reweight all models using synthetic data to compute new expectation, as if we had sampled from new posterior.
    log_model_weights = zeros(nmodels,1);
    for model_index = 1:nmodels
      % Compute model weight due to likelihood of new data.
      log_posterior = compute_log_posterior(bhmm_models(model_index), synthetic_data, options, false);
      log_model_weights(model_index) = log_posterior - log_sampled_probabilities(model_index);
    end    
    max_log_weight = max(log_model_weights);
    model_weights = exp(log_model_weights - max_log_weight);
    model_weights = model_weights / sum(model_weights);
    disp(sprintf('%.4f models contribute', sum(model_weights / max(model_weights))));

    % Compute property statistics.    
    EP = sum(model_weights .* property_n);
    varP = sum(model_weights .* (property_n - EP).^2);
    bhmm_predicted_estimates(synthetic_index) = EP;
    bhmm_predicted_stddevs(synthetic_index) = sqrt(varP)
  end

  % Store information.
  predicted_uncertainty(ntraces) = mean(bhmm_predicted_stddevs);
  predicted_uncertainty_error(ntraces) = std(bhmm_predicted_stddevs);
end

save predicted_uncertainty.mat predicted_uncertainty predicted_uncertainty_error;

% Plot.
clf;
errorbar(1:max_traces, predicted_uncertainty(1:max_traces), 2*predicted_uncertainty_error(1:max_traces), 2*predicted_uncertainty_error(1:max_traces), 'r.');
legend('actual', 'BHMM predicted');
xlabel('number of traces collected');
ylabel('uncertainty');
axis([0 (max_traces+1) 0 max(predicted_uncertainty+2*predicted_uncertainty_error)*1.01]);


