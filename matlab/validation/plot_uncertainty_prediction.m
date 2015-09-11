% Plot uncertainty prediction.

nmodels = 50

%clear;
%load true_uncertainties.mat;

% Define property to examine.
%property = @(model) model.Tij(1,2) / model.Tij(1,3);
property = @(model) model.Tij(1,2);
%property = @(model) model.Pi(2);
%property = @(model) 1.0 / (1 - model.Tij(2,2)); % lifetime of intermediate state

% Predict uncertainty with BHMM modeling.
actual_uncertainty = zeros(max_traces,1);
actual_uncertainty_error = zeros(max_traces,1);
for ntraces = 1:max_traces
  % Sample a certain number of synethetic experiments for prediction.
  bhmm_estimates = zeros(nsynthetic,1); 
  bhmm_stddevs = zeros(nsynthetic,1); 
  for synthetic_index = 1:nsynthetic
    [ntraces, synthetic_index];

    % Load models.
    bhmm_models = replicates{ntraces,synthetic_index}.bhmm_models;

    % Store statistics.
    property_n = zeros(nmodels,1); % storage for property
    for n = 1:nmodels
      property_n(n) = property(bhmm_models(n));
    end
    bhmm_estimates(synthetic_index) = mean(property_n);
    bhmm_stddevs(synthetic_index) = std(property_n);        
  end

  % Store information.
  actual_uncertainty(ntraces) = mean(bhmm_stddevs);
  actual_uncertainty_error(ntraces) = std(bhmm_stddevs) / sqrt(nsynthetic);
end

% Plot.
clf;
plot(1:max_traces, actual_uncertainty(1:max_traces), 'k.');
hold on;
errorbar(1:max_traces, predicted_uncertainty(1:max_traces), 2*predicted_uncertainty_error(1:max_traces), 2*predicted_uncertainty_error(1:max_traces), 'r.');
%hold on;
%plot(1:max_traces, predicted_uncertainty(1) ./ sqrt(1:max_traces), 'k-');
legend('actual', 'BHMM predicted');
xlabel('number of traces collected');
ylabel('uncertainty');
axis([0 (max_traces+1) 0 max(predicted_uncertainty+2*predicted_uncertainty_error)*1.01]);


