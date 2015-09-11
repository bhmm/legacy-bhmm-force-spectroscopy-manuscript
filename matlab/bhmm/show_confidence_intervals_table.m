function show_confidence_intervals_table(models, options)
% Generate LaTeX table showing mean and confidence intervals for all properties.
%
% ARGUMENTS
%   models - a set of BHMM-sampled models

% Augment models with mu and sigma vectors.
nmodels = length(models);
nstates = models(1).nstates;
clear augmented_models;
for model_index = 1:nmodels
  model = models(model_index);
  model.mu_i = zeros(model.nstates,1);
  model.sigma_i = zeros(model.nstates,1);
  for i = 1:model.nstates
    model.mu_i(i) = model.states{i}.mu;
    model.sigma_i(i) = model.states{i}.sigma;      
  end
  augmented_models(model_index) = model;
end
models = augmented_models;

% Augment models with rate matrices and lifetimes.
clear augmented_models;
for model_index = 1:nmodels
  model = models(model_index);
  model.Kij = real(logm(model.Tij)) / options.tau;
  model.tau_i = options.tau ./ (1 - diag(model.Tij));
  augmented_models(model_index) = model;
end
models = augmented_models;

% Compute mean and confidence intervals.
ci = 0.95;
properties = {'Pi', 'Tij', 'Kij', 'tau_i', 'mu_i', 'sigma_i'};
nproperties = length(properties);
% Compute model parameters.
property_statistics = struct();
for property_index = 1:nproperties
  property = properties{property_index};
  statistics = struct();
  
  % Collect statistic into vector.
  nmodels = length(models);
  x = getfield(models(1), property);
  property_size = size(x);
  property_dim = length(property_size);
  if property_size(end) == 1
    property_dim = property_dim - 1;
  end
  x_n = zeros([nmodels size(x)]);
  for n = 1:nmodels
    model = models(n);
    x = getfield(model, property);
    if (property_dim==1)
      x_n(n,:) = x;
    else
      x_n(n,:,:) = x;
    end
  end    
  
  % Compute statistics.
  statistics.mean = squeeze(mean(x_n, 1));
  if (property_dim == 1)
    statistics.low = zeros(1,model.nstates);
    statistics.high = zeros(1,model.nstates);
    for i = 1:model.nstates
      [low, high] = empirical_confidence_interval(squeeze(x_n(:,i)), ci);
      statistics.low(i) = low;
      statistics.high(i) = high;
    end
  elseif (property_dim == 2)
    statistics.low = zeros(model.nstates,model.nstates);
    statistics.high = zeros(model.nstates,model.nstates);
    for i = 1:model.nstates
      for j = 1:model.nstates
	[low, high] = empirical_confidence_interval(squeeze(x_n(:,i,j)), ci);
	statistics.low(i,j) = low;
	statistics.high(i,j) = high;
      end      
    end
  end
  
  property_statistics = setfield(property_statistics, property, statistics);
end

% Generate LaTeX table.

% Header
disp('\begin{tabular}{|c|c|c|}');
disp('\hline');
disp('property & value & units \\ \hline');

% \pi_i
property = 'Pi';
for i = 1:nstates
  statistics = getfield(property_statistics, property);
  string = sprintf('$\\pi_{%d}$ & $%.3f_{\\:%.3f}^{\\:%.3f}$ & \\\\', i, statistics.mean(i), statistics.low(i), statistics.high(i));
  disp(string)
end  
disp('\hline');

% T_ij
property = 'Tij';
for i = 1:nstates
  for j = 1:nstates
    statistics = getfield(property_statistics, property);
    string = sprintf('$T_{%d%d}$ & $%.3f_{\\:%.3f}^{\\:%.3f}$ & \\\\', i, j, statistics.mean(i,j), statistics.low(i,j), statistics.high(i,j));
    disp(string);
  end
end  
disp('\hline');

% \mu_i
property = 'mu_i';
for i = 1:nstates
  statistics = getfield(property_statistics, property);
  string = sprintf('$\\mu_{%d}$ & $%.3f_{\\:%.3f}^{\\:%.3f}$ & %s \\\\', i, statistics.mean(i), statistics.low(i), statistics.high(i), options.observable_units);
  disp(string);
end  
disp('\hline');

% \sigma_i
property = 'sigma_i';
for i = 1:nstates
  statistics = getfield(property_statistics, property);
  string = sprintf('$\\sigma_{%d}$ & $%.3f_{\\:%.3f}^{\\:%.3f}$ & %s \\\\', i, statistics.mean(i), statistics.low(i), statistics.high(i), options.observable_units);
  disp(string);
end  
disp('\hline');

% K_ij
property = 'Kij';
for i = 1:nstates
  for j = 1:nstates
    if i ~= j
      statistics = getfield(property_statistics, property);
      string = sprintf('$K_{%d%d}$ & $%.2f_{\\:%.2f}^{\\:%.2f}$ & %s$^{-1}$ \\\\', i, j, statistics.mean(i,j), statistics.low(i,j), statistics.high(i,j), options.time_units);
      disp(string);
    end
  end
end  
disp('\hline');

% \sigma_i
property = 'tau_i';
for i = 1:nstates
  statistics = getfield(property_statistics, property);
  string = sprintf('$\\tau_{%d}$ & $%.1f_{\\:%.1f}^{\\:%.1f}$ & ms \\\\', i, 1000*statistics.mean(i), 1000*statistics.low(i), 1000*statistics.high(i));
  disp(string);
end  
disp('\hline');


disp('\end{tabular}');

return
