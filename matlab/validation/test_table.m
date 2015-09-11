% Generate LaTeX table showing mean and confidence intervals for all properties, for different trajectory lengths.

clear;

% Load three-state example data.
load test-example-data.mat;

% Extract mu_i and sigma_i for true model.
model = true_model;
model.mu_i = zeros(model.nstates,1);
model.sigma_i = zeros(model.nstates,1);
for i = 1:model.nstates
  model.mu_i(i) = model.states{i}.mu;
  model.sigma_i(i) = model.states{i}.sigma;      
end
true_model = model;

% Collect mu and sigma into vectors for models.
nlengths = length(L_n);
for length_index = 1:NL
  % Extract models.
  models = model_store{length_index};  
  % Process models.
  nmodels = length(models);
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
  model_store{length_index} = augmented_models;  
end

% Compute mean and confidence intervals.
ci = 0.95;
properties = {'Pi', 'Tij', 'mu_i', 'sigma_i'};
nproperties = length(properties);
nlengths = length(L_n);
property_statistics = cell(nlengths,1);
for length_index = 1:NL
  % Extract models.
  models = model_store{length_index};

  % Compute model parameters.
  property_statistics{length_index} = struct();
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
    
    property_statistics{length_index} = setfield(property_statistics{length_index}, property, statistics);
  end
end

% Generate table that goes the other way.
disp(' ');
disp(' ');

string = '\begin{tabular}{|c|c|';
for length_index = 1:NL
  string = sprintf('%sc', string);
end
string = sprintf('%s|}',string);
disp(string)
disp('\hline');

% header line
disp(sprintf('property & true value & \\multicolumn{%d}{c|}{observation length} \\\\', nlengths));
string = ' &';
for length_index = 1:NL
  string = sprintf('%s & %d', string, L_n(length_index));
end
string = sprintf('%s \\\\ \\hline', string);
disp(string)

% \pi_i
property = 'Pi';
for i = 1:nstates
  string = sprintf('$\\pi_{%d}$ & $%.3f$', i, true_model.Pi(i));
  for length_index = 1:NL
    statistics = getfield(property_statistics{length_index}, property);
    string = sprintf('%s & $%.3f_{\\:%.3f}^{\\:%.3f}$', string, statistics.mean(i), statistics.low(i), statistics.high(i));
  end
  string = sprintf('%s \\\\', string);
  disp(string)
end  

disp('\hline');

% T_ij
property = 'Tij';
for i = 1:nstates
  for j = 1:nstates
    string = sprintf('$T_{%d%d}$ & $%.3f$', i, j, true_model.Tij(i,j));
    for length_index = 1:NL
      statistics = getfield(property_statistics{length_index}, property);
      string = sprintf('%s & $%.3f_{\\:%.3f}^{\\:%.3f}$', string, statistics.mean(i,j), statistics.low(i,j), statistics.high(i,j));
    end
    string = sprintf('%s \\\\', string);
    disp(string)
  end
end  

disp('\hline');

% mu_i
property = 'mu_i';
for i = 1:nstates
  string = sprintf('$\\mu_{%d}$ & $%.3f$', i, true_model.mu_i(i));
  for length_index = 1:NL
    statistics = getfield(property_statistics{length_index}, property);
    string = sprintf('%s & $%.3f_{\\:%.3f}^{\\:%.3f}$', string, statistics.mean(i), statistics.low(i), statistics.high(i));
  end
  string = sprintf('%s \\\\', string);
  disp(string)
end  

disp('\hline');

% sigma_i
property = 'sigma_i';
for i = 1:nstates
  string = sprintf('$\\sigma_{%d}$ & $%.3f$', i, true_model.sigma_i(i));
  for length_index = 1:NL
    statistics = getfield(property_statistics{length_index}, property);
    string = sprintf('%s & $%.3f_{\\:%.3f}^{\\:%.3f}$', string, statistics.mean(i), statistics.low(i), statistics.high(i));
  end
  string = sprintf('%s \\\\', string);
  disp(string)
end  

disp('\hline');

disp('\end{tabular}');

