function plot_correlated_model_parameters(models, options)
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

properties = {'Pi', 'Tij', 'Kij', 'tau_i', 'mu_i', 'sigma_i'};
nproperties = length(properties);
% Compute model parameters.
model_parameters = struct()
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

  model_parameters = setfield(model_parameters, property, x_n);

end

% Generate plot.
clf;
%plot(squeeze(model_parameters.Kij(:,3,1)), squeeze(model_parameters.Kij(:,3,2)), '.');
%plot(squeeze(model_parameters.Pi(:,1)), squeeze(model_parameters.Pi(:,3)), '.');
%plot(squeeze(model_parameters.mu_i(:,1)), squeeze(model_parameters.mu_i(:,2)), '.');
%plot(squeeze(model_parameters.sigma_i(:,1)), squeeze(model_parameters.sigma_i(:,3)), '.');
%plot(squeeze(model_parameters.mu_i(:,2)), squeeze(model_parameters.Pi(:,2)), '.');
%plot(squeeze(model_parameters.mu_i(:,2)), squeeze(model_parameters.sigma_i(:,2)), '.');
%plot(squeeze(model_parameters.Kij(:,2,3)), squeeze(model_parameters.sigma_i(:,2)), '.');
hist(squeeze(model_parameters.Kij(:,3,2)) ./ squeeze(model_parameters.Kij(:,3,1)),40);
%plot(squeeze(model_parameters.Kij(:,3,2)) ./ squeeze(model_parameters.Kij(:,3,1)), squeeze(model_parameters.Kij(:,1,2)) ./ squeeze(model_parameters.Kij(:,1,3)), '.');
axis square;

return
