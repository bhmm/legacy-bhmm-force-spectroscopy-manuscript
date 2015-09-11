function model = generate_segmentation_model(data, nstates, options)
% Construct model based on segmentation of the data using Gaussian mixture model EM fit.
%
% model = generate_segmentation_model(data, nstates, options, [initial_model])
%
% ARGUMENTS
%   data (cell array of 1D arrays) - observed trajectories of some real-valued signal (e.g. trap extension, FRET efficiency)
%   nstates (int) - number of states to fit
%   options - options structure to control model generation
%
% RETURNS
%   model (structure) - segmentation model
%
% NOTES
%   Trajectories provided in 'data' cell array may be of different lengths. 
%
% TODO
%   Add ability for user to pass in optional 'options' structure.

% Import fast Java helper code.
import bhmm_helper;

% Create initial emission model from MLE of Gaussian mixture.
model = em_gaussian_mixture(data, nstates, options);

% Find points where Gaussian probabilities cross.
divisions = zeros(model.nstates-1, 1);
for i = 1:(nstates-1)
  mu1 = model.states{i}.mu;
  mu2 = model.states{i+1}.mu;
  sigma1 = model.states{i}.sigma;
  sigma2 = model.states{i+1}.sigma;
  % Find intersection.
  f = @(x) (-log(sigma1)-0.5*((x-mu1)/sigma1).^2) - (-log(sigma2)-0.5*((x-mu2)/sigma2).^2);
  alpha = sigma1^(-2) / (sigma1^(-2) + sigma2^(-2));
  beta =  sigma2^(-2) / (sigma1^(-2) + sigma2^(-2));
  [xcross, fval] = fzero(f, alpha*mu1 + beta*mu2);
  divisions(i) = xcross;
end
divisions

% Determine state trajectories.
model.state_trajectories = {};
for trajectory_index = 1:length(data)
  o_t = data{trajectory_index}; % o_t(t) is the observation at time t \in 1...T.
  s_t = ones(size(o_t));
  for i = 1:(nstates-1)
    indices = find(o_t >= divisions(i));
    s_t(indices) = i+1;
  end
  model.state_trajectories{trajectory_index} = s_t;
end

% Compute new state properties.
for i = 1:nstates
  o_n = [];
  for trajectory_index = 1:length(data)    
    o_t = data{trajectory_index};
    s_t = model.state_trajectories{trajectory_index};
    indices = find(s_t == i);
    o_n = [o_n o_t(indices)];
  end
  model.states{i}.mu = mean(o_n);
  model.states{i}.sigma = std(o_n);
end
clear s_t o_n;

% Compute empirical transition matrix.
Nij = zeros(nstates, nstates); % Nij(i,j) is fractional expected transition counts from i to j
% Construct transition count matrix from state trajectories.
for trajectory_index = 1:length(data)
  s_t = model.state_trajectories{trajectory_index};
  for i = 1:model.nstates
    for j = 1:model.nstates
      indices = find(s_t(1:end-1)==i & s_t(2:end)==j);
      Nij(i,j) = Nij(i,j) + length(indices);
    end
  end
end
% Update transition matrix estimate.
model.Tij = transition_matrix_mle(Nij, options);

% Compute stationary probability.
model.Pi = stationary_probability(model.Tij);
model.logPi = log(model.Pi);

return


