function model = bayesian_gaussian_mixture(data, nstates, options)
% Fit a Gaussian mixture model to observed data using the EM algorithm using random initialization.
%
% model = bayesian_gaussian_mixture(data, nstates, options)
%
% ARGUMENTS
%   data (cell array of 1D arrays) - observed trajectories of some real-valued signal
%   nstates - number of Gaussian components to fit
%   options - options data structure
%
% RETURNS
%   model (structure) - initial model parameters
%
% NOTES
%
%
% REFERENCES
%   

if (options.verbosity >= 1)
  disp('******************************************************************************');
  disp('Gaussian mixture model fitting by expectation maximization (EM) algorithm');
end

% DEBUG
options.enforce_state_ordering = true; % enforce ordering of states at this stage
RELATIVE_TOLERANCE = options.convergenceTolerance;
  
% Create model structure.
model = struct();

% Store number of states.
model.nstates = nstates;

% Determine number of trajectories.
nTrajectories = length(data);

% Concatenate all observations.
observations = [];
for trajectoryIndex = 1:nTrajectories
  % Extract trajectory.
  o_t = data{trajectoryIndex};
  % Make into a row vector.
  o_t = reshape(o_t, 1, length(o_t));
  % Concatenate observations.
  observations = [observations o_t];  
end
nobservations = length(observations);

%
% Generate initial guess by randomly sampling from the data
%

nobservations;
random_permutation = randperm(nobservations);
random_centers = random_permutation(1:nstates);
distances = abs(observations'*ones(1,nstates) - ones(nobservations,1)*observations(random_centers))
[mindistance, initialstatelabels] = min(distances')

% Allocate storage for model parameters.
mu = zeros(nstates,1); % mu(i) is the Gaussian mean for state i
sigma = zeros(nstates,1); % sigma(i) is the Guassian standard deviation for state i
Pi = zeros(nstates,1); % Pi(i) is the normalized weight for Gaussian state i

% Determine mean and std dev of each cluster.
for i = 1:nstates
  indices = find(initialstatelabels == i);
  mu(i)    = mean(observations(indices)); % mean
  sigma(i) = std(observations(indices));  % standard deviation
  Pi(i)    = length(indices) / double(length(observations)); % normalized weight
end

if (options.verbosity >= 1)
  disp('initial guess:');
  for i = 1:nstates
    disp(sprintf('state %5d : Pi = %8.6f mu = %8.3f  sigma = %8.3f', i, Pi(i), mu(i), sigma(i)));
  end
  disp('******************************************************************************');
end

% DEBUG: Initialize for example to give segmentation model best chance of working.
%mu(1) = 3;
%mu(2) = 4.7;
%mu(3) = 5.6;
%sigma(1) = 1.0;
%sigma(2) = 0.3;
%sigma(3) = 0.2;

%
% Carry out EM iterations.
%

niterations = 1000; % TODO: Set this to some large maximum once convergence assessment is added.

for iteration = 1:niterations
  if (options.verbosity >= 2)
    disp(sprintf('EM iteration %d', iteration));
  end

  % Store old model.
  old_mu = mu;
  old_sigma = sigma;

  % Compute p(l | x_i, mu, sigma) for each sample.
  p_ti = bhmm_helper.computeStateProbabilities(observations, mu, sigma, Pi);

  % Update Gaussian parameters.
  for i = 1:model.nstates
    Pi(i) = mean(p_ti(:,i));
    mu(i) = sum(p_ti(:,i)' .* observations) / sum(p_ti(:,i));
    sigma(i) = sqrt(sum(p_ti(:,i)' .* (observations - mu(i)).^2) / sum(p_ti(:,i)));
  end

  % Reorder states if necessary.
  if options.enforce_state_ordering
    % Sort by state means.
    [sorted_mu, indices] = sort(mu);
    % Permute states.
    Pi = Pi(indices);
    mu = mu(indices);
    sigma = sigma(indices);
  end
    
  if (options.verbosity >= 2)
    for i = 1:nstates
      disp(sprintf('state %5d : Pi = %8.6f mu = %8.3f  sigma = %8.3f', i, Pi(i), mu(i), sigma(i)));
    end
    likelihood = compute_gmm_log_likelihood(mu, sigma, log(Pi), data);
    disp(sprintf('Log Likelihood %8.6f:', likelihood));
    disp('******************************************************************************');
  end

  % Assess convergence.
  if (norm(mu - old_mu) / norm(mu) < RELATIVE_TOLERANCE) && (norm(sigma - old_sigma) / norm(sigma) < RELATIVE_TOLERANCE)
    if (options.verbosity >= 2)
      disp(sprintf('state means and variances converged to specified relative tolerance %e', RELATIVE_TOLERANCE));
    end    
    break;
  end
end

if (options.verbosity >= 1)
  disp('Converged results:');
  for i = 1:nstates
    disp(sprintf('state %5d : Pi = %8.6f mu = %8.3f  sigma = %8.3f', i, Pi(i), mu(i), sigma(i)));
  end
  disp('******************************************************************************');
end

% Store in model.
model.states = cell(nstates,1);
model.Pi = Pi;
for i = 1:nstates
  % Compute parameters that define emission probabilities of this state.
  state = struct();
  state.mu = mu(i);
  state.sigma = sigma(i);
  % Store this state.
  model.states{i} = state;  
end

return
