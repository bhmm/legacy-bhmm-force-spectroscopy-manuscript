function model = update_emission_probabilities(data, model, options)
% Update state emission probability densities given state trajectories.
%
% ARGUMENTS
%   data (cell array of 1D arrays) - observed trajectories of some real-valued signal
%   model (structure) - model parameters, containing state trajectories
%
% RETURNS
%   model (structure) - model with updated state emission probability densities

% Store old model in case we need to reject state observable update.
old_model = model; 

% Update emission parameters for each state.
for state_index = 1:model.nstates
  % Extract state for convenience.
  state = model.states{state_index};
  
  % Collect all samples assigned to this state.
  observations = [];
  for trajectory_index = 1:length(data)    
    % Identify indices of trajectory samples within this state.
    indices = find(model.state_trajectories{trajectory_index} == state_index);
    % Collect trajectory samples and append them to current list.
    observations = [observations data{trajectory_index}(indices)];
  end
  if (options.verbosity >= 4)
    disp(sprintf('state %3d: %8d observations', state_index, length(observations)));
  end

  % Update emission probability distribution parameters for this state.
  switch (options.updateMethod)
   case 'bayesian'
    % Update estimate of emission probability density for this state.
    % TODO: Sample from posterior distribution.  Can we sample from joint posterior?
    N = length(observations); % number of samples in this state
    
    if (N <= 1)
      % don't update empty states
      continue;
    end
    
    % Sample new mu.
    state.mu = randn()*state.sigma/sqrt(N) + mean(observations);    

    % Slow implementation of chi-squared distribution for those Matlab installations without stats toolbox.
    % This scheme uses the improper Jeffreys prior on sigma^2, P(mu, sigma^2) \propto 1/sigma
    if (N == 1)
      chisquared = 0;
    else
      chisquared = sum(randn(N-1,1).^2);    
    end
    sigmahat2 = mean((observations - state.mu).^2);
    state.sigma = sqrt(sigmahat2) / sqrt(chisquared / N);
   
   case 'maximum-likelihood'
    % Determine sample statistics.
    state.mu = mean(observations);
    state.sigma = std(observations);
    
   otherwise
    error(sprintf('options.updateMethod = "%s" is unknown', options.updateMethod));    
  end  
  model.states{state_index} = state; % store update
end

%% Enforce any state emission probability restrictions, rejecting state update if any fail.

% Enforce state ordering so that mu are in ascending order.
if (options.enforce_state_ordering == true)
  for i = 1:(model.nstates-1)
    if (model.states{i}.mu > model.states{i+1}.mu)
      disp('State ordering violated; rejecting.');
      model = old_model;
      return
    end    
  end
end

% Enforce separation between state mu.
if (options.enforce_state_separation > 0)
  for i = 1:(model.nstates-1)
    state_separation = abs(model.states{i}.mu - model.states{i+1}.mu);
    if (state_separation < options.enforce_state_separation)
      disp(sprintf('State separation between states %d and %d is %f > %f; rejecting.', i, i+1, state_separation, options.enforce_state_separation));
      model = old_model;
      return
    end    
  end
end

return
