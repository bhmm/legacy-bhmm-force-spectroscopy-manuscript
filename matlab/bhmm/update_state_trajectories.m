function model = update_state_trajectories(data, model, options)
% Update state trajectories given emission probabilities and transition matrix.
%
% model = update_state_trajectories(data, model, options)
%
% ARGUMENTS
%   data (cell array of 1D arrays) - observed trajectories of some real-valued signal
%   model (structure) - model parameters, containing state trajectories
%   options (struct) - options for run - see bhmm_default_options()
%
% RETURNS
%   model (structure) - model with updated state trajectories.
%
% NOTE
%   We compute the probability of a trajectory conditioned on the first observation, so that the trajectory may not necessarily be in equilibrium.
%   A modified version of the forward-backward algorithm (due to N. S. Hinrichs) is used in the Bayesian sampling.

% Associate emission probabilities with states.
% TODO: Move this elsewhere?
for i = 1:model.nstates
  state = model.states{i};
  % Define normal density for emission probability.
  model.states{i}.emission_probability = @(o) (2*pi)^(-1/2) * state.sigma^(-1) * exp(-(1/2)*((o - state.mu) / state.sigma)^2); 
  model.states{i}.log_emission_probability = @(o) - (1/2)*log(2*pi) - log(state.sigma) - (1/2)*((o - state.mu) / state.sigma)^2;
end

% Calculate the stationary probability distribution from dominant eigenvector of transition matrix.
model.logTij = log(model.Tij); % elementwise logarithm - can we use a more numerically stable approach later?
model.logPi = log(model.Pi);

% Update state trajectories.
for trajectory_index = 1:length(data)
  % Extract trajectory of observations.
  o_t = data{trajectory_index}; % o_t(t) is the observation at time t \in 1...T.

  switch (options.updateMethod)
   case 'bayesian'
    % Sample a state trajectory.

    % Use fast Java helper code.
    mu = zeros(model.nstates,1);
    sigma = zeros(model.nstates,1);
    for i = 1:model.nstates
      mu(i) = model.states{i}.mu;
      sigma(i) = model.states{i}.sigma;
    end  
    s_t = bhmm_helper.sampleStateTrajectory(o_t, mu, sigma, log(model.Pi), log(model.Tij));

    % Use slow Matlab code.
    %s_t = sample_state_trajectory(o_t, model);
   case 'maximum-likelihood'
    % Determine the most likely state sequence using the Viterbi algorithm.
    [s_t, log_trajectory_probability] = viterbi(o_t, model);
   otherwise
    error(sprintf('options.updateMethod = "%s" is unknown', options.updateMethod));
  end
  
  % Store trajectory.
  model.state_trajectories{trajectory_index} = s_t; 
end

% Remove temporary functions from states.
for i = 1:model.nstates
  state = model.states{i};
  state = rmfield(state, 'emission_probability');
  state = rmfield(state, 'log_emission_probability');
  model.states{i} = state;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s_t = sample_state_trajectory(o_t, model)
% Sample a state trajectory given the observations and model.
%
% ARGUMENTS
%   o_t (T array) - trajectory of observations
%   model - model structure
%
% NOTES
%   We sample state trajectory S ~ P(S | O, T, E)
%     where O is observed emissions
%           T is transition matrix
%           E is emission probabilities
%
% TODO
%  * This algorithm currently assumes initial state s_1 is samples from the stationary density.  Add an option for nonequilibrium initial conditions?

% Determine trajectory length.
T = length(o_t);

%
% Forward part
%

log_alpha_it = zeros(model.nstates, T, 'double'); 

for i = 1:model.nstates
  log_alpha_it(i,1) = model.logPi(i) + model.states{i}.log_emission_probability(o_t(1));
end

for t = 2:T
  for j = 1:model.nstates
    log_alpha_it(j,t) = logsum(log_alpha_it(:,t-1) + model.logTij(:,j)) + model.states{j}.log_emission_probability(o_t(t));    
  end
end

%
% Sample trajectory
% 

s_t = zeros(1,T,'int32'); % state trajectory - s_t(t) is state at time t

log_p_i = log_alpha_it(:,T);
p_i = exp(log_p_i - logsum(log_alpha_it(:,T)));
s_t(T) = draw(p_i);

for t = (T-1):-1:1
  % Compute P(s_t = i | s_{t+1}..s_T).
  log_p_i = log_alpha_it(:,t) + model.logTij(:,s_t(t+1));  
  p_i = exp(log_p_i - logsum(log_p_i));
  
  % Draw from this distribution.
  s_t(t) = draw(p_i);
end

%show_trajectory(s_t)

return
