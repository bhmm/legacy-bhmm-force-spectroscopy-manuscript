function [log_xi_tij, log_gamma_ti] = baum_welch(o_t, model, options)
% Baum-Welch algorithm for computing expected transition and state observation counts.
% 
% ARGUMENTS
%   o_t (T) [real array] - trajectory of real-valued observations
%   model [structure] - model structure
%   options [structure] - options structure
%
% RETURNS
%   log_xi_tij (T,N,N) [real array] - log-probability of observing transition  - xi_tij(t,i,j) = P(X_t = i, X_{t+1} = j | O, model)
%   log_gamma_ti (T,N) [real array] - log-probability observing state i at time t - gamma_ti(t,i) = P(X_t = i | O, model)

% Constants.
T = length(o_t); % number of observations in this trajectory
N = model.nstates; % number of states

if options.equilibrium
  % Experiment starts in equilibrium.
  log_rho_i = model.logPi;
elseif isfield(options, 'log_rho_i')
  % Use specified nonequilibrium starting conditions.
  log_rho_i = options.log_rho_i;
else
  % Use uniform initial starting conditions, which should be equivalent to complete ignorance.
  log_rho_i = model.logPi * 0.0;
  log_rho_i = log_rho_i - log(sum(exp(log_rho_i)));
end

% Use Java code.
if options.use_java
  % Prepare data for fast Java helper code.
  mu = zeros(model.nstates,1);
  sigma = zeros(model.nstates,1);
  for i = 1:model.nstates
    mu(i) = model.states{i}.mu;
    sigma(i) = model.states{i}.sigma;
  end  
  % Update model parameters.
  model.logTij = log(model.Tij); % elementwise logarithm - can we use a more numerically stable approach later?
  model.logPi = log(model.Pi);

  % Perform Baum-Welch.
  % log-probability of observing transition i to j at time t : xi_tij(t,i,j) = P(X_t = i, X_{t+1} = j | O, model)
  log_alpha_ti = bhmm_helper.forwardAlgorithm(o_t, mu, sigma, model.logPi, model.logTij, log_rho_i);  
  log_beta_ti = bhmm_helper.backwardAlgorithm(o_t, mu, sigma, model.logPi, model.logTij, log_rho_i);  
  log_xi_tij = bhmm_helper.baumWelch_xi(o_t, mu, sigma, model.logPi, model.logTij, log_alpha_ti, log_beta_ti);
  log_gamma_ti = bhmm_helper.baumWelch_gamma(o_t, mu, sigma, model.logPi, model.logTij, log_alpha_ti, log_beta_ti);
  return
end

% Compute forward-backward parameters.
[log_alpha_ti, log_beta_ti] = forward_backward_algorithm(o_t, model, log_rho_i);

% First, compute log-probability of observation sequence o_t given model (used subsequently as a normalizing constant)
log_O = logsum(log_alpha_ti(T,:));

% Compute desired log-probabilities by Baum-Welch.
% TODO This could potentially be speeded up by using vector notation
log_xi_tij = zeros(T,N,N); % log-probability of observing transition i to j at time t : xi_tij(t,i,j) = P(X_t = i, X_{t+1} = j | O, model)
log_gamma_ti = zeros(T,N); % log-probability observing state i at time t : gamma_ti(t,i) = P(X_t = i | O, model)

for t = 1:(T-1)
  for i = 1:N
    for j = 1:N
      log_xi_tij(t,i,j) = log_alpha_ti(t,i) + model.states{i}.log_emission_probability(o_t(t)) + model.logTij(i,j) + log_beta_ti(t+1,j) - log_O;
    end
    log_gamma_ti(t,i) = logsum(log_xi_tij(t,i,:));
  end
end  

return
