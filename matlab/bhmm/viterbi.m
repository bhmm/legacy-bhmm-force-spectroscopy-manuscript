function [s_t, log_trajectory_probability] = viterbi(o_t, model)
% Viterbi algorithm for computing most probable state trajectory given trajectory of observations.
%
% ARGUMENTS
%   o_t (T array) - trajectory of observations
%   model - model structure
%
% RETURNS
%   s_t (T array of int32) - most probable state trajectory (s_t(t) \in 1..model.nstates)
%   log_trajectory_probability - log trajectory probability

% Determine trajectory length.
T = length(o_t);

%
% Initialization
%

s_t = zeros(1,T,'int32'); % s_t(t) is the most likely state trajectory time t
%delta_it = zeros(model.nstates, T, 'double');
log_delta_it = zeros(model.nstates, T, 'double');
Phi_it = zeros(model.nstates, T, 'int32');

for i = 1:model.nstates
  %delta_it(i,1) = model.Pi(i) .* model.states{i}.emission_probability(o_t(1));
  log_delta_it(i,1) = model.logPi(i) + model.states{i}.log_emission_probability(o_t(1));  
  Phi_it(i,1) = 0;
end

%
% Recursion
%

for t = 2:T
  for j = 1:model.nstates
    %[valmax, argmax] = max(delta_it(:,t-1) .* model.Tij(:,j));
    %delta_it(j,t) = valmax * model.states{j}.emission_probability(o_t(t));
    
    [log_valmax, argmax] = max(log_delta_it(:,t-1) + model.logTij(:,j));
    log_delta_it(j,t) = log_valmax + model.states{j}.log_emission_probability(o_t(t));    
    
    Phi_it(j,t) = argmax;
  end
end

%
% Termination
%

%[valmax, argmax] = max(delta_it(:,T));
%trajectory_probability = valmax;
[log_valmax, argmax] = max(log_delta_it(:,T));
log_trajectory_probability = log_valmax;
s_t(T) = argmax;

%
% Reconstruction
%

for t = (T-1):-1:1
  s_t(t) = Phi_it(s_t(t+1),t+1);
end

return
