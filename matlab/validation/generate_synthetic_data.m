function data = generate_synthetic_data(model, trajectory_lengths)
% Generate synthetic test data from the given model.
%
%  [data] = generate_synthetic_data(model, trajectory_lengths)
%
% ARGUMENTS
%  model - the model to use to generate synthetic data
%  trajectory_lengths - vector of trajectory lengths
%
% RETURNS
%  data  - the synthetic test data generated from the model
%
% NOTE
%  The random number seed is set deterministically.

% Get number of trajectories.
ntrajectories = length(trajectory_lengths);

% Allocate storage for data.
data = cell(ntrajectories,1);

% Generate synthetic trajectories.
for trajectory_index = 1:ntrajectories
  % Get trajectory length.
  T = trajectory_lengths(trajectory_index);

  % Generate state trajectory.
  s_t = zeros(1,T,'int32'); % s_t(t) is state visited at time t
  s_t(1) = draw(model.Pi);
  for t = 2:T
    s_t(t) = draw(model.Tij(s_t(t-1),:));
  end
  model.state_trajectories{1} = s_t; 

  % Generate observation trajectory.
  o_t = zeros(1,T,'single'); % o_t(t) is observation at time t
  for t = 1:T
    state = model.states{s_t(t)};
    o_t(t) = randn() * state.sigma + state.mu;
  end

  % Sotre observation trajectory.
  data{trajectory_index} = o_t;
end

return
