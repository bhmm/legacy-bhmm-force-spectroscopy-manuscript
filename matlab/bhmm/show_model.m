function show_model(model)
% Print the parameters of the specified model.
%
% show_model(model)
%
% ARGUMENTS
%  model (struct) - model to be shown

% Display transition matrix.
disp('Tij = ');
disp(model.Tij);

% Display states.
for stateIndex = 1:model.nstates
  state = model.states{stateIndex};
  disp(sprintf('state %4d : Pi = %8.6f  mu = %12.7f  sigma = %12.7f', stateIndex, model.Pi(stateIndex), state.mu, state.sigma));
end

return
