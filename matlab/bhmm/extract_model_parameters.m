function [Tij, mu, sigma] = extract_model_parameters(model)
% Extract model parameters, returning them as a list of matrices and vectors.
%
% [Tij, mu, sigma] = extract_model_parameters(model)

Tij = model.Tij;
nstates = model.nstates;
mu = zeros(nstates,1);
sigma = zeros(nstates,1);
for i = 1:nstates
  mu(i) = model.states{i}.mu;
  sigma(i) = model.states{i}.sigma;
end

return

  