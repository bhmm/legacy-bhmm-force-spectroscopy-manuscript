function Tij = transition_matrix_sample_nonreversible(Tij, Nij, options)
% Generate an updated nonreversible transition matrix sample T ~ P(T | N) using Dirichlet sampling.
%
% Tij = sample_transition_matrix(Tij, Nij, tau)
%
% ARGUMENTS
%  Tij (MxM matrix) - current transition matrix sample 
%    This transition matrix must be reversible.
%  Nij (MxM matrix) - transition conditional count matrix
%    Nij(i,j) is number of observed counts terminating in j at time tau given that they initially started in i at time zero.
%  options - optional options structure
%    if options.diagonally_dominant is true, then only diagonally-dominant transition matrices will be accepted
%    if options.verbosity sets level of verbose output (0 = none, 1 = summary report, 2 = output each iteration)
%    if options.equilibrium is true, trajectories are initially in equilibrium
%    if options.prior_transition_pseudocounts exists (a KxK matrix where the (i,j) element represents transitions from i to j), these pseudocounts are included; all ones is uniform prior
%
% RETURN VALUES
%  Tij_new (MxM matrix) - updated sample of nonreversible transition matrix, may be correlated with previous sample
%
% NOTES
%  If options.reversible = true, the transition matrix is guaranteed to satisfy detailed balance.
%  The transition matrix will not necessarily correspond to a rate matrix.
%
%  The algorithm described in [1] and modified in [2] (especially with regard to the prior) is implemented.
%  Transitions between all states are allowed -- no constraints on the sparisty of Kij are imposed.
% REFERENCES
%  [1] Singhal

% PARAMETERS
ROW_SUM_TOLERANCE = 1.0e-3; % row sum tolerance

% Include prior transition pseudocounts.
if isfield(options, 'prior_transition_pseudocounts')
  Nij = double(Nij) + double(options.prior_transition_pseudocounts);
else
  % DEBUG: Add pseudocounts.
  Nij = double(Nij) + 1.0;  
end

% Determine size of transition matrix.
M = size(Tij,1);

% Promote Nij to double.
Nij = double(Nij);

% Compute total number of observed transitions out of each state.
Ni = sum(Nij,2);

%% Sanity checks.

% Check initial transition matrix does not contain negative values.
if (any(any(Tij < 0.0)))
  step
  Tij
  error('Tij has negative elements.')
end

% Check initial row sums.
row_sums = sum(Tij,2);
if any(abs(row_sums - 1) > ROW_SUM_TOLERANCE)
  step
  Tij
  error('Tij has some row sums that differ from unity.')
end

% Sample until we have an accepted sample.
accepted = false;
while (accepted == false)  
  %% Sample rows from Dirichlet posterior.
  Tij_new = Tij;
  for i = 1:M
    Tij_new(i,:) = gamrnd(Nij(i,:), 1);
    Tij_new(i,:) = Tij_new(i,:) / sum(Tij_new(i,:));
  end

  %% Reject disallowed transition matrices.
  accepted = true;

  % Check for diagonal dominance.
  if (options.diagonally_dominant && any(diag(Tij_new) < 0.5))
    accepted = false;
    continue;
  end

  % Reject transition matrices that have lost row normalization.
  if (any(abs(sum(Tij_new,2) - ones(M,1)) > ROW_SUM_TOLERANCE))
    accepted = false;
    continue;
  end

  % Reject transition matrices containing any negative values.
  if (any(any(Tij_new < 0.0)))
    accepted = false;
    continue;
  end
end

% Transition matrix has been accepted.
Tij = Tij_new;

return
