function Tij = transition_matrix_sample(Tij, Nij, options)
% Generate an updated transition matrix sample T ~ P(T | N).
%
% Tij = sample_transition_matrix(Tij, Nij, tau)
%
% ARGUMENTS
%  Tij (MxM matrix) - current transition matrix sample 
%    This transition matrix must be reversible.
%  Nij (MxM matrix) - transition conditional count matrix
%    Nij(i,j) is number of observed counts terminating in j at time tau given that they initially started in i at time zero.
%  options - optional options structure
%    if options.reversible is true, then reversible transition matrix sampling will be used 
%    if options.diagonally_dominant is true, then only diagonally-dominant transition matrices will be accepted
%    if options.verbosity sets level of verbose output (0 = none, 1 = summary report, 2 = output each iteration)
%    if options.equilibrium is true, trajectories are initially in equilibrium
%
% RETURN VALUES
%  Tij_new (MxM matrix) - updated sample of transition matrix, may be correlated with previous sample
%
% NOTES
%  If options.reversible = true, the transition matrix is guaranteed to satisfy detailed balance.
%  The transition matrix will not necessarily correspond to a rate matrix.
%
%  The algorithm described in [1] and modified in [2] (especially with regard to the prior) is implemented.
%  Transitions between all states are allowed -- no constraints on the sparisty of Kij are imposed.

if (options.reversible)
  Tij = transition_matrix_sample_reversible(Tij, Nij, options);
else
  Tij = transition_matrix_sample_nonreversible(Tij, Nij, options);
end

return
