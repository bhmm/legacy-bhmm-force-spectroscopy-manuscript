function Tij = generate_transition_matrix_guess(Nij, tau)
% Generate an initial row-stochastic transition matrix from the specified count matrix.
% The transition matrix is guaranteed to satisfy
%
% Tij = exp(Kij * tau)
%
% where Kij is some rate matrix that satisfies:
%
% Kij(i,j) > 0 for i \ne j
% sum(K(i,:)) = 0 for all i
%
% The maximum-likelihood rate matrix according to Gerhard Hummer's formalism [1,2] is used to determine Kij.
%
% ARGUMENTS
%  Nij (MxM matrix) - transition conditional count matrix (row form)
%    Nij(i,j) is number of observed counts terminating in j at time tau given that they initially started in i at time zero.
%  tau - evolution time to relate rate matrix to transition matrix through Tij = expm(Kij * tau)
%
% RETURN VALUES
%  Tij (MxM matrix) - maximum likelihood row-stochastic transition matrix
%
% REFERENCES
%  [1] Sriraman S, Kevrekidis IG, and Hummer G. Coarse master equation from Bayesian analysis of replica molecular dynamics simulations. JPC B 109(14):6479-6484, 2005.
%  [2] Buchete N-V and Hummer G. Coarse master equations for peptide folding dynamics. JPC B 112(19):6057-6069, 2008.
%
% TODO
%  Have pack_vector and unpack_vector convert directly between unique log probabilities and rate elements and transition matrix?
%  Work with rate matrices here instead of transition matrices?

% Determine size of transition matrix.
M = size(Nij,1);

% Construct an initial rate matrix that satisfies detailed balance.
% We add a pseudocount to each transition (to avoid divergent rate estimates), 
% symmetrize the transition matrix by adding the transpose (which screws things up a lot for transition counts not taken from global equilibrium), 
% and estimate the short-time generator.  This should at least provide a crude valid initial rate matrix estimate that will be improved
% by subsequent optimization step.
Tij = zeros(M,M,'double');
Nij = double(Nij) + ones(M,M);
for i = 1:M
  for j = 1:M
    Tij(i,j) = (Nij(i,j) + Nij(j,i)) / (sum(Nij(i,:)) + sum(Nij(:,i)));
  end  
end
% Estimate short-time generator from transition matrix.
Kij = Tij / tau;
for i = 1:M
  Kij(i,i) = 0.0;
  Kij(i,i) = - sum(Kij(i,:));
end
% Exponentiate to obtain a transition matrix corresponding to positive rates.
Tij = expm(Kij * tau);

% Pack into a vector for use in optimization.
x0 = pack_vector(Tij, tau);

% Find maximum likelihood rate estiamte.
options = optimset;
options.MaxFunEvals = 50000;
options.MaxIter = 50000;
xmin = fminsearch(@(x) objective(x, Nij, tau), x0, options);

% Unpack vector.
Tij = unpack_vector(xmin, tau);

return

%%%%%%%%%%%%%%%%%%

function obj = objective(x, Nij, tau)
% Objective function used in likelihood maximization.
%
% ARGUMENTS
%  x (M-1 + M(M-1)/2 vector) - vector containing unique elements of logPi and logKij in packed form (see pack_vector and unpack_vector)
%  Nij (MxM) - transition counts - Nij(i,j) is number of times system was observed in state j a time tau after it was given to be in state i
%  tau - lag time for transition counts
%
% RETURN VALUE
%  obj - negative log likelihood, to be minimized in optimization scheme

% Determine number of states.
M = size(Nij,1);

% Unpack vector containing unique log probabilities and rates.
Tij = unpack_vector(x, tau);

% Compute logarithm of elements of transition matrix.
logTij = log(Tij);

% Compute log likelihood.
logL = sum(sum(Nij .* logTij));

% Compute objective.
obj = - logL;

return
