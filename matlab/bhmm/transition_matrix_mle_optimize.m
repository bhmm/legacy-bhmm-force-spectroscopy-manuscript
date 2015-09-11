function Tij = transition_matrix_mle_optimize(Nij, options)
% Generate maximum-likelihood estimate of row-stochastic transition matrix by optimization procedure.
%
% WARNING
%  There is currently no constraint ensuring sum(Pi(1:N-1)) < 1 during likelihood maximization.  This needs to be fixed before this code is usable.
%
% Tij = transition_matrix_mle(Nij, options)
%
% ARGUMENTS
%  Nij (MxM matrix) - transition conditional count matrix (row form)
%    Nij(i,j) is number of observed counts terminating in j at time tau given that they initially started in i at time zero.
%  options - optional options structure
%    if options.reversible is true, then reversible transition matrix sampling will be used 
%    if options.diagonally_dominant is true, then only diagonally-dominant transition matrices will be accepted
%    if options.verbosity sets level of verbose output (0 = none, 1 = summary report, 2 = output each iteration)
%    if options.equilibrium is true, trajectories are initially in equilibrium
%
% RETURN VALUES
%  Tij (MxM matrix) - maximum likelihood row-stochastic transition matrix
%
% NOTES
%
%  Here, we use the representation of Diaconis, where the element xij(i,j) represents the probability of observing a transition
% from either i to j or j to i at equilibrium.  Note that this is not a *conditional* probability, like the transition matrix
% element Tij, but rather an unconditioned probability of observing this transition, such that
%
%    sum_{i >= j} xij(i,j) = 1
%
% The transition matrix Tij is determined from
%
%    Tij(i,j) = xij(i,j) / sum(xij(i,:))
%
% The likelihood function is
%
%    L = \prod_{i,j} Tij(i,j)^Nij(i,j)
%
% where the log-likelihood is
%
%    logL = \sum_{i,j} Nij \log Tij(i,j)
%
% REFERENCES
%
% TODO
%  * Enforce positivity in xij elements.

if (options.verbosity >= 2)
  disp('Computing maximum-likelihood estimate of transition matrix');
  disp('Nij = ');
  disp(Nij);
end

% Determine size of transition matrix.
M = size(Nij,1);

% Cast count matrix to doubles.
Nij = double(Nij);

% Construct an initial rate matrix that satisfies detailed balance.
N = sum(sum(Nij)); % total number of observed transitions
xij = (Nij + Nij' - diag(diag(Nij))) / N; % xij(i,j) is fraction of total counts observed in transition from i to j or j to i

% TODO: If requested, modify initial xij to ensure self-transitions have probability of at least 0.5.

% Pack into a vector for use in optimization.
x0 = xij2x(xij);

% Find maximum likelihood rate estiamte.
% TODO: Use some scheme to guarantee xij >= 0 and, if requested, xij(i,i) >= 0.5.
minimization_options = optimset;
minimization_options.MaxFunEvals = 50000;
minimization_options.MaxIter = 50000;
xmin = fminsearch(@(x) objective(x, Nij, options), x0, minimization_options);

% Unpack vector.
xij = x2xij(xmin);

% Compute transition matrix.
Tij = zeros(M,M);
for i = 1:M
  Tij(i,:) = xij(i,:) / sum(xij(i,:));
end

if (options.verbosity >= 2)
  disp('Final transition matrix');
  disp(Tij);
end

return

%%%%%%%%%%%%%%%%%%

function x = xij2x(xij)
% Pack unique elements of xij into vector x.
%
% x = xij2x(xij)

M = size(xij,1);

x = zeros(M*(M+1)/2, 1);
index = 1;
for m = 1:M
  x(index:(index+M-m)) = xij(m,m:end);
  index = index + M-m+1;
end

return 

%%%%%%%%%%%%%%%%%%

function xij = x2xij(x)

N = length(x);
M = round(sqrt(1/4 + 2*N) - 1/2);
xij = zeros(M,M);
index = 1;
for m = 1:M
  xij(m,m:end) = x(index:(index+M-m));
  index = index + M-m+1;
end
xij = xij + triu(xij,1)';

return

%%%%%%%%%%%%%%%%%%

function obj = objective(x, Nij, options)
% Objective function used in likelihood maximization.
%
% ARGUMENTS
%  x (M-1 + M(M-1)/2 vector) - vector containing unique elements of matrix X = (x_ij) in packed form (see pack_vector and unpack_vector)
%  Nij (MxM) - transition counts - Nij(i,j) is number of times system was observed in state j a time tau after it was given to be in state i
%  options - options structure
%
% RETURN VALUE
%  obj - negative log likelihood, to be minimized in optimization scheme

% Determine number of states.
M = size(Nij,1);

% Reject if x < 0.
if (any(x < 0))
  LARGE = 1.0e10;
  obj = LARGE;
  if (options.verbosity >= 3)
    disp('x = ');
    disp(x);
    disp('rejecting');
  end
  return
end

% Normalize x.
x = x / sum(x);

% Unpack vector containing symmetrix matrix of i to j transition probabilities.
% xij(i,j) is the equilibrium probability of observing the system in state i and then state j the following time.
xij = x2xij(x);

% Compute transition matrix.
Tij = zeros(M,M);
for i = 1:M
  Tij(i,:) = xij(i,:) / sum(xij(i,:));
end

% Compute logarithm of elements of x_ij.
log_Tij = log(Tij);

% Compute log likelihood.
logL = sum(sum(Nij .* log_Tij));

% Compute objective.
obj = - logL;

if (options.verbosity >= 3)
  disp('x = ');
  disp(x);
  disp(sprintf('objective = %f', obj));
end

return
