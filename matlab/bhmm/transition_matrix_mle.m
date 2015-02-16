function Tij = transition_matrix_mle(cij, options)
% Generate maximum-likelihood estimate of row-stochastic transition matrix using iterative procedure of Ref. [1].
%
% WARNING
%  There is currently no constraint ensuring sum(Pi(1:N-1)) < 1 during likelihood maximization.  This needs to be fixed before this code is usable.
%
% Tij = transition_matrix_mle(cij, options)
%
% ARGUMENTS
%  cij (MxM matrix) - transition conditional count matrix (row form), including any prior pseudocounts
%    cij(i,j) is number of observed counts terminating in j at time tau given that they initially started in i at time zero.
%  options - optional options structure
%    if options.reversible is true, then reversible transition matrix sampling will be used 
%    if options.diagonally_dominant is true, then only diagonally-dominant transition matrices will be accepted
%    if options.verbosity sets level of verbose output (0 = none, 1 = summary report, 2 = output each iteration)
%    if options.equilibrium is true, trajectories are initially in equilibrium
%
% RETURN VALUES
%  Tij (MxM matrix) - maximum likelihood row-stochastic transition matrix
%
% REFERENCES
%
% [1] Prinz JH, Wu H, Sarich M, Keller B, Fischbach M, Held M, Chodera JD, Schuette C, and Noe F. 
% Markov models of molecular kinetics: Generation and validation. Submitted.
%
% TODO
%  * Enforce positivity in xij elements.

% Parameters.
MAXITS = 1000; % maximum number of iterations
RELATIVE_TOLERANCE = 1.0e-5; % relative convergence tolerance
EPSILON = 1.0e-3; % fractional counts to add to each transition to avoid numerical issues

if (options.verbosity >= 2)
  disp('Computing maximum-likelihood estimate of transition matrix');
  disp('cij = ');
  disp(cij);
end

% Determine size of transition count matrix.
M = size(cij,1);

% Cast count matrix to doubles.
cij = double(cij);

if (~options.reversible)
  % Dispatch non-reversible estimate.
  if (options.verbosity >= 2)
    disp('Using non-reversible estimate of transition matrix.');
  end
  Tij = cij;
  for i = 1:M
    Tij(i,:) = cij(i,:) / sum(cij(i,:));
  end
  return;
end

%% Initialize.

% Add small EPSILON to achieve numerical stability, if necessary.
if any(cij < EPSILON)
  cij = cij + EPSILON;
end

ci = zeros(1,M);
for i = 1:M
  ci(i) = sum(cij(i,:));
end

xij = cij + cij';
xi = zeros(1,M);
for i = 1:M
  xi(i) = sum(xij(i,:));
end

%% Iterate

for iteration = 1:MAXITS
  xij_old = xij;

  % update step 1.1
  for i = 1:M
    xi(i) = cij(i,i) * (xi(i) - xij(i,i)) / (ci(i) - cij(i,i));
  end
  for i = 1:M
    xi(i) = sum(xij(i,:));
  end
  
  % update step 1.2
  for i = 1:(M-1)
    for j = (i+1):M
      a = (ci(i) - cij(i,j) + ci(j) - cij(j,i));
      b = ci(i)*(xi(j) - xij(j,i)) + ci(j)*(xi(i) - xij(i,j)) - (cij(i,j) + cij(j,i))*(xi(i) - xij(i,j) + xi(j) - xij(j,i));
      c = - (cij(i,j) + cij(j,i))*(xi(i)-xij(i,j))*(xi(j)-xij(j,i));
      xij(i,j) = (-b + sqrt(b^2 - 4*a*c)) / (2*a);
      xij(j,i) = xij(i,j);
      for i = 1:M
	xi(i) = sum(xij(i,:));
      end
    end
  end

  % TODO: Check convergence.
  delta = norm(xij_old - xij, 'fro') / norm(xij, 'fro');
  if (delta < RELATIVE_TOLERANCE)
    if (options.verbosity >= 2)
      disp(sprintf('Converged to relative tolerance of %e in %d iterations.', delta, iteration));
    end
    break;
  end
end

%% step 2: Compute Tij from xij:
% TODO: Check this.
Tij = zeros(M,M);
for i = 1:M
  Tij(i,:) = xij(i,:) / xi(i); 
end

if (options.verbosity >= 2)
  disp('Final transition matrix');
  disp(Tij);
end

return
