function model = rate_matrix_sample(model, Nij, mode)
% Generate an uncorrelated transition matrix sample using the method of Gerhard Hummer.
%
% Tij = sample_transition_matrix(Tij, Nij, tau)
%
% ARGUMENTS
%  Tij (MxM matrix) - current transition matrix sample 
%    This transition matrix must correspond to a valid rate matrix, such that Tij = expm(Kij * tau)
%  Nij (MxM matrix) - transition conditional count matrix
%    Nij(i,j) is number of observed counts terminating in j at time tau given that they initially started in i at time zero.
%  tau - evolution time to relate rate matrix to transition matrix through Tij = expm(Kij * tau)
%
% OPTIONAL ARGUMENTS
%  mode - if mode is 'sample', sampling will take place; if mode is 'optimize', then MC optimization is used
%
% RETURN VALUES
%  Tij_new (MxM matrix) - updated sample of transition matrix, may be correlated with previous sample
%
% NOTES
%  The transition matrix is guaranteed to satisfy detailed balance and correspond to a rate matrix with positive off-diagonal entries.
%
%  The algorithm described in [1] and modified in [2] (especially with regard to the prior) is implemented.
%  Transitions between all states are allowed -- no constraints on the sparisty of Kij are imposed.
%
% REFERENCES
%  [1] Sriraman S, Kevrekidis IG, and Hummer G. Coarse master equation from Bayesian analysis of replica molecular dynamics simulations. JPC B 109(14):6479-6484, 2005.
%  [2] Buchete N-V and Hummer G. Coarse master equations for peptide folding dynamics. JPC B 112(19):6057-6069, 2008.
%
% TODO
%  * Check to be sure we don't need a correction for asymmetric proposal probability.
%  * Automatically adjust delta to provide ~ 50% acceptance before fixing delta.
%  * Automatically determine 'nsteps' to provide an uncorrelated sample.
%  * Speed up sampling with Java methods?
%  * Consider starting sampling from maximum-likelihood estimate if we're going to produce a completely decorrelated estimate.

% Handle optional arguments
if (nargin == 2)
  mode = 'sample';
end

% PARAMETERS
nsteps = 1000; % number of update steps (set to give approximately uncorrelated new sample)
%delta = 0.1; % move width in logarithmic parameterization (set to give acceptance probability of ~ 50%)
delta = 0.05; % move width in logarithmic parameterization (set to give acceptance probability of ~ 50%)

% Create copies of things.
logKij = model.logKij;
logPi = model.logPi;
Kij = model.Kij;
Tij = model.Tij;

% Determine size of transition matrix.
M = size(Tij,1);

% Promote Nij to double.
Nij = double(Nij);

% Compute compressed vector
x = pack_vector2(logKij, logPi);

% Compute initial log-likelihood.
logL = log_likelihood(model.Tij, Nij);

% Storage for determining correlation time.
x_ni = zeros(nsteps,M-1 + M*(M-1)/2);

% Update sample using Hummer scheme.
naccept = 0; % number of accepted moves
for step = 1:nsteps

  % Choose an element of x with uniform probability.
  i = floor(rand * length(x)) + 1;
  
  % Perturb this element.
  x_new = x;
  x_new(i) = x_new(i) + (2*rand - 1) * delta;
  
  % Reject if any state probabilities exceed 1.
  if (sum(exp(x_new(1:(M-1)))) < 1)
    % Compute new transition matrix.
    [logKij_new, logPi_new, Kij_new, Tij_new] = unpack_vector2(x_new, model.tau);

    % Ensure all eigenvalues are positive
    smallest_eigenvalue = eigs(Tij_new, 1, 'SR', struct('disp', 0));   
    if (smallest_eigenvalue <= eps * M)
      %smallest_eigenvalue
      % reject this transition matrix if any eigenvalues are nonpositive
      continue;
    end
    
    % Check Tij.
    if (any(~isreal(Tij_new)))
      step
      i
      delta
      [x x_new x_new-x]
      Tij
      Tij_new
      error('proposed Tij is not strictly real');
    end      
  
%    % Ensure all states are metastable.
%    if any(diag(Tij_new) < 0.5)
%      %diag(Tij_new)
%      continue;
%    end
    
    % Compute new log-likelihood function.
    logL_new = log_likelihood(Tij_new, Nij);
    
    % Compute change.
    delta_logL = logL_new - logL;
    
    % DEBUG
    %disp(sprintf('step %6d : logL_new = %16.1f : logL = %16.1f : delta_logL = %16.1f', step, logL_new, logL, delta_logL));
    
    % Accept or reject by Metropolis-Hastings procedure.
    accept = false;
    if (delta_logL >= 0.0)
      accept = true;
    elseif (strcmp(mode, 'sample') && (rand < exp(delta_logL)))
      accept = true;
    end
    if (accept)
      % Store new parameters.
      x = x_new;
      logKij = logKij_new;
      logPi = logPi_new;
      Kij = Kij_new;
      Tij = Tij_new;      
      logL = logL_new;
      naccept = naccept + 1;
    end
  end
  
  % store
  x_ni(step,:) = x;
end

% Show summary statistics.
MIN_ACCEPT = 0.2; % minimum acceptance fraction before warning is issued
if (naccept/nsteps < MIN_ACCEPT)
  disp(sprintf('\nWARNING: accepted %d / %d move attempts (%.1f %%)', naccept, nsteps, naccept / nsteps * 100));
end

% DEBUG make sure x is strictly real
if ~isreal(x)
  x
  smallest_eigenvalue = eigs(Tij, 1, 'SR', struct('disp', 0)); 
  smallest_eigenvalue
  error('final x is not strictly real');
end

% Check Tij.
if (any(~isreal(Tij)))
  Tij
  smallest_eigenvalue = eigs(Tij, 1, 'SR', struct('disp', 0)); 
  smallest_eigenvalue
  error('final Tij is not strictly real');
end    

% Store.
model.logKij = logKij;
model.logPi = logPi;
model.Kij = Kij;
model.Tij = Tij;

return

%%%%%%%%%%%%%%

function logL = log_likelihood(Tij, Nij)
% Compute log-likelihood
%
% ARGUMENTS
%  Tij (MxM) - row-stochastic transition matrix - Tij(i,j) is probability of observing system in state j given that it was initial in state i.
%   We require all Tij(i,j) to be strictly positive.
%  Nij (MxM) - row transition count matrix - Nij(i,j) is number of transitions observed terminating in j given that the system was initially in state i
%
% RETURN VALUE
%  logL - log-likelihood

% Compute elementwise logarithm of Tij.
logTij = log(Tij); 

% Compute log-likelihood
logL = sum(sum(Nij .* logTij));

return
