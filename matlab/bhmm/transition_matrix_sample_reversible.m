function Tij = transition_matrix_sample_reversible(Tij, Nij, options)
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
%    if options.prior_transition_pseudocounts exists (a KxK matrix where the (i,j) element represents transitions from i to j), these pseudocounts are included; all ones is uniform prior
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
%
% REFERENCES
%  [1] Noe F. Probability distributions of molecular observables computed from Markov models. JCP 128:244103, 2008.
%  [2] Chodera JD and Noe F. Probability distributions of molecular observables computed from Markov models. 
%      II. Uncertainties in observables and their time-evolution. Submitted, 2009.
%
% TODO
%  * Automatically adjust delta to provide ~ 50% acceptance before fixing delta.
%  * Automatically determine 'nsteps' to provide an uncorrelated sample.
%  * Speed up sampling with Java methods?
%  * Consider starting sampling from maximum-likelihood estimate if we're going to produce a completely decorrelated estimate.

% PARAMETERS
nsteps = 10000; % number of update steps (set to give approximately uncorrelated new sample)
nupdate = 100; % number of accepted steps between updating pi_i.
ROW_SUM_TOLERANCE = 1.0e-3; % row sum tolerance

if (options.reversible == false)
  error('Non-reversible sampling not yet implemented.');
end

% Determine size of transition matrix.
M = size(Tij,1);

% Include prior transition pseudocounts.
if isfield(options, 'prior_transition_pseudocounts')
  Nij = double(Nij) + double(options.prior_transition_pseudocounts-1);
end

% Promote Nij to double.
Nij = double(Nij);

% Compute total number of observed transitions out of each state.
Ni = sum(Nij,2);

% Compute stationary probability.
Pi = stationary_probability(Tij);

% Keep track of how many steps have been accepted.
naccepted = 0;

%% Sanity checks.

% Check initial transition matrix does not contain negative values.
if (any(any(Tij < 0.0)))
  Tij
  error('Initial Tij has negative elements.')
end

% Check initial row sums.
row_sums = sum(Tij,2);
if any(abs(row_sums - 1) > ROW_SUM_TOLERANCE)
  Tij
  error('Initial Tij has some row sums that differ from unity.')
end

% Perform a number of sampling steps.
for step = 1:nsteps
  initial_Tij = Tij;

  %% Generate proposal move.

  if (rand() < 0.5)
    % Reversible element shift (Section II.B of [1], Section III.A of [2])
    move_type = 'reversible element shift';

    % Choose transition element (i,j) where i \ne j, to perturb.
    i = randi(M);
    j = randi(M-1); 
    if (j >= i)
      j = j + 1;
    end

    % Choose element perturbation uniformly over allowed range.
    min_delta = max(-Tij(i,i), -Pi(j)/Pi(i)*Tij(j,j)); 
    max_delta = Tij(i,j);
    delta = rand() * (max_delta - min_delta) + min_delta;

    % Generate perturbed transition matrix.
    Tij_new = Tij;
    Tij_new(i,j) = Tij(i,j) - delta;
    Tij_new(i,i) = Tij(i,i) + delta;
    Tij_new(j,i) = Tij(j,i) - Pi(i) / Pi(j) * delta;
    Tij_new(j,j) = Tij(j,j) + Pi(i) / Pi(j) * delta;

    % Stationary probabilities do not need updating.
    Pi_new = Pi;

    if (options.verbosity >= 6)
      disp(sprintf('  T(%d,%d) = %8.6f -> %8.6f', i, j, Tij(i,j), Tij_new(i,j)));
    end

    % Compute log proposal probability ratio (Eq. 22 of [2])
    log_proposal_ratio = 0.5*log( (Tij_new(i,j)^2 + Tij_new(j,i)^2) / (Tij(i,j)^2 + Tij(j,i)^2) );
    
    % Compute log posterior probability ratio (Eq. 23 of [2])
    log_posterior_ratio = Nij(i,i)*log(Tij_new(i,i) / Tij(i,i)) + Nij(i,j)*log(Tij_new(i,j)/Tij(i,j)) + Nij(j,j)*log(Tij_new(j,j)/Tij(j,j)) + Nij(j,i)*log(Tij_new(j,i)/Tij(j,i));

    if (options.verbosity >= 6)
      disp(sprintf('step %d / %d : element shift (%d,%d) | delta %f | proposal %f | posterior %f', step, nsteps, i, j, delta, log_proposal_ratio, log_posterior_ratio))
    end
  else
    % Row shift (Section II.C of Ref. [1], Section III.B of [2])
    move_type = 'row shift';

    % Choose row to shift.
    i = randi(M);
    j = -1; % dummy
    
    % Choose perturbation factor uniformly over allowed range.
    min_alpha = 0.0;
    max_alpha = 1.0/(1.0-Tij(i,i));
    alpha = rand() * (max_alpha - min_alpha) + min_alpha;

    % Generate perturbed transition matrix.
    Tij_new = Tij;
    Tij_new(i,:) = alpha * Tij(i,:);
    Tij_new(i,i) = alpha * (Tij(i,i) - 1.0) + 1.0;

    if (options.verbosity >= 6)
      disp('  old row:');
      disp(Tij(i,:));
      disp('  new row:');
      disp(Tij_new(i,:));
      disp(sprintf('  row sum: %f', sum(Tij_new(i,:),2)));
    end

    % Update stationary probabilities.
    Pi_new = Pi;
    factor = 1.0 / (Pi(i) + alpha * (1.0 - Pi(i)));
    Pi_new(:) = alpha * Pi(:) * factor;
    Pi_new(i) = Pi(i) * factor;
    %Pi_new = stationary_probability(Tij_new); % DEBUG: To avoid row normalization failure

    if (options.verbosity >= 6)
      disp('  old Pi:');
      disp(Pi');
      disp('  new Pi:');
      disp(Pi_new');
    end

    % Compute log proposal probability ratio (Eq. 24 of [2])
    log_proposal_ratio = double(M-2) * log(alpha);
    
    % Compute log posterior probability ratio (Eq. 25 of [2])
    log_posterior_ratio = (Ni(i) - Nij(i,i))*log(alpha) + Nij(i,i)*log(Tij_new(i,i)/Tij(i,i));

    if (options.verbosity >= 6)
      disp(sprintf('step %d / %d : row shift | row %d | alpha  %f | proposal %f | posterior %f', step, nsteps, i, alpha, log_proposal_ratio, log_posterior_ratio));
    end
  end

  if any(any(isnan(Tij_new)))
    if (options.verbosity >= 6)
      disp('rejected Tij_new for NaNs');
    end
    continue;
  end

  if any(isnan(Pi_new))
    if (options.verbosity >= 6)
      disp('rejected Pi_new for NaNs');
    end
    continue;
  end

  % Reject disallowed transition matrices.
  if (options.diagonally_dominant && any(diag(Tij_new) < 0.5))
    if (options.verbosity >= 6)
      disp('rejected disallowed transition matrix for not being diagonally dominant.');
    end
    continue;
  end

  % Reject transition matrices that have lost probability normalization.
  if (abs(sum(Pi_new) - 1.) > ROW_SUM_TOLERANCE)
%    step
%    move_type
%    [i,j]    
%    [min_alpha, max_alpha, alpha]
%    Tij
%    Tij_new
%    disp(sum(Pi_new));
%    error('Pi_new sum differs from unity.');
    continue;
  end

  % Reject transition matrices that have lost row normalization.
  if (any(abs(sum(Tij_new,2) - ones(M,1)) > ROW_SUM_TOLERANCE))
%    step
%    move_type
%    [i,j]    
%    [min_alpha, max_alpha, alpha]
%    Tij
%    Tij_new
%    disp(sum(Tij_new));
%    error('Proposed Tij has row sum that differs from unity.');
    continue;
  end

  % Reject transition matrices containing any negative values.
  % TODO: Figure out why this fails so often
  if (any(any(Tij_new < 0.0)))
%    step
%    move_type
%    [i,j]
%    if strcmp(move_type, 'row shift')
%      [min_alpha, max_alpha, alpha]
%    else
%      [min_delta, max_delta, delta]
%    end
%    Pi
%    Pi_new
%    Tij
%    Tij_new    
%    error('Proposed transition matrix contains negative values.');    
    continue;
  end

  log_P_accept = log_proposal_ratio + log_posterior_ratio;
  if ((log_P_accept >= 0.0) || (rand() < exp(log_P_accept)))
    % Accept.
    naccepted = naccepted + 1;
    Tij = Tij_new;
    Pi = Pi_new;
  else
    % Reject.
    continue;
  end

  % Update stationary probabilities and re-symmetrize matrix every so often.
  if (mod(naccepted, nupdate) == 0)
    Pi = stationary_probability(Tij);
    % TODO: Re-symmetrize matrix in case elements have drifted.    
  end
end

% Report on acceptance.
if (options.verbosity >= 5)
  disp(sprintf('Accepted %d / %d attempts (%.1f %%)', naccepted, nsteps, double(naccepted) / double(nsteps) * 100.0));
end

return
