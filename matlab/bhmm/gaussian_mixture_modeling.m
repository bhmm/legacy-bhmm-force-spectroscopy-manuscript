% Gaussian mixture modeling with variable number of normal components
%
% On Bayesian Analysis of Mixtures with an Unknown Number of Components 
% Sylvia Richardson; Peter J. Green 
% Journal of the Royal Statistical Society. Series B (Methodological), Vol. 59, No. 4. (1997), pp. 731-792. 
% http://www.jstor.org/pss/2985194


clear;

% TEST PROBLEM
mu_k = [0.0 2.0 4.0]; % means
sigma_k = [0.5 0.5 0.5]; % standard deviations
N_k = [500 500 500]; % number of samples
x_n = zeros(sum(N_k), 1); % x_n(n) is sample n
n = 0;
for k = 1:length(mu_k)  
  x_n(n+1:n+N_k(k)) = sigma_k(k) * randn(N_k(k),1) + mu_k(k);
  n = n + N_k(k);
end
clear mu_k sigma_k N_k n k;

% PARAMETERS
niterations = 1000; % number of sampler sweeps to run
K_max = 20; % maximum number of normal components

% PRIOR 
delta = 1.0; % pseudocounts for mixture proportion
R = max(x_n) - min(x_n); % range of data
xi = min(x_n) + R/2; % mean of prior on normal means
kappa = 1 / R^2; 
alpha = 2;
g = 0.2; % gamma parameter on beta
h = 10 / R^2; % gamma parameter on beta
p_k = ones(1,K_max) / K_max; % p_k(k) is the prior on k states
b_k = ones(1,K_max) * 0.5; b_k(1) = 1; b_k(K_max) = 0; % b_k(k) is the probability of a birth process when there are K components
d_k = ones(1,K_max) * 0.5; d_k(1) = 0; d_k(K_max) = 1; % d_k(k) is the probability of a death process when there are K components

% Count number of samples.
N = length(x_n);

% INITIAL SAMPLER CONDITIONS
K = 1; % number of normal component
w_k = 1; % w_k(k) is the fraction or proportion of normal component k
mu_k = mean(x_n); % mu_k(k) is the mean of normal component k
sigma_k = std(x_n); % sigma_k(k) is the std dev of normal component k
z_n = ones(N, 1); % assign all samples to state 1
betaprm = 1.0;
N_k = N;

% HISTORY
history = cell(niterations,1);

% GIBBS SAMPLING
for iteration = 1:niterations
  % DEBUG
  fprintf('iteration %d\n', iteration);
  disp(w_k);
  disp(mu_k);
  disp(sigma_k);
  disp(N_k);
  
  % Update weights
  % w ~ Dirichlet(N_k + delta)
  % We use a gamma random number generator because a Dirichlet sampler is not available.
  N_k = hist(z_n, 1:K); % N_k(k) is the number of samples assigned to state k
  w_k = gamrnd(N_k + delta, 1); w_k = w_k / sum(w_k); % w_k(k) is the fraction of normal component k
  
  % Update normal parameters.
  for k = 1:K
    indices = find(z_n == k); % find samples assigned to state k    
    mu = (sigma_k(k)^(-2) * sum(x_n(indices)) + kappa*xi) / (sigma_k(k)^(-2) * N_k(k) + kappa);
    sigma = (sigma_k(k)^(-2) * N_k(k) + kappa)^(-1/2);
    mu_k(k) = sigma * randn() + mu;    
    sigma_k(k) = gamrnd(alpha + 0.5*N_k(k), (betaprm + 0.5*sum( (x_n(indices) - mu_k(k)).^2 ))^(-1))^(-1/2);
  end
  
  % Update the assignment of samples to states.
  for n = 1:N    
    log_P_k = log(w_k ./ sigma_k) -0.5 * ((x_n(n) - mu_k)./sigma_k).^2; 
    P_k = exp(log_P_k - logsum(log_P_k));
    z_n(n) = draw(P_k);
  end

  % Update the scale hyperparameter beta.
  betaprm = gamrnd(g + K * alpha, (h + sum(sigma_k.^(-2)))^(-1));

  % SANITY CHECK
  if(any(z_n < 0) || any(z_n > K))
    error('error with indices z_n');
  end
  
  % Attempt to split or combine normal components.
  if (rand < b_k(K)) 
    % Attempt to split a component.
    % Choose component to split.
    k = unidrnd(K);
    fprintf('Attempting to split component %d\n', k);
    % Copy original state parameters.
    w_j = w_k(k);
    mu_j = mu_k(k);
    sigma2_j = sigma_k(k)^2;    
    % Choose degrees of freedom in split.
    u1 = betarnd(2,2);
    u2 = betarnd(2,2);
    u3 = betarnd(1,1);
    % Propose split state parameters;
    w_j1 = w_j * u1;
    w_j2 = w_j * (1-u1);
    mu_j1 = mu_j - u2*sqrt(sigma2_j) * sqrt(w_j2 / w_j1);
    mu_j2 = mu_j + u2*sqrt(sigma2_j) * sqrt(w_j1 / w_j2);
    sigma2_j1 =    u3 *(1-u2^2)*sigma2_j * w_j/w_j1;
    sigma2_j2 = (1-u3)*(1-u2^2)*sigma2_j * w_j/w_j2;
    % Reject if two new normal distributions are not adjacent.    
    [mu_k_sorted, indices] = sort([mu_j1 mu_j2 mu_k(1:k-1) mu_k(k+1:end)]);
    new_states_are_adjacent = any((indices(1:end-1)==1 & indices(2:end)==2)) || (any(indices(1:end-1)==2 & indices(2:end)==1));
    if(new_states_are_adjacent)    
      % Update the assignment of samples to states.
      l1 = 0; % number of samples assigned to j1
      l2 = 0; % number of samples assigned to j2
      log_P_alloc = 0.0; % log proabability of this particular allocation of samples to j1 and j2
      log_likelihood_ratio = 0.0; % log likelihood ratio for observing data in j1 and j2 rather than j
      indices = find(z_n == k)'; % indices of those samples in original state j
      znew_n = z_n; % new allocation
      for n = indices
	log_P_1 = log(w_j1 ./ sqrt(sigma2_j1)) - 0.5 * ((x_n(n) - mu_j1)^2/sigma2_j1);
	log_P_2 = log(w_j2 ./ sqrt(sigma2_j2)) - 0.5 * ((x_n(n) - mu_j2)^2/sigma2_j2);
	log_Z = logsum([log_P_1 log_P_2]);
	log_P_1 = log_P_1 - log_Z;
	log_P_2 = log_P_2 - log_Z;
	
	if (rand <= exp(log_P_1))
	  l1 = l1 + 1; % increment count of samples allocated to j1
	  znew_n(n) = k; % j1 will eventually replace index k
	  log_P_alloc = log_P_alloc + log_P_1;
	  log_likelihood_ratio = log_likelihood_ratio - 0.5*log(sigma2_j1) - 0.5*((x_n(n) - mu_j1)^2/sigma2_j1) + 0.5*log(sigma2_j) + 0.5*((x_n(n) - mu_j)^2/sigma2_j);	  
	else
	  l2 = l2 + 1; % increment count of samples allocated to j2
	  znew_n(n) = K+1; % j2 will eventually replace index K+1
	  log_P_alloc = log_P_alloc + log_P_2;
	  log_likelihood_ratio = log_likelihood_ratio - 0.5*log(sigma2_j2) - 0.5*((x_n(n) - mu_j2)^2/sigma2_j2) + 0.5*log(sigma2_j) + 0.5*((x_n(n) - mu_j)^2/sigma2_j);	  
	end
      end
      
      log_P_alloc
      log_likelihood_ratio
      
      % Accept or reject
      log_A = log_likelihood_ratio + log(K+1) + log(p_k(K+1)) - log(p_k(K)) + log(K+1) + (delta-1+l1)*log(w_j1) + (delta-1+l2)*log(w_j2) - (delta-1+l1+l2)*log(w_j) - betaln(delta,K*delta) + (1/2)*log(kappa/(2*pi)) - (1/2)*kappa*((mu_j1-xi)^2+(mu_j2-xi)^2-(mu_j-xi)^2) + alpha*log(betaprm) - gammaln(alpha) - (alpha+1)*log(sigma2_j1*sigma2_j2/sigma2_j) - betaprm*(1/sigma2_j1 + 1/sigma2_j2 - 1/sigma2_j) + log(d_k(K+1)) - log(b_k(K)) - log_P_alloc - log(betapdf(u1,2,2)*betapdf(u2,2,2)*betapdf(u3,1,1)) + log(w_j) + log(abs(mu_j1-mu_j2)) + log(sigma2_j1*sigma2_j2/sigma2_j) - log(u2) - log(1-u2^2) - log(u3) - log(1-u3);
      log_A
      if(isnan(log_A))
	error('log_A is NaN');
      end
      
      if((log_A > 0) || (rand < exp(log_A)))
	% accept
	disp('Split accepted');
	K = K + 1;
	w_k(k) = w_j1;
	w_k(K) = w_j2;
	mu_k(k) = mu_j1;
	mu_k(K) = mu_j2;
	sigma_k(k) = sqrt(sigma2_j1);
	sigma_k(K) = sqrt(sigma2_j2);
	z_n = znew_n;
      else
	% reject
	disp('Split rejected.');
      end
    else
      disp('Split is not adjacent.')
    end
  else
    % Attempt to combine two components.
    % Sort means
    [mu_k_sorted, indices] = sort(mu_k);
    i = unidrnd(K-1);
    j1 = indices(i);
    j2 = indices(i+1);
    % Duplicate
    w_j1 = w_k(j1); mu_j1 = mu_k(j1); sigma2_j1 = sigma_k(j1)^2;
    w_j2 = w_k(j2); mu_j2 = mu_k(j2); sigma2_j2 = sigma_k(j2)^2;
    % Combine        
    w_j = w_j1 + w_j2;
    mu_j = (w_j1*mu_j1 + w_j2*mu_j2) / w_j;
    sigma2_j = (w_j1*(mu_j1^2 + sigma2_j1) + w_j2*(mu_j2^2 + sigma2_j2))/w_j - mu_j^2;
    % Compute allocation probability.
    l1 = 0;
    l2 = 0;
    log_P_alloc = 0.0;
    log_likelihood_ratio = 0.0;
    indices = find((z_n == j1) | (z_n == j2))';
    for n = indices
      log_P_1 = log(w_j1 ./ sqrt(sigma2_j1)) - 0.5 * ((x_n(n) - mu_j1)^2/sigma2_j1);
      log_P_2 = log(w_j2 ./ sqrt(sigma2_j2)) - 0.5 * ((x_n(n) - mu_j2)^2/sigma2_j2);
      log_Z = logsum([log_P_1 log_P_2]);
      log_P_1 = log_P_1 - log_Z;
      log_P_2 = log_P_2 - log_Z;
	
      if (z_n(n) == j1)
	l1 = l1 + 1;
	log_P_alloc = log_P_alloc + log_P_1;
	log_likelihood_ratio = log_likelihood_ratio - 0.5*log(sigma2_j1) - 0.5*((x_n(n) - mu_j1)^2/sigma2_j1) + 0.5*log(sigma2_j) + 0.5*((x_n(n) - mu_j)^2/sigma2_j);	  	  
      else
	l2 = l2 + 1;
	log_P_alloc = log_P_alloc + log_P_2;	
	log_likelihood_ratio = log_likelihood_ratio - 0.5*log(sigma2_j2) - 0.5*((x_n(n) - mu_j2)^2/sigma2_j2) + 0.5*log(sigma2_j) + 0.5*((x_n(n) - mu_j)^2/sigma2_j);	  
      end
    end
    
    % Accept or reject
    K = K - 1;
    log_A = log_likelihood_ratio + log(K+1) + log(p_k(K+1)) - log(p_k(K)) + log(K+1) + (delta-1+l1)*log(w_j1) + (delta-1+l2)*log(w_j2) - (delta-1+l1+l2)*log(w_j) - betaln(delta,K*delta) + (1/2)*log(kappa/(2*pi)) - (1/2)*kappa*((mu_j1-xi)^2+(mu_j2-xi)^2-(mu_j-xi)^2) + alpha*log(betaprm) - gammaln(alpha) - (alpha+1)*log(sigma2_j1*sigma2_j2/sigma2_j) - betaprm*(1/sigma2_j1 + 1/sigma2_j2 - 1/sigma2_j) + log(d_k(K+1)) - log(b_k(K)) - log_P_alloc - log(betapdf(u1,2,2)*betapdf(u2,2,2)*betapdf(u3,1,1)) + log(w_j) + log(abs(mu_j1-mu_j2)) + log(sigma2_j1*sigma2_j2/sigma2_j) - log(u2) - log(1-u2^2) - log(u3) - log(1-u3);
    K = K + 1;
    log_A
    if((-log_A > 0) || (rand < exp(-log_A)))
      % Accept combine move.
      % Swap higher index state with state K so we can discard it.
      j = max(j1,j2);      
      [w_k(j), w_k(K)] = deal(w_k(K), w_k(j));      
      [mu_k(j), mu_k(K)] = deal(mu_k(K), mu_k(j));
      [sigma_k(j), sigma_k(K)] = deal(sigma_k(K), sigma_k(j));
      indices = find(z_n == K);
      z_n(indices) = j;
      % Discard state K, shrinking parameter vector.
      K = K - 1;
      w_k = w_k(1:K);
      mu_k = mu_k(1:K);
      sigma_k = sigma_k(1:K);
      % Replace lower index state with new combined state.
      j = min(j1,j2);            
      w_k(j) = w_j;
      mu_k(j) = mu_j;
      sigma_k(j) = sqrt(sigma2_j);
      % Reallocate its members.
      indices = find((z_n == j1) | (z_n == j2))';
      z_n(indices) = j;
      disp('Merge accepted.');
    else
    end
    
  end
  
  % SANITY CHECK
  if(any(z_n < 0) || any(z_n > K))
    error('error with indices z_n');
  end
  
  % Handle birth/death process.  
  N_k = hist(z_n, 1:K); % N_k(k) is the number of samples assigned to state k
  indices = find(N_k == 0); % indices of empty components
  k0 = length(indices); % number of empty components (components with no samples assigned to them)
  if (rand < b_k(K)) 
    % Attempt birth process.
    k = K+1;
    w_k(K+1) = betarnd(1,K)
    mu_k(K+1) = randn*kappa^(-1/2) + xi
    sigma_k(K+1) = gamrnd(alpha, 1/betaprm)^(-1/2)
    % Compute acceptance probability.
    log_A = log(p_k(K+1)) - log(p_k(K)) - betaln(K*delta,delta) + (delta-1)*log(w_k(k)) + (N+K*delta-K)*log(1-w_k(k)) + log(K+1) + log(d_k(K+1)) - log(k0+1) - log(b_k(K)) - log(betapdf(w_k(k),1,K)) + K*log(1-w_k(k))
    if((log_A > 0) || (rand < exp(log_A)))
      % Accept.
      w_k(1:K) = w_k(1:K) * (1 - w_k(k));
      K = K + 1;
      disp('Birth process accepted.');
    else
      % Reject.
      disp('Birth process rejected.');
      % expunge trial state
      w_k = w_k(1:K);
      mu_k = mu_k(1:K);
      sigma_k = sigma_k(1:K);
    end   
  else
    % Attempt death process if there are any empty components.
    if (k0 > 0)
      % Select an empty component at random.
      k = indices(unidrnd(k0)); % state to remove
      % Compute acceptance probability.
      K = K - 1;
      log_A = log(p_k(K+1)) - log(p_k(K)) - betaln(K*delta,delta) + (delta-1)*log(w_k(k)) + (N+K*delta-K)*log(1-w_k(k)) + log(K+1) + log(d_k(K+1)) - log(k0+1) - log(b_k(K)) - log(betapdf(w_k(k),1,K)) + K*log(1-w_k(k));
      K = K + 1;
      if((-log_A > 0) || (rand < exp(-log_A)))
	% Accept death process.
	% Swap empty state k with state K so we can contract parameter vector.
	[w_k(k), w_k(K)] = deal(w_k(K), w_k(k));      
	[mu_k(k), mu_k(K)] = deal(mu_k(K), mu_k(k));
	[sigma_k(k), sigma_k(K)] = deal(sigma_k(K), sigma_k(k));	
	indices = find(z_n == K); % switch indices
	z_n(indices) = k;
	% Delete state K.
	K = K - 1;
	w_k = w_k(1:K);
	mu_k = mu_k(1:K);
	sigma_k = sigma_k(1:K);	
	% Renormalize
	w_k = w_k / sum(w_k);
	disp('Death process accepted.');	
      else
	% Reject.
	% Restore number of states.
	disp('Death process rejected.');
      end      
    end    
  end  

  % SANITY CHECK
  if(any(z_n < 0) || any(z_n > K))
    error('error with indices z_n');
  end

  % Store history.
  history{iteration} = struct('K', K, 'w_k', w_k, 'mu_k', mu_k, 'sigma_k', sigma_k);
  K_t(iteration) = K;
end
