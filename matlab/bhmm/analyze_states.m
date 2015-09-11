function estimates = analyze_states(models, options)
% Analyze state parameter distributions from sample of models.
%
%  analyze_states(models, options)
%
% ARGUMENTS
%   models - set of models
%   options - options structure
%     options.sameaxis - if true, will put all state plots on same scale; if false, each plot will have different scale
%
% RETURNS
%   estimates - struct with means and std devs
%               means: mu (average observable), sigma (std dev of observable), tau (lifetimes)
%               std devs: dmu, dsigma, dtau

% PARAMETERS
nbins = 20; % number of histogram bins

% Determine number of states.
nstates = models(1).nstates;

% Determine number of models.
nmodels = length(models);

% Extract state parameters into arrays for analysis.
mu = zeros(nstates, nmodels);
sigma = zeros(nstates, nmodels);
lifetimes = zeros(nstates, nmodels);
for n = 1:nmodels
  model = models(n);
  for state = 1:nstates
    mu(state,n) = model.states{state}.mu;
    sigma(state,n) = model.states{state}.sigma;
    lifetimes(state,n) = options.tau * 1.0 / (1.0 - model.Tij(state,state));
  end
end

% Compute means and std devs of various properties.
estimates = struct();
estimates.mu = mean(mu');
estimates.dmu = std(mu');
estimates.sigma = mean(sigma');
estimates.dsigma = std(sigma');
estimates.tau = mean(lifetimes');
estimates.dtau = std(lifetimes');

% Plot distributions for each state.
clf;
ny = nstates;
nx = 3;

mu_min = min(min(mu));
mu_max = max(max(mu));

sigma_min = min(min(sigma));
sigma_max = max(max(sigma));

tau_min = min(min(lifetimes));
tau_max = max(max(lifetimes));

for state = 1:nstates
  % plot mu
  subplot(ny,nx,3*(state-1)+1)
  [N,X] = hist(mu(state,:), nbins);
  h = bar(X,N,'grouped','r');
  set(h, 'edgecolor', 'none');
  if (isfield(options, 'observable_units'))
    xlabel(sprintf('\\mu / %s', options.observable_units));
    title(sprintf('\\mu = %.2f \\pm %.2f %s', estimates.mu(state), estimates.dmu(state), options.observable_units));
  else
    xlabel('\mu');
    title(sprintf('\\mu = %.2f \\pm %.2f', estimates.mu(state), estimates.dmu(state)));
  end
  ylabel('P(\mu)');
  set(gca, 'YTick', []);
  oldaxis = axis;
  if options.sameaxis
    axis([mu_min mu_max oldaxis(3:4)]);
  end

  subplot(ny,nx,3*(state-1)+2)
  [N,X] = hist(sigma(state,:), nbins);
  h = bar(X,N,'grouped','k');
  set(h, 'edgecolor', 'none');
  if (isfield(options, 'observable_units'))
    xlabel(sprintf('\\sigma / %s', options.observable_units));
    title(sprintf('\\sigma = %.2f \\pm %.2f %s', estimates.sigma(state), estimates.dsigma(state), options.observable_units));
  else
    xlabel('\sigma');
    title(sprintf('\\sigma = %.2f \\pm %.2f', estimates.sigma(state), estimates.dsigma(state)));
  end
  ylabel('P(\sigma)');
  set(gca, 'YTick', []);    
  oldaxis = axis;
  if options.sameaxis
    axis([sigma_min sigma_max oldaxis(3:4)]);
  end

  subplot(ny,nx,3*(state-1)+3)
  [N,X] = hist(lifetimes(state,:), nbins);
  h = bar(X,N,'grouped','b');
  set(h, 'edgecolor', 'none');
  if (isfield(options, 'observable_units'))
    xlabel(sprintf('\\tau / %s', options.time_units));
    title(sprintf('\\tau = %.2f \\pm %.2f %s', estimates.tau(state), estimates.dtau(state), options.time_units));
  else
    xlabel('\tau');
    title(sprintf('\\tau = %.2f \\pm %.2f', estimates.tau(state), estimates.dtau(state)));
  end
  ylabel('P(\tau)');
  set(gca, 'YTick', []);    
  oldaxis = axis;
  if options.sameaxis
    axis([tau_min tau_max oldaxis(3:4)]);
  end
end

return

