function analyze_rates(models, options)
% Analyze transition probability or rate constant distribution.

nstates = models(1).nstates;
nmodels = length(models);
nbins = 30; % number of histogram bins for plotting
ci = 0.95; % confidence interval

Tijn = zeros(nstates, nstates, nmodels);
for n = 1:nmodels
  Tijn(:,:,n) = models(n).Tij;
  Kijn(:,:,n) = logm(models(n).Tij) / options.tau;
end

tau_max = 0.0;
k_max = 0.0;

% Determine maxes.
for i = 1:nstates
  for j = 1:nstates
    subplot(nstates,nstates,nstates*(i-1)+j);
    
    if i == j
      % show lifetimes
      tau_n = options.tau ./ (1.0 - squeeze(Tijn(i,j,:)));
      tau_max = max(max(tau_n), tau_max);
    else
      % show rates
      Kn = squeeze(Kijn(i,j,:));
      k_max = max(max(Kn), k_max);
    end
  end
end


% Plot.
clf;
for i = 1:nstates
  for j = 1:nstates
    subplot(nstates,nstates,nstates*(i-1)+j);
    
    if i == j
      % show lifetimes
      tau_n = options.tau ./ (1.0 - squeeze(Tijn(i,j,:)));
      tau_mean = mean(tau_n);
      [tau_low, tau_high] = confidence_interval(tau_n, 0.5-ci/2, 0.5+ci/2);
      [N,X] = hist(tau_n, nbins);
      bar(X,N,'grouped','r', 'edgecolor', 'r');
      title(sprintf('%.4f [%.4f, %.4f]', tau_mean, tau_low, tau_high), 'FontSize', 7); % DEBUG
      if (isfield(options, 'time_units'))
	xlabel(sprintf('%s', options.time_units));
      end
      ylabel(sprintf('\\tau_{%d}', i));      
      axis([min(tau_n) max(tau_n) 0 1.15*max(N)]);
      axis([0 tau_max 0 1.15*max(N)]); % override
      set(gca, 'YTick', []);    
    else
      % show rates
      Kn = squeeze(Kijn(i,j,:));
      Kmean = mean(Kn);
      [Klow, Khigh] = confidence_interval(Kn, 0.5-ci/2, 0.5+ci/2);
      [N,X] = hist(Kn, nbins);
      bar(X,N,'grouped','k', 'edgecolor', 'k');
      title(sprintf('%.4f [%.4f, %.4f]', Kmean, Klow, Khigh), 'FontSize', 7); % DEBUG
      if (isfield(options, 'time_units'))
	xlabel(sprintf('%s^{-1}', options.time_units));
      end
      ylabel(sprintf('K_{%d %d}', i, j));      
      axis([min(Kn) max(Kn) 0 1.15*max(N)]);
      axis([0 k_max 0 1.15*max(N)]); % override
      set(gca, 'YTick', []);    
    end
  end
end

return
