function analyze_transition_matrix(models, options)
% Analyze transition probability or rate constant distribution.
%
% TODO: Allow options.nbins to override hard-coded nbins default below.

nstates = models(1).nstates;
nmodels = length(models);
nbins = 40; % number of histogram bins for plotting
ci = 0.95; % confidence interval to report

Tijn = zeros(nstates, nstates, nmodels);
for n = 1:nmodels
  Tijn(:,:,n) = models(n).Tij;
end

% Plot.
clf;
for i = 1:nstates
  for j = 1:nstates
    subplot(nstates,nstates,nstates*(i-1)+j);
    Tn = squeeze(Tijn(i,j,:));
    Tmean = mean(Tn);
    [Tlow, Thigh] = confidence_interval(Tn, 0.5 - ci/2, 0.5 + ci/2);
    [N,X] = hist(Tn, nbins);
    bar(X,N,'grouped','k');

    % Label plot.
    title(sprintf('%.5f [%.5f, %.5f]', Tmean, Tlow, Thigh), 'FontSize', 7);
    xlabel('');
    ylabel(sprintf('T_{%d %d}', i, j));
    set(gca, 'YTick', []);    

    % Adjust axes.
    axis([-0.05 1.05 0 1.15*max(N)]);
  end
end

return
