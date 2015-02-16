% Plot confidence interval data.

clf;

% Pi
subplot(2,2,1);
plot([0 1], [0 1], 'k-');
hold on;
[Plow, Phigh] = beta_confidence_intervals(ci_Pi, ntrials, 0.95);
%errorbar(repmat(cis',1,nstates), ci_Pi, sqrt(ci_Pi .* (1-ci_Pi))/sqrt(ntrials),'.');
errorbar(repmat(cis',1,nstates), ci_Pi, ci_Pi - Plow, Phigh - ci_Pi, '.');
axis square;
set(gca, 'XTick', [0 1]);
set(gca, 'YTick', [0 1]);
%title('\pi_i');
text(0.1, 0.85, '\pi_i');

% Tij

ci_X = zeros(ncis, nstates*nstates);
for i = 1:nstates
  for j = 1:nstates
    ci_X(:,i+nstates*(j-1)) = ci_Tij(:,i,j);
  end
end

% TODO: Compute more accurate confidence intervals for these estimates using Bayesian inference.  See earlier work.

subplot(2,2,2);
plot([0 1], [0 1], 'k-');
hold on;
[Plow, Phigh] = beta_confidence_intervals(ci_X, ntrials, 0.95);
%errorbar(repmat(cis',1,nstates*nstates), ci_X, sqrt(ci_X .* (1-ci_X))/sqrt(ntrials),'.');
errorbar(repmat(cis',1,nstates*nstates), ci_X, ci_X - Plow, Phigh - ci_X,'.');
axis square;
set(gca, 'XTick', [0 1]);
set(gca, 'YTick', [0 1]);
%title('T_{ij}');
text(0.1, 0.85, 'T_{ij}');

% mu
subplot(2,2,3);
plot([0 1], [0 1], 'k-');
hold on;
[Plow, Phigh] = beta_confidence_intervals(ci_mu_i, ntrials, 0.95);
%errorbar(repmat(cis',1,nstates), ci_mu_i, sqrt(ci_mu_i .* (1-ci_mu_i))/sqrt(ntrials),'.');
errorbar(repmat(cis',1,nstates), ci_mu_i, ci_mu_i - Plow, Phigh - ci_mu_i, '.');
axis square;
set(gca, 'XTick', [0 1]);
set(gca, 'YTick', [0 1]);
%title('\mu_i');
text(0.1, 0.85, '\mu_i');
xlabel('expected');
ylabel('observed');

% sigma
subplot(2,2,4);
plot([0 1], [0 1], 'k-');
hold on;
[Plow, Phigh] = beta_confidence_intervals(ci_sigma_i, ntrials, 0.95);
%errorbar(repmat(cis',1,nstates), ci_sigma_i, sqrt(ci_sigma_i .* (1-ci_sigma_i))/sqrt(ntrials),'.');
errorbar(repmat(cis',1,nstates), ci_sigma_i, ci_sigma_i - Plow, Phigh - ci_sigma_i, '.');
axis square;
set(gca, 'XTick', [0 1]);
set(gca, 'YTick', [0 1]);
%title('\sigma_i');
text(0.1, 0.85, '\sigma_i');

% Export figure.
exportfig(gcf, 'confidence-intervals.eps', 'color', 'cmyk', 'width', 3.5, 'height', 3.5);
system('epstopdf confidence-intervals.eps');