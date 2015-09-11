% Plot confidence interval data for random model test.

clear;
load random_model_validation.mat;

clf;

% Pi
subplot(2,2,1);
plot([0 1], [0 1], 'k-');
hold on;
[Plow, Phigh] = beta_confidence_intervals(ci_Pi, ntrials, 0.95);
%errorbar(repmat(cis',1,nstates), ci_Pi, sqrt(ci_Pi .* (1-ci_Pi))/sqrt(ntrials),'.');
errorbar(cis', ci_Pi, ci_Pi - Plow, Phigh - ci_Pi, 'k.');
axis square;
set(gca, 'XTick', [0 1]);
set(gca, 'YTick', [0 1]);
%title('\pi_i');
text(0.1, 0.85, '\pi_i');

% Tij

% TODO: Compute more accurate confidence intervals for these estimates using Bayesian inference.  See earlier work.

subplot(2,2,2);
plot([0 1], [0 1], 'k-');
hold on;
[Plow, Phigh] = beta_confidence_intervals(ci_Tij, ntrials, 0.95);
%errorbar(repmat(cis',1,nstates*nstates), ci_Tij, sqrt(ci_Tij .* (1-ci_Tij))/sqrt(ntrials),'.');
errorbar(cis', ci_Tij, ci_Tij - Plow, Phigh - ci_Tij,'k.');
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
errorbar(cis', ci_mu_i, ci_mu_i - Plow, Phigh - ci_mu_i, 'k.');
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
errorbar(cis', ci_sigma_i, ci_sigma_i - Plow, Phigh - ci_sigma_i, 'k.');
axis square;
set(gca, 'XTick', [0 1]);
set(gca, 'YTick', [0 1]);
%title('\sigma_i');
text(0.1, 0.85, '\sigma_i');

% Export figure.
addpath exportfig
exportfig(gcf, 'confidence-intervals-random.eps', 'color', 'cmyk', 'width', 3.5, 'height', 3.5);
system('epstopdf confidence-intervals-random.eps');
