% Make plots for JDC's talk

% Load data.
load actual_parameters_2;
load trajectory_2;
load data_set_2;

% Clear figure.
clf;

% Determine length of trajectory.
T = length(realstates);

% Plot actual state trajectory.
subplot('position', [0.1 0.70 0.85 0.29]);
plot(1:T, realstates, 'r.');
%xlabel('time');
ylabel('state');
set(gca,'XTick', []);
set(gca, 'YTick', [1 2 3 4]);
axis([1 T 0 5]);

% Plot FRET data
subplot('position', [0.1 0.40 0.85 0.29]);
plot(1:T, fretdata, 'k.');
%xlabel('time');
set(gca,'XTick', []);
ylabel('EFRET');
set(gca,'YTick', [0 1]);
axis([1 T 0 1]);

% Plot inferred trajectories.
subplot('position', [0.1 0.1 0.85 0.29]);
plot(1:T, allstates(:,1000), 'k.');
xlabel('time');
ylabel('state');
set(gca, 'YTick', [1 2 3 4]);
axis([1 T 0 5]);

% Export figure.
exportfig(gcf, 'jdcplot.eps', 'width', 4.0, 'height', 6.0, 'color', 'rgb', 'fontmode', 'fixed');
system('epstopdf jdcplot.eps');

