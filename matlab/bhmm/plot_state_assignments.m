% Plot traces of observables with state assignments.
function plot_state_assignments(data, model, options)

% PARAMETERS
markersize = 5; % size of points in plot
colors = hsv(model.nstates); % colors for states
stddev_alpha = 0.2;
mean_alpha = 0.4;

% Get number of traces.
ntraces = length(data);

% Determine number of plot panels.
%if (ntraces > 2)
%  nx = 2;
%  ny = ceil(ntraces/2);
%else
%  nx = 1;
%  ny = 1;
%end

nx = 1;
ny = ntraces;
  
% Make new figure.
clf;
  
% Plot all traces.
for trace = 1:ntraces
  % Extract trace data.
  o_t = data{trace};
  T = length(o_t);

  % Extract state assignments.
  s_t = model.state_trajectories{trace};  
  
  % Determine data extents.
  xmin = 0;
  xmax = T * options.tau;
  ymin = min(o_t);
  ymax = max(o_t);

  dx = 0.015 * (xmax - xmin); % width of state marker  

  % Plot each state.
  %subplot(ny, nx, trace);
  %subplot('position', [0.05, 0.80 - 0.60*trace/ny, 0.85, 0.80 - 0.60*(trace-1)/ny]);
  subplot('position', [0.05, 0.80 - 0.60*trace/ny, 0.85, 0.60/(ny+1)]);
  box on;
  hold on;
  for state_index = 1:model.nstates
    state = model.states{state_index}; % state data

    % Plot observations belonging to this state.
    indices = find(s_t == state_index);
    h = plot(options.tau * indices, o_t(indices), '.');
    set(h, 'markersize', markersize);
    set(h, 'color', colors(state_index,:));

    % Shade one standard deviation about state mean.
    h = fill([xmin xmin xmax+dx xmax+dx], [state.mu-state.sigma state.mu+state.sigma state.mu+state.sigma state.mu-state.sigma], colors(state_index,:));
    set(h, 'facealpha', stddev_alpha);
    set(h, 'edgealpha', 0);

    % Plot state means.
    h = fill([xmin xmax+dx xmax+dx xmin], state.mu*[1 1 1 1], colors(state_index,:));
    set(h, 'facecolor',  colors(state_index,:));
    set(h, 'edgecolor',  colors(state_index,:));
    set(h, 'facealpha', stddev_alpha);
    set(h, 'edgealpha', stddev_alpha);
    set(h, 'linewidth', 2);

    % Show a triangle on side of plot to indicate state mean.
    %dy = 0.03 * (ymax - ymin); % 
    %h = fill([xmax xmax+dx xmax+dx], [state.mu state.mu+dy state.mu-dy], colors(state_index,:));    
    h = fill([xmax xmax+dx xmax+dx], [state.mu state.mu+state.sigma state.mu-state.sigma], colors(state_index,:));    
    set(h, 'facecolor',  colors(state_index,:));
    set(h, 'edgecolor',  colors(state_index,:));
    %h = text(xmax+dx, state.mu, '\bf\leftarrow', 'fontsize', 10);
    %set(h, 'color',  colors(state_index,:));
    

  end

  if isfield(options, 'hide_observable_ticks')
    if options.hide_observable_ticks
      set(gca, 'YTick', []);
    end
  end
  if isfield(options, 'hide_time_ticks')
    if options.hide_time_ticks
      set(gca, 'XTick', []);
    end
  end

  %axis([xmin xmax ymin ymax]);
  plot([xmax xmax], [ymin ymax], 'k-'); % black line on right of plot
  axis([xmin xmax+dx ymin ymax]);
  box on

  if (trace == ntraces)
    if (isfield(options, 'time_units'))
      xlabel(sprintf('time / %s', options.time_units));  
    else
      xlabel('time');
    end

    if (isfield(options, 'observable_units'))
      ylabel(sprintf('%s / %s', options.observable_name, options.observable_units));
    else
      ylabel(options.observable_name);
    end
  end

  % Histogram
  nbins = 50;
  %subplot('position', [0.9, 0.80 - 0.60*trace/ny, 0.1, 0.80 - 0.60*(trace-1)/ny]);
  subplot('position', [0.9, 0.80 - 0.60*trace/ny, 0.1, 0.60/(ny+1)]);
  box on;
  hold on;
  [N,x] = hist(o_t, nbins);
  dx = x(2) - x(1);
  N = N / sum(N) / dx;
  %plot(N,x,'k-','linewidth',2);
  h = fill([N zeros(size(N))], [x flipud(x)], [0.5 0.5 0.5]);
  set(h, 'edgealpha', 0.0);
  % TODO: Inscribe state Gaussians.
  hold on
  for state_index = 1:model.nstates
    state = model.states{state_index}; % state data
    phi = @(x) model.Pi(state_index) * (2.*pi)^(-1./2.) * state.sigma^(-1) * exp(-0.5 *((x-state.mu)/state.sigma).^2);
    h = plot(phi(x), x, 'k-', 'linewidth', 1.5);
    set(h, 'color', colors(state_index,:));
  end
  % Plot sum of state Gaussians.
  sum_of_phi = 0.0 * x;
  for state_index = 1:model.nstates
    state = model.states{state_index}; % state data
    phi = @(x) model.Pi(state_index) * (2.*pi)^(-1./2.) * state.sigma^(-1) * exp(-0.5 *((x-state.mu)/state.sigma).^2);
    sum_of_phi = sum_of_phi + phi(x);
  end
  h = plot(sum_of_phi, x, 'k-', 'linewidth', 1.5);

  set(gca,'XTick',[]);
  set(gca,'YTick',[]);
  axis([0 max(N) ymin ymax]);
  box on;
end

return
