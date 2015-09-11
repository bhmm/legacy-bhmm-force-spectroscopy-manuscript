% Plot traces of observables with state assignments.
function plot_fractional_state_assignments(data, models, options)

% PARAMETERS
markersize = 5; % size of points in plot
colors = hsv(models(1).nstates); % colors for states
alpha = 0.2;

% Get number of traces.
ntraces = length(data);

% Determine number of plot panels.
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

  % Plot each state.
  subplot(ny, nx, trace);
  hold on;
  for state_index = 1:model.nstates
    state = model.states{state_index}; % state data
    indices = find(s_t == state_index);
    h = plot(options.tau * indices, o_t(indices), '.');
    set(h, 'markersize', markersize);
    set(h, 'color', colors(state_index,:));

    % Color states    
    h = fill([xmin xmin xmax xmax], [state.mu-state.sigma state.mu+state.sigma state.mu+state.sigma state.mu-state.sigma], colors(state_index,:));
    set(h, 'facealpha', alpha);
    set(h, 'edgealpha', 0);

  end
  set(gca, 'YTick', []);
  set(gca, 'XTick', []);

  axis([xmin xmax ymin ymax]);
  box on;


  if (trace == ntraces)
    xlabel('time');  
    ylabel('obs');
  end
end

return
