function model = update_number_of_states(data, model, options)
% Update number of hidden states using reversible-jump scheme.
%
% model = update_number_of_states(data, model, options)
%
% ARGUMENTS
%   data (cell array of 1D arrays) - observed trajectories of some real-valued signal
%   model (structure) - model parameters, containing state trajectories
%   options (structure) - options structure
%
% RETURNS
%   model (structure) - model with updated number of states
%
% NOTES
%   The hidden state assignment trajctory must be redrawn immediately afterwards.
%   If options.minstates is specified, the number of states will not be allowed to fall below this.
%   If options.maxstates is specified, the number of states will not be allowed to rise above this.
%
% TODO
% - Instead of just rejecting moves that create too few or too many states, only propose
%   new models with valid numbers of states.

% DEBUG
options.verbosity = 10;

% Do nothing if sample_states is not set, or if minstates == maxstates.
if (options.sample_states == 0) || (options.minstates == options.maxstates)
  return
end

% DEBUG
show_model(model)

% Choose split or merge move.
% Here, we split or merge with equal probability, unless we're at min or max states.
move_type = '';
if (model.nstates == options.maxstates)
  move_type = 'merge';
elseif (model.nstates == options.minstates)
  move_type = 'split';
elseif (rand < 0.5)
  move_type = 'merge';
else
  move_type = 'split';
end

% Compute model likelihood, ignoring hidden state assignments.
if (options.verbosity >= 5)
  disp('Computing log-likelihood...');
end
log_likelihood = compute_log_likelihood(model, data, options, 0)
if (options.verbosity >= 6)
  disp(sprintf('log P(O | \\Theta_old) = %f', log_likelihood));
end

% Attempt move
if (move_type == 'merge')
  %
  % Merge model states a and b into newmodel state i.
  %
  % Create a proposed model.
  newmodel = model;
  % Decrement number of states.
  newmodel.nstates = model.nstates - 1;
  % Select states to merge.
  if options.enforce_state_ordering
    a = randi(model.nstates-1);
    b = a+1;
  else
    P = randperm(model.nstates);
    a = P(1); % state to survive
    b = P(2); % state to be removed
  end
  if (options.verbosity >= 4) 
    disp(sprintf('Attempting to merge states %d and %d', a, b));
  end
  
  % Construct state indices for mapping.
  source_states = setdiff(1:model.nstates, [a b]);
  dest_states = setdiff(1:newmodel.nstates, [newmodel.nstates]);
  i = newmodel.nstates; % (a b) in old model become state (i) in new model
  
  % Merge probabilities.
  newmodel.Pi = zeros(newmodel.nstates,1);
  newmodel.Pi(dest_states) = model.Pi(source_states);
  newmodel.Pi(i) = model.Pi(a) + model.Pi(b);

  % Determine fraction coming from states a and b.
  f_a = model.Pi(a) / newmodel.Pi(i);
  f_b = model.Pi(b) / newmodel.Pi(i);

  % Construct new state emission parameters.
  newmodel.states = cell(newmodel.nstates,1);
  for index = 1:length(dest_states)
    newmodel.states{dest_states(index)} = model.states{source_states(index)};
  end
  newmodel.states{i} = model.states{a};
  newmodel.states{i}.mu = f_a*model.states{a}.mu + f_b*model.states{b}.mu;
  newmodel.states{i}.sigma = sqrt( f_a*(model.states{a}.mu^2 + model.states{a}.sigma^2) + f_b*(model.states{b}.mu^2 + model.states{b}.sigma^2) - newmodel.states{i}.mu^2 );
  
  % Construct merged transition matrix.
  newmodel.Tij = zeros(newmodel.nstates,newmodel.nstates);
  newmodel.Tij(dest_states,dest_states) = model.Tij(source_states,source_states); % transitions among unchanged states
  newmodel.Tij(i,dest_states) = f_a*model.Tij(a,source_states) + f_b*model.Tij(b,source_states); % weighted sum of probability out of states
  newmodel.Tij(dest_states,i) = model.Tij(source_states,a) + model.Tij(source_states,b); % total sum of probability going into states a and b
  newmodel.Tij(i,i) = 1 - sum(newmodel.Tij(i,dest_states)); % row normalization

  % Compute r_k and s_k required to unmerge.
  r_k = model.Tij(source_states,a) ./ newmodel.Tij(dest_states,i);
  s_k = 1 - r_k;
  % TODO: Compute r_ab.

  % DEBUG: Sanity check
  if (any(r_k > 1))
    show_model(model)
    show_model(newmodel)
    r_k
    error('r_k out of range');
  end

  % Compute proposal probability ratio.
  log_proposal_ratio = log(model.nstates) - log(2*newmodel.nstates) - 0.5*log(2*pi) - (model.states{a}.mu-newmodel.states{i}.mu)^2/(newmodel.states{i}.sigma^2) + log(newmodel.states{i}.sigma) - (model.states{a}.sigma/newmodel.states{i}.sigma) + sum(f_a*log(r_k) + f_b*log(s_k)) - newmodel.nstates*betaln(1+f_a,1+f_b);
elseif (move_type == 'split')
  % 
  % Split model state i into newmodel states a and b.
  %
  if (options.verbosity >= 4)
    disp('Attempting split move...');
  end
  % Create a proposed model.
  newmodel = model;
  % Increment number of states.
  newmodel.nstates = model.nstates + 1;
  % Select state to split.
  i = randi(model.nstates); % state to split
  
  % Construct state indices for mapping.
  source_states = setdiff(1:model.nstates, [i]);
  dest_states = 1:(newmodel.nstates-2);
  a = newmodel.nstates-1; % state to split to 
  b = newmodel.nstates; % state to split to
  
  % Determine fraction coming from states a and b.
  f_a = rand(); % fraction of probability split from i -> a
  f_b = 1.0 - f_a; % fraction of probability split from i -> b

  % Split probabilities.
  newmodel.Pi = zeros(newmodel.nstates,1);
  newmodel.Pi(dest_states) = model.Pi(source_states);
  newmodel.Pi(a) = f_a*model.Pi(i);
  newmodel.Pi(b) = f_b*model.Pi(i);

  % Construct new state emission parameters.
  newmodel.states = cell(newmodel.nstates,1);
  for index = 1:length(dest_states)
    newmodel.states{dest_states(index)} = model.states{source_states(index)};
  end
  newmodel.states{a} = model.states{i};
  newmodel.states{a}.mu = model.states{i}.sigma*randn() + model.states{i}.mu;
  newmodel.states{a}.sigma = exprnd(1.0 / model.states{i}.sigma);
  newmodel.states{b} = model.states{i};
  newmodel.states{b}.mu = (1.0/f_b)*model.states{i}.mu - (f_a/f_b)*newmodel.states{a}.mu;
  newmodel.states{b}.sigma = sqrt( (1.0/f_b)*(model.states{i}.mu^2 + model.states{i}.sigma^2) - (f_a/f_b)*(newmodel.states{a}.mu^2 + newmodel.states{b}.sigma^2) - newmodel.states{b}.mu^2 );
  
  % Construct split transition matrix.
  newmodel.Tij = zeros(newmodel.nstates,newmodel.nstates);
  newmodel.Tij(dest_states,dest_states) = model.Tij(source_states,source_states);
  r_k = betarnd((1+f_a)*ones(model.nstates-1,1), (1+f_b)*ones(model.nstates-1,1));
  s_k = 1 - r_k;
  r_ab = rand;
  r_ab = rand * min(1.0 - sum(newmodel.Tij(a,dest_states)), (1.0 - sum(newmodel.Tij(b,dest_states))) * (newmodel.Pi(b)/newmodel.Pi(a)))
  newmodel.Tij(dest_states,a) = r_k .* model.Tij(source_states,i);
  newmodel.Tij(dest_states,b) = s_k .* model.Tij(source_states,i);
  newmodel.Tij(a,dest_states) = (newmodel.Pi(dest_states)/newmodel.Pi(a))' .* newmodel.Tij(dest_states,a)';
  newmodel.Tij(b,dest_states) = (newmodel.Pi(dest_states)/newmodel.Pi(b))' .* newmodel.Tij(dest_states,b)';
  newmodel.Tij(a,b) = r_ab * (1.0 - sum(newmodel.Tij(a,dest_states)))
  newmodel.Tij(b,a) = (newmodel.Pi(a)/newmodel.Pi(b))*newmodel.Tij(a,b);
  newmodel.Tij(a,a) = 1.0 - sum(newmodel.Tij(a,dest_states)) - newmodel.Tij(a,b);
  newmodel.Tij(b,b) = 1.0 - sum(newmodel.Tij(b,dest_states)) - newmodel.Tij(b,a);

  % Compute proposal probability ratio.
  %log_proposal_ratio = log(2*model.nstates) - log(newmodel.nstates);
  log_proposal_ratio = log(2*model.nstates) - log(newmodel.nstates) + 0.5*log(2*pi) + (model.states{a}.mu-newmodel.states{i}.mu)^2/(newmodel.states{i}.sigma^2) - log(newmodel.states{i}.sigma) + (model.states{a}.sigma/newmodel.states{i}.sigma) - sum(f_a*log(r_k) - f_b*log(s_k)) + model.nstates*betaln(1+f_a,1+f_b);
elseif (mode == 'birth')
  %
  % Birth mode: Add a new state.
  %
  if (options.verbosity >= 4)
    disp('Attempting birth move...');
  end
  % Create a proposed model.
  newmodel = model;
  % Increment number of states.
  newmodel.nstates = model.nstates + 1;

  % Construct state indices for mapping.
  source_states = 1:model.nstates;
  dest_states = source_states;

  % Split probabilities.
  newmodel.Pi = zeros(newmodel.nstates,1);
  rP = rand; sP = 1.0 - rP;
  newmodel.Pi(newmodel.nstates) = rP;
  newmodel.Pi(dest_states) = sP * model.Pi(source_states); 

  % Construct new state emission parameters.
  newmodel.states = cell(newmodel.nstates,1);
  for index = 1:length(dest_states)
    newmodel.states{dest_states(index)} = model.states{source_states(index)};
  end
  i = randi(model.nstates);
  mu = model.states{i}.sigma * randn + model.states{i}.mu;
  sigma = exprnd(model.states{i}.sigma^(-1));
  newmodel.states{newmodel.nstates} = model.states{1};
  newmodel.states{newmodel.nstates}.mu = mu;
  newmodel.states{newmodel.nstates}.sigma = sigma;

  % Construct split transition matrix.
  newmodel.Tij = zeros(newmodel.nstates,newmodel.nstates);
  newmodel.Tij(dest_states,dest_states) = model.Tij(source_states,source_states);
  newmodel.Tij(newmodel.nstates,dest_states) = rand(1,model.nstates-1,1);
  newmodel.Tij(dest_states,newmodel.nstates) = (newmodel.Pi(newmodel.nstates)./newmodel.Pi(dest_states)).*newmodel.Tij(newmodel.nstates,dest_states);
  newmodel.Tij(newmodel.nstates,newmodel.nstates) = 1.0;  


elseif (mode == 'death')
  %
  % Death mode: Remove a state.
  %
  if (options.verbosity >= 4)
    disp('Attempting death move...');
  end
  % Create a proposed model.
  newmodel = model;
  % Increment number of states.
  newmodel.nstates = model.nstates - 1;
  % Choose which state to eliminate.
  i = randi(model.nstates);

  % Construct state indices for mapping.
  source_states = setdiff(1:model.nstates, [i]);
  dest_states = 1:newmodel.nstates;

  % Assign new state probabilities.
  newmodel.Pi = zeros(newmodel.nstates,1);
  rP = model.Pi(i); sP = 1.0 - rP;
  newmodel.Pi(dest_states) = model.Pi(source_states) / sP; 

  % Construct new state emission parameters.
  newmodel.states = cell(newmodel.nstates,1);
  for index = 1:length(dest_states)
    newmodel.states{dest_states(index)} = model.states{source_states(index)};
  end

  % Compute probability of creating killed state from remaining states.  
  log_proposal_probability_k = zeros(newmodel.nstates,1);
  for k = 1:newmodel.nstates
    log_proposal_probability_k(k) = log(newmodel.Pi(k)) - 0.5*log(2*pi) - 0.5*(model.states{i}.mu-newmodel.states{k}.mu)^2/(newmodel.states{k}.sigma^2) - log(newmodel.states{k}.sigma) - (model.states{i}.sigma/newmodel.states{k}.sigma); 
  end
  log_proposal_probability = log(0.5) - log(model.nstates) + logsum(log_proposal_probability_k);

  % Construct transition matrix.
  newmodel.Tij = zeros(newmodel.nstates,newmodel.nstates);
  newmodel.Tij(dest_states,dest_states) = model.Tij(source_states,source_states);
  for k = 1:newmodel.nstates
    newmodel.Tij(k,:) = newmodel.Tij(k,:) / sum(newmodel.Tij(k,:));
  end
  % TODO: Check to make sure detailed balance is still satisfied and Pi is still accurate.
end

% Catch problematic transition matrices.
if(any(any(newmodel.Tij < 0)) || any(any(newmodel.Tij > 1)))
%  show_model(model)
%  show_model(newmodel)
  disp('Tij out of range')
  return
end

% Catch problematic variances.
for k = 1:newmodel.nstates
  if (newmodel.states{k}.sigma^2 <= 0.0)
%    show_model(model)
%    show_model(newmodel)
    disp(sprintf('sigma^2 is negative for state %d', k));
    return
  end
end

% Bookkeeping updates.
newmodel.logTij = log(newmodel.Tij);
newmodel.Pi = stationary_probability(newmodel.Tij);
newmodel.logPi = log(newmodel.Pi);

% Compute proposed model likelihood, ignoring hidden state assignments.
if (options.verbosity >= 5)
  disp('Computing log-likelihood...');
end
new_log_likelihood = compute_log_likelihood(model, data, options, 0)
if (options.verbosity >= 6)
  disp(sprintf('log P(O | \\Theta_new) = %f', new_log_likelihood));
end

% Compute log prior ratio.
log_prior = 0.0;
for i = 1:model.nstates
  log_prior = log_prior - log(model.states{i}.sigma);
end
new_log_prior = 0.0;
for i = 1:newmodel.nstates
  new_log_prior = new_log_prior - log(newmodel.states{i}.sigma);
end
log_prior_ratio = (new_log_prior - log_prior);

% Compute log acceptance probability.
log_acceptance_probability = (new_log_likelihood - log_likelihood) + log_proposal_ratio + log_prior_ratio;

disp('model')
show_model(model)
disp('newmodel')
show_model(newmodel)

if options.enforce_state_ordering
  % TODO: Re-order states.
end

% Accept or reject.
disp(sprintf('log_acceptance_probability = %f', log_acceptance_probability));
if (log_acceptance_probability >= 0.0) || (rand < exp(log_acceptance_probability))
  % Accept.
  %options.number_of_state_moves_accepted = options.number_of_state_moves_accepted + 1;
  disp('Accepted');
  % Save new model.
  model = newmodel;
else
  disp('Rejected');
end

% Always draw new state trajectories after attempt.
if (options.verbosity >= 5)
  disp('Drawing new state trajectories...');
end
model = update_state_trajectories(data, model, options);

if (options.verbosity >= 5)
  disp('State update move completed.');
end
return

