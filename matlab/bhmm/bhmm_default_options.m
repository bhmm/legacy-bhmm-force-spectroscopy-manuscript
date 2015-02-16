function options = bhmm_default_options()
% Generate structure containing default options for Bayesian HMM.
%
% RETURNS
%  options (struct) - structure containing default options for Bayesian HMM
%
% NOTES
%
% The options structure contains the following fields:
%
% updateMethod : {'bayesian', 'maximum-likelihood'} Choice of method to update HMM each iteration

options = struct();

options.updateMethod = 'bayesian'; % method for updates - one of {'bayesian', 'maximum-likelihood'}
options.tau = 1.0; % time between observations

options.maximumIterations = 1000; % maximum number of allowed iterations
options.convergenceTolerance = 1.0e-4; % relative convergence tolerance for EM

% Options for transition matrix estimation.
options.reversible = true; % infer reversible transition matrices
options.diagonally_dominant = true; % enforce diagonally-dominant transition matrices
options.verbosity = 10; % set verbosity level (0 = silent, 1 = minimal, 2 = more output, ..., 10)
options.equilibrium = true; % trajectory data is initially drawn from equilibrium

options.observable_name = 'force'; % label for obserable axis
options.observable_units = 'pN'; % units for observable axis
options.time_units = 's'; % units for time axis

options.sample_states = 0; % controls whether number of hidden states is to be sampled
options.minstates = 2; % minimum number of states if sample_states is set
options.maxstates = 100; % maximum number of states if sample_states is set

% Options for state observables.
options.enforce_state_ordering = true; % enforce states are ordered in increasing order
options.enforce_state_separation = false; % if nonzero, will enforce state means to be separated by specified amount

options.use_java = 1; % use Java acceleration

return

