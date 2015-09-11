function model_samples = bhmm(data, model, nsamples, options)
% Bayesian hidden Markov model (HMM) analysis of single-molecule data
%
% The data is assumed to come from a fixed number of states each emitting real-valued signals corresponding to 
% one-dimensional i.i.d. normal random variates.  Each state is characterized by a mean and variance, inferred 
% from the data.  Transitions between the states are governed by a transition matrix, and dynamics is assumed
% to be Markovian.
%
% model_samples = bhmm(data, model, nsamples)
%
% ARGUMENTS
%   data (cell array of 1D arrays) - observed trajectories of some real-valued signal
%   model (structure) - initial model parameters
%   nsamples (integer) - number of (potentially correlated) samples to generate
%   options (struct) - options for run - see bhmm_default_options()
%
% RETURNS
%   model_samples (array of structures) - (potentially correlated) samples of model parameters
%
% NOTES
%   Trajectories provided in 'data' may be of different lengths. 
%   The model must have already been initialized by a call to mlhmm or bhmm_guess.

% Import fast Java helper code.
import bhmm_helper;

if (options.verbosity >= 1)
  disp('******************************************************************************');
  disp('Sampling from BHMM posterior using Gibbs sampling...');
end

% Generate model samples.
for sample = 1:nsamples
  
  if (options.verbosity >= 2)
    disp('******************************************************************************');
    disp(sprintf('Bayesian HMM iteration %d', sample));
  end
  
  %
  % Update model using Gibbs sampling procedure.
  % TODO: We should add an inner loop to allow the samples to decorrelate between storage.
  %
  
  % Update state trajectories.
  if (options.verbosity >= 3)
    disp('Updating state trajectories...');
  end
  model = update_state_trajectories(data, model, options);

  % Update state emission probability densities given the current state trajectories.
  if (options.verbosity >= 3)
    disp('Updating state emission probabilities...');  
  end
  model = update_emission_probabilities(data, model, options);

  % Update transition matrix given observed transitions.
  if (options.verbosity >= 3)
    disp('Updating transition rates...');
  end
  model = update_transition_matrix(data, model, options);

  % Update number of hidden states, if requested.
  if (options.verbosity >= 3)
    disp('Updating number of hidden states...');
  end
  if (options.sample_states)
    % Attempt reversible-jump move in number of states.
    model = update_number_of_states(data, model, options)
  end  

  % Store updated model sample.
  model_samples(sample) = model;

  % DEBUG
  if (options.verbosity >= 3)
    disp('Model summary:');
    show_model(model);
  end
end


if (options.verbosity >= 1)
  disp('******************************************************************************');
end
