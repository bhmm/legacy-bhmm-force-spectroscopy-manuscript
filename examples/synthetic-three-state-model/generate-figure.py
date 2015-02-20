#!/usr/bin/env python
"""
Generate plots for synthetic three-state force spectroscopy model.

"""

from bhmm import testsystems
from bhmm import BHMM
from bhmm import plots

# Create model.
true_model = testsystems.force_spectroscopy_model()
nstates = true_model.nstates

# Generate synthetic data.
[O, S] = true_model.generate_synthetic_observation_trajectories(ntrajectories=1, length=50000)

# DEBUG
print "synthetic observation trajectories:"
print O
print "Total state visits, min_state, max_state:"
print testsystems.total_state_visits(nstates, S)

# Initialize BHMM.
tau = 0.001 # time interval per observation
bhmm = BHMM(O, nstates)

# Sample models.
models = bhmm.sample(nsamples=10, save_hidden_state_trajectory=False)

# Generate a sample saving a hidden state trajectory.
models = bhmm.sample(nsamples=1, save_hidden_state_trajectory=True)

# Plot.
model = models[0]
s_t = model.hidden_state_trajectories[0]
o_t = O[0]
plots.plot_state_assignments(model, s_t, o_t, time_units='s', obs_label='force / pN', tau=tau, pdf_filename='synthetic-three-state-model.pdf')


