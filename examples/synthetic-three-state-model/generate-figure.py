#!/usr/bin/env python
"""
Generate plots for synthetic three-state force spectroscopy model.

"""

from bhmm import testsystems
from bhmm import BHMM

# Create model.
true_model = testsystems.force_spectroscopy_model()

# Generate synthetic data.
[O, S] = true_model.generate_synthetic_observation_trajectories(ntrajectories=1, length=50000)



# Initialize BHMM.
tau = 0.001 # time interval per observation
nstates = 3
bhmm = BHMM(nstates, O)

# Sample models.
models = bhmm.sample(nsamples=10, save_hidden_state_trajectory=False)

# Extract hidden state trajectories and observations.
s_t = S[0]
o_t = O[0]

# Generate a sample saving a hidden state trajectory.
models = bhmm.sample(nsamples=1, save_hidden_state_trajectory=True)

# Plot.
plot_state_assignments(model, s_t, o_t, time_units='s', obs_label='force / pN', tau=tau, pdf_filename='synthetic-three-state-model.pdf')


