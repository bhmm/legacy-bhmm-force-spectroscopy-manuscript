#!/usr/bin/env python
"""
Generate plots for synthetic three-state force spectroscopy model.

"""

from bhmm import testsystems
from bhmm import MLHMM, BHMM
from bhmm import plots

# Create model.
true_model = testsystems.force_spectroscopy_model()
nstates = true_model.nstates
tau = 0.001 # time interval per observation

# Generate synthetic data.
print "Generating synthetic data..."
[O, S] = true_model.generate_synthetic_observation_trajectories(ntrajectories=1, length=50000)

# DEBUG
print "synthetic observation trajectories:"
print O
print "Total state visits, min_state, max_state:"
print testsystems.total_state_visits(nstates, S)

# Generate MLHMM.
print "Generating MLHMM..."
mlhmm = MLHMM(O, nstates, verbose=True)

print "Initial guess:"
print str(mlhmm.model.output_model)
print mlhmm.model.Tij
print mlhmm.model.Pi

# Plot final BHMM sample.
model = mlhmm.model
s_t = None
o_t = O[0]
plots.plot_state_assignments(model, s_t, o_t, time_units='s', obs_label='force / pN', tau=tau, pdf_filename='synthetic-three-state-model-guess.pdf')

print "Fitting HMM..."
mlhmm_model = mlhmm.fit()

# Plot.
s_t = mlhmm.hidden_state_trajectories[0]
o_t = O[0]
plots.plot_state_assignments(mlhmm_model, s_t, o_t, time_units='s', obs_label='force / pN', tau=tau, pdf_filename='synthetic-three-state-model-mlhmm.pdf')

# Initialize BHMM with MLHMM model.
print "Sampling models from BHMM..."
bhmm = BHMM(O, nstates, initial_model=mlhmm_model, verbose=True)
bhmm_models = bhmm.sample(nsamples=10, save_hidden_state_trajectory=False)

# Generate a sample saving a hidden state trajectory.
final_models = bhmm.sample(nsamples=1, save_hidden_state_trajectory=True)

# Plot final BHMM sample.
model = final_models[0]
s_t = model.hidden_state_trajectories[0]
o_t = O[0]
plots.plot_state_assignments(model, s_t, o_t, time_units='s', obs_label='force / pN', tau=tau, pdf_filename='synthetic-three-state-model-bhmm.pdf')



