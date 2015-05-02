#!/usr/bin/env python
"""
Generate plots for synthetic three-state force spectroscopy model.

"""

from bhmm import testsystems
from bhmm import BHMM
from bhmm import MaximumLikelihoodEstimator
from bhmm import plots

# Load force data.
from netCDF4 import Dataset
ncfile = Dataset('rnase-h-d10a-trace47.nc', 'r')
tau = 0.001 # 1 kHz
obs_label = 'force / pN'
time_units = 's' # seconds
o_t = ncfile.variables['force'] # load trace

# copy data
o_t = o_t[:]
O = [o_t] # form list of traces

# Initialize MLHMM.
print "Initializing MLHMM..."
nstates = 6
mlhmm = MaximumLikelihoodEstimator(O, nstates, verbose=True)

# Plot initial guess.
plots.plot_state_assignments(mlhmm.model, None, O[0], time_units=time_units, obs_label=obs_label, tau=tau, pdf_filename='RNAseH_trace47-guess-nstates'+str(nstates)+'-stateassignments.pdf')

# Fit HMM.
mle = mlhmm.fit()

# Plot.
plots.plot_state_assignments(mle, mlhmm.hidden_state_trajectories[0], o_t, time_units=time_units, obs_label=obs_label, tau=tau, pdf_filename='RNAseH_trace47-mlhmm-nstates'+str(nstates)+'-stateassignments.pdf')

# Initialize BHMM, using MLHMM model as initial model.
print "Initializing BHMM..."
bhmm = BHMM(O, nstates, initial_model=mle, verbose=True)

# Sample models.
models = bhmm.sample(nsamples=10, save_hidden_state_trajectory=False)

# Generate a sample saving a hidden state trajectory.
final_models = bhmm.sample(nsamples=1, save_hidden_state_trajectory=True)

# Plot.
model = final_models[0]
s_t = model.hidden_state_trajectories[0]
o_t = O[0]
plots.plot_state_assignments(model, s_t, o_t, time_units=time_units, obs_label=obs_label, tau=tau, pdf_filename='RNAseH_trace47-bhmm-nstates'+str(nstates)+'-stateassignments.pdf')



