#!/usr/bin/env python
"""
Generate plots for synthetic three-state force spectroscopy model.

"""

from bhmm import testsystems
from bhmm import BHMM
from bhmm import plots

# Load force data.
from netcdf4 import Dataset
ncfile = Dataset('fiber3-trace011', 'r')
tau = 0.001 # 1 kHz
obs_label = 'force / pN'
time_units = 's' # seconds
o_t = ncfile.variables['force'] # load trace
O = [o_t] # form list of traces

# Initialize MLHMM.
print "Initializing MLHMM..."
nstates = 3
mlhmm = MLHMM(O, nstates)

# Plot.
plots.plot_state_assignments(model, mlhmm.S[0], o_t, time_units=time_unit, obs_label=obs_label, tau=tau, pdf_filename='fiber3-trace11-mlhmm-stateassignments.pdf')

# Initialize BHMM, using MLHMM model as initial model.
print "Initializing BHMM..."
bhmm = BHMM(O, nstates, intial_model=mlhmm)

# Sample models.
models = bhmm.sample(nsamples=10, save_hidden_state_trajectory=False)

# Generate a sample saving a hidden state trajectory.
models = bhmm.sample(nsamples=1, save_hidden_state_trajectory=True)



