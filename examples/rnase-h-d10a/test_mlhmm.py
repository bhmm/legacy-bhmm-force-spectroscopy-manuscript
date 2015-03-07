#!/usr/bin/env python
"""
Generate plots for synthetic three-state force spectroscopy model.

"""

from bhmm import testsystems
from bhmm import BHMM
from bhmm import MLHMM
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
nstates = 3
mlhmm = MLHMM(O, nstates)
mle = mlhmm.fit()

# Plot.
plots.plot_state_assignments(mle, mlhmm.hidden_state_trajectories[0], o_t, time_units=time_units, obs_label=obs_label, tau=tau, pdf_filename='RNAseH_trace47-mlhmm-nstates'+str(nstates)+'-stateassignments.pdf')


