#!/usr/bin/env python
"""
Generate plots for synthetic three-state force spectroscopy model.

"""

import bhmm
import argparse
from bhmm.util.analysis import generate_latex_table

# dynamically import plotting tools
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import plots

def run(nstates, nsamples):
    # Load force data.
    from netCDF4 import Dataset
    ncfile = Dataset('fiber3-trace011.nc', 'r')
    tau = 0.001 # 1 kHz
    obs_label = 'force / pN'
    time_units = 's' # seconds
    o_t = ncfile.variables['force'] # load trace

    # force to make a copy because netCDF appears to cause problems
    nsubsample = 50 # subsampling rate
    o_t = o_t[::nsubsample]
    tau *= nsubsample
    # -------------------

    O = [o_t] # form list of traces

    # Initialize MLHMM.
    print "Initializing MLHMM with "+str(nstates)+" states."
    estimator = bhmm.MLHMM(O, nstates)

    # Plot initial guess.
    plots.plot_state_assignments(estimator.hmm, None, O[0], time_units=time_units, obs_label=obs_label, tau=tau,
                                 pdf_filename='fiber3-trace11-guess-stateassignments-nstate'+str(nstates)+'.pdf')

    # Fit MLHMM
    mle = estimator.fit()

    # Plot.
    plots.plot_state_assignments(mle, mle.hidden_state_trajectories[0], o_t, time_units=time_units,
                                 obs_label=obs_label, tau=tau,
                                 pdf_filename='fiber3-trace11-mlhmm-stateassignments-nstate'+str(nstates)+'.pdf')

    # Initialize BHMM, using MLHMM model as initial model.
    print "Initializing BHMM and running with "+str(nsamples)+" samples."
    sampler = bhmm.BHMM(O, nstates, initial_model=mle)

    # Sample models.
    bhmm_models = sampler.sample(nsamples=nsamples, save_hidden_state_trajectory=False)

    # Generate a sample saving a hidden state trajectory.
    final_models = sampler.sample(nsamples=1, save_hidden_state_trajectory=True)

    # Plot.
    model = final_models[0]
    s_t = model.hidden_state_trajectories[0]
    o_t = O[0]
    plots.plot_state_assignments(model, s_t, o_t, time_units=time_units, obs_label=obs_label, tau=tau,
                                 pdf_filename='fiber3-trace11-bhmm-stateassignments-nstate'+str(nstates)+'.pdf')

    # write latex table with sample statistics
    conf = 0.95
    sampled_hmm = bhmm.SampledGaussianHMM(mle, bhmm_models)
    generate_latex_table(sampled_hmm, conf=conf, dt=tau, time_unit='s',
                         caption='BHMM model estimates for RNA hairpin data.',
                         outfile='p5ab-bhmm-statistics-table'+str(nstates)+'.tex')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Maximum-likelihood and Bayesian HMM estimation from p5ab hairpin data')
    parser.add_argument('--nstates', default=3, type=int, help='number of states')
    parser.add_argument('--nsamples', default=100, type=int, help='number of samples in Bayesian estimator')
    parser.add_argument('--verbose', dest='verbose', action='store_true', default=True, help='be loud and noisy')
    args = parser.parse_args()

    # be verbose?
    bhmm.config.verbose = args.verbose

    # go
    run(args.nstates, args.nsamples)
