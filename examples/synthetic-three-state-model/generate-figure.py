#!/usr/bin/env python
"""
Generate plots for synthetic three-state force spectroscopy model.

"""

import argparse

import bhmm
from bhmm.util import testsystems
from bhmm.util.analysis import generate_latex_table


# dynamically import plotting tools
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import plots

def run(nstates, nsamples):
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
    estimator = bhmm.MLHMM(O, nstates)

    print "Initial guess:"
    print str(estimator.hmm.output_model)
    print estimator.hmm.transition_matrix
    print estimator.hmm.stationary_distribution

    # Plot initial guess.
    s_t = None
    o_t = O[0]
    plots.plot_state_assignments(estimator.hmm, s_t, o_t, time_units='s', obs_label='force / pN', tau=tau,
                                 pdf_filename='synthetic-three-state-model-guess-nstates'+str(nstates)+'.pdf')

    print "Fitting HMM..."
    mle = estimator.fit()

    # Plot.
    s_t = mle.hidden_state_trajectories[0]
    import numpy as np
    o_t = O[0]
    plots.plot_state_assignments(mle, s_t, o_t, time_units='s', obs_label='force / pN', tau=tau,
                                 pdf_filename='synthetic-three-state-model-mlhmm-nstates'+str(nstates)+'.pdf')

    # Initialize BHMM with MLHMM model.
    print "Sampling models from BHMM..."
    sampler = bhmm.BHMM(O, nstates, initial_model=mle)
    bhmm_models = sampler.sample(nsamples=nsamples, save_hidden_state_trajectory=False)

    # Generate a sample saving a hidden state trajectory.
    final_models = sampler.sample(nsamples=1, save_hidden_state_trajectory=True)

    # Plot final BHMM sample.
    model = final_models[0]
    s_t = model.hidden_state_trajectories[0]
    o_t = O[0]
    plots.plot_state_assignments(model, s_t, o_t, time_units='s', obs_label='force / pN', tau=tau,
                                 pdf_filename='synthetic-three-state-model-bhmm-nstates'+str(nstates)+'.pdf')

    # write latex table with sample statistics
    conf = 0.95
    sampled_hmm = bhmm.SampledGaussianHMM(mle, bhmm_models)
    generate_latex_table(sampled_hmm, conf=conf, dt=1, time_unit='step',
                         caption='Bayesian HMM parameter estimates for synthetic three-state model.',
                         outfile='synthetic-three-state-model-bhmm-statistics.tex')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Maximum-likelihood and Bayesian HMM estimation for synthetic data')
    parser.add_argument('--nstates', default=3, type=int, help='number of states')
    parser.add_argument('--nsamples', default=100, type=int, help='number of samples in Bayesian estimator')
    parser.add_argument('--verbose', dest='verbose', action='store_true', default=True, help='be loud and noisy')
    args = parser.parse_args()

    # be verbose?
    bhmm.config.verbose = args.verbose

    # go
    run(args.nstates, args.nsamples)




