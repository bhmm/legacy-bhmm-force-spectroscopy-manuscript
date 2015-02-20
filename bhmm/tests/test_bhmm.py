#!/usr/local/bin/env python

"""
Test BHMM using simple analytical models.

"""

from functools import partial
from bhmm import testsystems
from bhmm import BHMM

def run_bhmm(nstates):
    """
    Run the BHMM on synthetic data with the given number of states.

    Parameters
    ----------
    nstates : int
        The number of states to test the BHMM with.

    """
    # Generate synthetic observations.
    [model, S, O] = model.generate_synthetic_observation_trajectories(nstates=nstates)
    # Initialize a BHMM.
    bhmm = bhmm.BHMM(O, nstates)
    # Sample from the posterior.
    models = bhmm.sample(nsamples=10)

    return

def test_bhmm_synthetic():
    """
    Test the BHMM model on synthetic datasets.

    """

    nstates_min = 2 # minimum number of states to test
    nstates_max = 8 # maximum number of states to test

    for nstates in range(nstates_min, nstates_max):
        f = partial(run_bhmm, nstates)
        f.description = "Testing BHMM on synthetic data for %d states"
        yield f

    return

