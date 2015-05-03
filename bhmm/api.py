__author__ = 'noe'

import numpy as _np

from util import types as _types
from hmm.generic_hmm import HMM as _HMM
from estimators.maximum_likelihood import MaximumLikelihoodEstimator as _MaximumLikelihoodEstimator
from estimators.bayesian_sampling import BHMM as _BHMM

def _guess_model_type(observations):
    o1 = _np.array(observations[0])

    # CASE: vector of int? Then we want a discrete HMM
    if _types.is_int_vector(o1):
        return 'discrete'

    # CASE: not int type, but everything is an integral number. Then we also go for discrete
    if _np.allclose(o1, _np.round(o1)):
        isintegral = True
        for i in range(1, len(observations)):
            if not _np.allclose(observations[i], _np.round(observations[i])):
                isintegral = False
                break
        if isintegral:
            return 'discrete'

    # CASE: vector of double? Then we want a gaussian
    if _types.is_float_vector(o1):
        return 'gaussian'

    # None of the above? Then we currently do not support this format!
    raise TypeError('Observations is neither sequences of integers nor 1D-sequences of floats. The current version'
                    'does not support your input.')

def _lag_observations(observations, lag):
    # create new trajectories that are subsampled at lag but shifted
    obsnew = []
    for obs in observations:
        for shift in range(0, lag):
            obsnew.append(obs[shift:][::lag])
    return obsnew

def init_hmm(observations, nstates, lag=1, type=None):
#def generate_initial_model(observations, nstates, output_model_type):
    """Use a heuristic scheme to generate an initial model.

    Parameters
    ----------
    observations : list of ndarray((T_i))
        list of arrays of length T_i with observation data
    nstates : int
        The number of states.
    type : str, optional, default=None
        Output model type from [None, 'gaussian', 'discrete']. If None, will automatically select an output
        model type based on the format of observations.

    Examples
    --------

    Generate initial model for a gaussian output model.

    >>> import bhmm
    >>> [model, observations, states] = bhmm.testsystems.generate_synthetic_observations(output_model_type='gaussian')
    >>> initial_model = init_hmm(observations, model.nstates, type='gaussian')

    Generate initial model for a discrete output model.

    >>> import bhmm
    >>> [model, observations, states] = bhmm.testsystems.generate_synthetic_observations(output_model_type='discrete')
    >>> initial_model = init_hmm(observations, model.nstates, type='discrete')

    """
    # select output model type
    if (type is None):
        type = _guess_model_type(observations)

    if type == 'discrete':
        from bhmm.init import discrete
        return discrete.initial_model_discrete(observations, nstates, lag=lag, reversible=True)
    elif type == 'gaussian':
        from bhmm.init import gaussian
        return gaussian.initial_model_gaussian1d(observations, nstates, reversible=True)
    else:
        raise NotImplementedError('output model type '+str(type)+' not yet implemented.')

def gaussian_hmm(P, means, sigmas, pi=None, stationary=True, reversible=True):
    """ Initializes a 1D-Gaussian HMM

    Parameters
    ----------
    P : ndarray(nstates,nstates)
        Hidden transition matrix
    means : ndarray(nstates, )
        Means of Gaussian output distributions
    sigmas : ndarray(nstates, )
        Standard deviations of Gaussian output distributions
    pi : ndarray(nstates, )
        Fixed initial (if stationary=False) or fixed stationary distribution (if stationary=True).
    stationary : bool, optional, default=True
        If True: initial distribution is equal to stationary distribution of transition matrix
    reversible : bool, optional, default=True
        If True: transition matrix will fulfill detailed balance constraints.

    """
    from hmm.gaussian_hmm import GaussianHMM
    from output_models.gaussian import GaussianOutputModel
    # count states
    nstates = _np.array(P).shape[0]
    # initialize output model
    output_model = GaussianOutputModel(nstates, means, sigmas)
    # initialize general HMM
    ghmm = _HMM(P, output_model, Pi=pi, stationary=stationary, reversible=reversible)
    # turn it into a Gaussian HMM
    ghmm = GaussianHMM(ghmm)
    return ghmm

def discrete_hmm(P, pout, pi=None, stationary=True, reversible=True):
    """ Initializes a discrete HMM

    Parameters
    ----------
    P : ndarray(nstates,nstates)
        Hidden transition matrix
    pout : ndarray(nstates,nsymbols)
        Output matrix from hidden states to observable symbols
    pi : ndarray(nstates, )
        Fixed initial (if stationary=False) or fixed stationary distribution (if stationary=True).
    stationary : bool, optional, default=True
        If True: initial distribution is equal to stationary distribution of transition matrix
    reversible : bool, optional, default=True
        If True: transition matrix will fulfill detailed balance constraints.

    """
    from hmm.discrete_hmm import DiscreteHMM
    from output_models.discrete import DiscreteOutputModel
    # initialize output model
    output_model = DiscreteOutputModel(pout)
    # initialize general HMM
    dhmm = _HMM(P, output_model, Pi=pi, stationary=stationary, reversible=reversible)
    # turn it into a Gaussian HMM
    dhmm = DiscreteHMM(dhmm)
    return dhmm

def estimate_hmm(observations, nstates, lag=1, initial_model=None, type=None,
                 reversible=True, stationary=True, p=None, accuracy=1e-3, maxit=1000):
    r""" Estimate maximum-likelihood HMM

    Generic maximum-likelihood estimation of HMMs

    Parameters
    ----------
    observations : list of numpy arrays representing temporal data
        `observations[i]` is a 1d numpy array corresponding to the observed trajectory index `i`
    nstates : int
        The number of states in the model.
    lag : int
        the lag time at which observations should be read
    initial_model : HMM, optional, default=None
        If specified, the given initial model will be used to initialize the BHMM.
        Otherwise, a heuristic scheme is used to generate an initial guess.
    type : str, optional, default=None
        Output model type from [None, 'gaussian', 'discrete']. If None, will automatically select an output
        model type based on the format of observations.
    reversible : bool, optional, default=True
        If True, a prior that enforces reversible transition matrices (detailed balance) is used;
        otherwise, a standard  non-reversible prior is used.
    stationary : bool, optional, default=True
        If True, the initial distribution of hidden states is self-consistently computed as the stationary
        distribution of the transition matrix. If False, it will be estimated from the starting states.
    p : ndarray (nstates), optional, default=None
        Initial or fixed stationary distribution. If given and stationary=True, transition matrices will be
        estimated with the constraint that they have p as their stationary distribution. If given and
        stationary=False, p is the fixed initial distribution of hidden states.
    accuracy : float
        convergence threshold for EM iteration. When two the likelihood does not increase by more than accuracy, the
        iteration is stopped successfully.
    maxit : int
        stopping criterion for EM iteration. When so many iterations are performed without reaching the requested
        accuracy, the iteration is stopped without convergence (a warning is given)

    Return
    ------
    hmm : :class:`HMM <bhmm.hmm.generic_hmm.HMM>`

    """
    # select output model type
    if (type is None):
        type = _guess_model_type(observations)

    if lag > 1:
        observations = _lag_observations(observations, lag)

    # construct estimator
    est = _MaximumLikelihoodEstimator(observations, nstates, initial_model=initial_model, type=type,
                                      reversible=reversible, stationary=stationary, p=p, accuracy=accuracy, maxit=maxit)
    # run
    est.fit()
    # set lag time
    est.hmm._lag = lag
    # return model
    return est.hmm

def bayesian_hmm(observations, estimated_hmm, nsample=100, store_hidden=False):
    r""" Bayesian HMM based on sampling the posterior

    Generic maximum-likelihood estimation of HMMs

    Parameters
    ----------
    observations : list of numpy arrays representing temporal data
        `observations[i]` is a 1d numpy array corresponding to the observed trajectory index `i`
    estimated_hmm : HMM
        HMM estimated from estimate_hmm or initialize_hmm
    nsample : int, optional, default=100
        number of Gibbs sampling steps
    store_hidden : bool, optional, default=False
        store hidden trajectories in sampled HMMs

    Return
    ------
    hmm : :class:`SampledHMM <bhmm.hmm.generic_sampled_hmm.SampledHMM>`

    """
    # construct estimator
    sampler = _BHMM(observations, estimated_hmm.nstates, initial_model=estimated_hmm,
                    reversible=estimated_hmm.is_reversible, transition_matrix_sampling_steps=1000,
                    type=estimated_hmm.output_model.model_type)

    # Sample models.
    sampled_hmms = sampler.sample(nsamples=nsample, save_hidden_state_trajectory=store_hidden)
    # return model
    from bhmm.hmm.generic_sampled_hmm import SampledHMM
    return SampledHMM(estimated_hmm, sampled_hmms)
