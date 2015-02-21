"""
Hidden Markov model

"""

import copy
import numpy as np
from numpy.linalg import norm
from bhmm import HMM

class MLHMM(object):
    """
    Maximum likelihood Hidden Markov model (HMM).

    This class is used to fit a maximum-likelihood HMM to data.

    Examples
    --------

    >>> import testsystems
    >>> [model, O, S] = testsystems.generate_synthetic_observations()
    >>> mlhmm = MLHMM(O, model.nstates)
    >>> model = mlhmm.fit()

    """
    def __init__(self, observations, nstates, initial_model=None, reversible=True, verbose=False, output_model_type='gaussian'):
        """Initialize a Bayesian hidden Markov model sampler.

        Parameters
        ----------
        observations : list of numpy arrays representing temporal data
            `observations[i]` is a 1d numpy array corresponding to the observed trajectory index `i`
        nstates : int
            The number of states in the model.
        initial_model : HMM, optional, default=None
            If specified, the given initial model will be used to initialize the BHMM.
            Otherwise, a heuristic scheme is used to generate an initial guess.
        reversible : bool, optional, default=True
            If True, a prior that enforces reversible transition matrices (detailed balance) is used;
            otherwise, a standard  non-reversible prior is used.
        verbose : bool, optional, default=False
            Verbosity flag.
        output_model_type : str, optional, default='gaussian'
            Output model type.  ['gaussian', 'discrete']

        """
        # Store options.
        self.verbose = verbose
        self.reversible = reversible

        # Store the number of states.
        self.nstates = nstates

        # Store a copy of the observations.
        self.observations = copy.deepcopy(observations)

        # Determine number of observation trajectories we have been given.
        self.ntrajectories = len(self.observations)

        if initial_model:
            # Use user-specified initial model, if provided.
            self.model = copy.deepcopy(initial_model)
        else:
            # Generate our own initial model.
            self.model = self._generateInitialModel(output_model_type)

        return

    def fit(self):
        """Fit a maximum-likelihood HMM model.

        Returns
        -------
        model : HMM
            The maximum likelihood HMM model.


        Examples
        --------

        >>> import testsystems
        >>> [model, O, S] = testsystems.generate_synthetic_observations()
        >>> mlhmm = MLHMM(O, model.nstates)
        >>> model = mlhmm.fit()

        """
        # TODO: Perform EM procedure.

        return self.model

    @classmethod
    def _transitionMatrixMLE(cls, Cij, reversible=True, maxits=1000, reltol=1e-5, epsilon=1e-3, verbose=False):
        """Compute maximum likelihood estimator of transition matrix from fractional transition counts.

        Parameters
        ----------
        Cij : np.array with size (nstates,nstates)
            The transition count matrix.
        reversible : bool, optional, default=True
            If True, the reversible transition matrix will be estimated.
        maxits : int, optional, default=1000
            The maximum number of iterations.
        reltol : float, optional, default=1e-5
            The relative tolerance.
        verbose : bool, optional, default=False
            If True, verbose output will be printed.

        References
        ----------
        [1] Prinz JH, Wu H, Sarich M, Keller B, Fischbach M, Held M, Chodera JD, Schuette C, and Noe F.
        Markov models of molecular kinetics: Generation and validation. JCP 134:174105, 2011.

        Examples
        --------

        >>> Cij = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

        Non-reversible estimate.

        >>> Tij = MLHMM._transitionMatrixMLE(Cij, reversible=False)

        Reversible estimate.

        >>> Tij = MLHMM._transitionMatrixMLE(Cij, reversible=True)

        """
        # Determine size of count matrix.
        nstates = Cij.shape[0]

        # Ensure count matrix is double precision.
        Cij = np.array(Cij, dtype=np.float64)

        if reversible:
            Tij = msm.linalg.transition_matrix_MLE_reversible(Cij)
        else:
            Tij = msm.linalg.transition_matrix_MLE_nonreversible(Cij)

        return Tij

    def _generateInitialModel(self, output_model_type):
        """Use a heuristic scheme to generate an initial model.

        Parameters
        ----------
        output_model_type : str, optional, default='gaussian'
            Output model type.  ['gaussian', 'discrete']

        TODO
        ----
        * Replace this with EM or MLHMM procedure from Matlab code.

        """
        nstates = self.nstates

        # Concatenate all observations.
        collected_observations = np.array([], dtype=np.float64)
        for o_t in self.observations:
            collected_observations = np.append(collected_observations, o_t, axis=0)

        # Fit a Gaussian mixture model to obtain emission distributions and state stationary probabilities.
        from sklearn import mixture
        gmm = mixture.GMM(n_components=nstates)
        gmm.fit(collected_observations)
        states = list()
        for state_index in range(nstates):
            state = { 'model' : 'gaussian', 'mu' : gmm.means_[state_index,0], 'sigma' : np.sqrt(gmm.covars_[state_index,0]) }
            states.append(state)

        # Extract stationary distributions.
        Pi = np.zeros([nstates], np.float64)
        Pi[:] = gmm.weights_[:]

        # Compute transition matrix that gives specified Pi.
        Tij = np.tile(Pi, [nstates, 1])

        # Construct simple model.
        model = HMM(nstates, Tij, states)

        # Compute fractional state memberships.
        from scipy.misc import logsumexp
        Nij = np.zeros([nstates, nstates], np.float64)
        for trajectory_index in range(self.ntrajectories):
            o_t = self.observations[trajectory_index] # extract trajectory
            T = o_t.shape[0]
            # Compute log emission probabilities.
            log_p_ti = np.zeros([T,nstates], np.float64)
            for i in range(nstates):
                log_p_ti[:,i] = model.log_emission_probability(i, o_t)
            # Exponentiate and normalize
            # TODO: Account for initial distribution.
            p_ti = np.zeros([T,nstates], np.float64)
            for t in range(T-1):
                p_ti[t,:] = np.exp(log_p_ti[t,:] - logsumexp(log_p_ti[t,:]))
                p_ti[t,:] /= p_ti[t,:].sum()
            # Accumulate fractional transition counts from this trajectory.
            for t in range(T-1):
                Nij[:,:] = Nij[:,:] + p_ti[t,:].T * p_ti[t+1,:]

        # Compute transition matrix maximum likelihood estimate.
        Tij = self._transitionMatrixMLE(Nij, reversible=self.reversible)

        # Update model.
        model = HMM(nstates, Tij, states, reversible=self.reversible)

        return model
