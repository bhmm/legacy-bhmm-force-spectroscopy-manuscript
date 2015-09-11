#!/usr/bin/env python
"""
Generate figures plots and tables for synthetic three-state force spectroscopy model for use in manuscript.

"""

import os, os.path
import numpy as np

import bhmm
from bhmm.util import testsystems
from bhmm.util.analysis import generate_latex_table

# dynamically import plotting tools
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import plots

def generate_synthetic_data(nobservations):
    # Create model.
    true_model = testsystems.force_spectroscopy_model()
    nstates = true_model.nstates

    # Generate synthetic data.
    print "Generating synthetic data..."
    [O, S] = true_model.generate_synthetic_observation_trajectories(ntrajectories=1, length=nobservations)

    # DEBUG
    print "synthetic observation trajectories:"
    print O
    print "Total state visits, min_state, max_state:"
    print testsystems.total_state_visits(nstates, S)

    return [true_model, O, S]

def analyze_data(O, nstates, nsamples=1000, nobservations=None):
    """
    Analyze the data with the specified number of states.

    Parameters
    ----------
    O : numpy.float
        observation trajectory
    nstates : int
        Number of states to use for analysis.
    nsamples : int, optional, default=1000
        Number of iterations to sample from the Bayesian posterior for the BHMM.
    nobservations : int, optional, default=None
        If specified, number of observations to use from O.

    """

    # Time interval.
    tau = 0.001 # time interval (s) for plotting

    # Truncate O to number of observations.
    if nobservations:
        print "Using only %d observations" % nobservations
        O = [ o_t[0:nobservations] for o_t in O ]
    else:
        nobservations = len(O[0])

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
    filename = os.path.join('figures', 'synthetic-three-state-model-guess-nstates%(nstates)d-nobs%(nobservations)d.pdf' % vars())
    plots.plot_state_assignments(estimator.hmm, s_t, o_t, time_units='s', obs_label='force / pN', tau=tau, pdf_filename=filename)

    print "Fitting HMM..."
    mle = estimator.fit()

    # Plot.
    s_t = mle.hidden_state_trajectories[0]
    import numpy as np
    o_t = O[0]
    filename = os.path.join('figures', 'synthetic-three-state-model-mlhmm-nstates%(nstates)d-nobs%(nobservations)d.pdf' % vars())
    plots.plot_state_assignments(mle, s_t, o_t, time_units='s', obs_label='force / pN', tau=tau, pdf_filename=filename)

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
    filename = os.path.join('figures', 'synthetic-three-state-model-bhmm-nstates%(nstates)d-nobs%(nobservations)d.pdf' % vars())
    plots.plot_state_assignments(model, s_t, o_t, time_units='s', obs_label='force / pN', tau=tau, pdf_filename=filename)

    return [mle, bhmm_models]

def compute_rate(Tij, dt):
    """
    Compute rate via the pseudogenerator.

    Properties
    ----------
    Tij : numpy.array of nstates x nstates
        The transition matrix.
    dt : float
        Time interval for transition matrix.

    Returns
    -------
    Kij : numpy.array of nstates x nstates
        Pseudogenerator for transition matrix.

    """
    nstates = Tij.shape[0]
    Kij = Tij - np.eye(nstates)
    Kij /= dt
    return Kij

def generate_latex_table(true_hmm, sampled_hmm_list, conf=0.95, dt=1, time_unit='ms', obs_name='force', obs_units='pN', outfile=None):
    """
    Generate a LaTeX two-column-wide table showing various true and computed properties and uncertainties.

    Parameters
    ----------
    true_model : bhmm.HMM
        True model parameters.
    sampled_hmm_list : list of bhmm.HMM
        List of three sampled HMMs
    conf : float
        confidence interval. Use 0.68 for 1 sigma, 0.95 for 2 sigma etc.

    """

    # confidence interval
    for sampled_hmm in sampled_hmm_list:
        sampled_hmm.set_confidence(conf)
    # dt
    dt = float(dt)
    # nstates
    nstates = sampled_hmm_list[0].nstates

    table = r"""
\begin{table*}
\caption{{\bf Estimated mean model parameters and confidence intervals for synthetic timeseries data}}
\label{table:synthetic-confidence-intervals}
\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lccccc}
\hline
&  &  & \multicolumn{3}{c}{\bf Estimated Model Parameters}  \\ \cline{4-6}
\multicolumn{2}{l}{\bf Property} & \bf True Value & \bf 1 000 observations & \bf 10 000 observations & \bf 100 000 observations\\ \hline
"""
    # Stationary probability.
    for i in range(nstates):
        if (i == 0):
            table += '\t\tEquilibrium probability '
        table += '\t\t& $\pi_{%d}$ & $%0.3f$' % (i+1, true_hmm.stationary_distribution[i])
        for sampled_hmm in sampled_hmm_list:
            p = sampled_hmm.stationary_distribution_mean
            p_lo, p_hi = sampled_hmm.stationary_distribution_conf
            table += ' & $%0.3f_{\:%0.3f}^{\:%0.3f}$ ' % (p[i], p_lo[i], p_hi[i])
        table += ' \\\\' + '\n'
    table += '\t\t\hline' + '\n'

    # Transition probabilities.
    for i in range(nstates):
        for j in range(nstates):
            if (i == 0) and (j==0):
                table += '\t\tTransition probability ($\Delta t = $%s) ' % (str(dt)+' '+time_unit)
            table += '\t\t& $T_{%d%d}$ & $%0.3f$' % (i+1, j+1, true_hmm.transition_matrix[i,j])
            for sampled_hmm in sampled_hmm_list:
                P = sampled_hmm.transition_matrix_mean
                P_lo, P_hi = sampled_hmm.transition_matrix_conf
                table += ' & $%0.3f_{\:%0.3f}^{\:%0.3f}$' % (P[i,j], P_lo[i,j], P_hi[i,j])
            table += ' \\\\' + '\n'
    table += '\t\t\hline' + '\n'
    table += '\t\t\hline' + '\n'

    # Transition rates via pseudogenerator.
    index = 0
    for i in range(nstates):
        for j in range(nstates):
            if (i != j):
                if (index==0):
                    table += '\t\tTransition rate (%s$^{-1}$) ' % time_unit
                Ktrue = compute_rate(true_hmm.transition_matrix, dt)
                table += '\t\t& $k_{%d%d}$ & $%2.3f$' % (i+1, j+1, Ktrue[i,j])
                for sampled_hmm in sampled_hmm_list:
                    P = sampled_hmm.transition_matrix_mean
                    P_lo, P_hi = sampled_hmm.transition_matrix_conf
                    K = compute_rate(P, dt)
                    K_lo = compute_rate(P_lo, dt)
                    K_hi = compute_rate(P_hi, dt)
                    table += ' & $%.3f_{\:%.3f}^{\:%.3f}$' % (K[i,j], K_lo[i,j], K_hi[i,j])
                index += 1
                table += ' \\\\' + '\n'
    table += '\t\t\hline' + '\n'

    # State mean lifetimes.
    for i in range(nstates):
        if (i == 0):
            table += '\t\tState mean lifetime (%s) ' % time_unit
        l = true_hmm.lifetimes
        l *= dt
        table += '\t\t& $t_{%d}$ & $%.3f$' % (i+1, l[i])
        for sampled_hmm in sampled_hmm_list:
            l = sampled_hmm.lifetimes_mean
            l *= dt
            l_lo, l_hi = sampled_hmm.lifetimes_conf
            l_lo *= dt; l_hi *= dt
            table += ' & $%.3f_{\:%.3f}^{\:%.3f}$' % (l[i], l_lo[i], l_hi[i])
        table += ' \\\\' + '\n'
    table += '\t\t\hline' + '\n'

    # State relaxation timescales.
    for i in range(nstates-1):
        if (i == 0):
            table += '\t\tRelaxation time (%s) ' % time_unit
        t = true_hmm.timescales
        t *= dt
        table += '\t\t& $\\tau_{%d}$ & $%.3f$' % (i+1, t[i])
        for sampled_hmm in sampled_hmm_list:
            t = sampled_hmm.timescales_mean
            t *= dt
            t_lo, t_hi = sampled_hmm.timescales_conf
            t_lo *= dt; t_hi *= dt
            table += ' & $%.3f_{\:%.3f}^{\:%.3f}$' % (t[i], t_lo[i], t_hi[i])
        table += ' \\\\' + '\n'
    table += '\t\t\hline' + '\n'

    if True:
        table += '\t\t\hline' + '\n'

        # State mean forces.
        for i in range(nstates):
            if (i == 0):
                table += '\t\tState %s mean (%s) ' % (obs_name, obs_units)
            m = true_hmm.output_model.means
            table += '\t\t& $\mu_{%d}$ & $%.3f$' % (i+1, m[i])
            for sampled_hmm in sampled_hmm_list:
                m = sampled_hmm.means_mean
                m_lo, m_hi = sampled_hmm.means_conf
                table += ' & $%.3f_{\:%.3f}^{\:%.3f}$' % (m[i], m_lo[i], m_hi[i])
            table += ' \\\\' + '\n'
        table += '\t\t\hline' + '\n'

        # State force standard deviations.
        for i in range(nstates):
            if (i == 0):
                table += '\t\tState %s std dev (%s) ' % (obs_name, obs_units)
            s = true_hmm.output_model.sigmas
            table += '\t\t& $s_{%d}$ & $%.3f$' % (i+1, s[i])
            for sampled_hmm in sampled_hmm_list:
                s = sampled_hmm.sigmas_mean
                s_lo, s_hi = sampled_hmm.sigmas_conf
                table += ' & $%.3f_{\:%.3f}^{\:%.3f}$' % (s[i], s_lo[i], s_hi[i])
            table += ' \\\\' + '\n'
        table += '\t\t\hline' + '\n'

    table += r"""\hline
\end{tabular*}
\end{table*}
"""

    # Write to file if desired.
    if outfile is not None:
        f = open(outfile,'w')
        f.write(table)
        f.close()

    return table


if __name__ == "__main__":
    # Set verbosity.
    bhmm.config.verbose = True

    # Make figures directory.
    directory = 'figures'
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Generate data.
    nobservations = 100000
    [true_hmm, O, S] = generate_synthetic_data(nobservations)

    # Analyze data.
    nstates = 3
    [mle1, bhmm_models1] = analyze_data(O, nstates, nobservations=1000)
    [mle2, bhmm_models2] = analyze_data(O, nstates, nobservations=10000)
    [mle3, bhmm_models3] = analyze_data(O, nstates, nobservations=100000)

    # Write latex table with true and inferred sample statistics.
    conf = 0.95 # confidence interval
    sampled_hmm_models = [ bhmm.SampledGaussianHMM(mle1, bhmm_models1), bhmm.SampledGaussianHMM(mle2, bhmm_models2), bhmm.SampledGaussianHMM(mle3, bhmm_models3)]
    generate_latex_table(true_hmm, sampled_hmm_models, conf=conf, dt=1, time_unit='ms',
                         outfile='synthetic-three-state-model-bhmm-statistics.tex')



