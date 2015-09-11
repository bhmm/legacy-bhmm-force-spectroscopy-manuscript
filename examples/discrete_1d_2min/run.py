import numpy as np
import pyemma.msm.io as msmio
import pyemma.msm.analysis as msmana
import bhmm
bhmm.config.verbose=True

# load observations
o = np.loadtxt('2well_traj_100K.dat', dtype=int)

# hidden states
nstates = 2

# multiple lags
lags = [1,5,10,50,100,500,1000]
#lags = [1,2,5,10,20,50,100,200,300,400,500,600,700,800,900,1000]
its      = np.zeros((len(lags)))
its_mean = np.zeros((len(lags)))
its_std  = np.zeros((len(lags)))
likelihoods = np.zeros((len(lags)))
for (i,lag) in enumerate(lags):
    print ("\n========================================================================")
    print ("LAG = ",lag)
    # prepare shifted lagged data
    observations = []
    for shift in range(0, lag):
        observations.append(o[shift:][::lag])

    # initial HMM
    hmm = bhmm.estimate_hmm(observations, nstates, type='discrete')
    its[i] = lag*hmm.timescales
    likelihoods[i] = hmm.likelihood

    # Sample models.
    print "BHMM for lag ",lag
    nsample = 10
    sampled_hmms = bhmm.bayesian_hmm(observations, hmm, nsample=nsample, store_hidden=False)
    print 'sampled timescales: ',sampled_hmms.timescales_samples

    # store sampled timescale moments
    its_mean[i] = lag*sampled_hmms.timescales_mean
    its_std[i] = lag*sampled_hmms.timescales_std

    print ("========================================================================")

print 'Reference:'
P = msmio.read_matrix('2well_P.dat', mode='sparse').toarray()
itsref = msmana.timescales(P, tau=0.1)[1]
print itsref

print 'Resulting timescales:'
for i in range(len(lags)):
    print lags[i], likelihoods[i], its[i], its_mean[i], '+-', its_std[i]

