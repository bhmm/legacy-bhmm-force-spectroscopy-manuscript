import numpy as np
import pyemma.msm.io as msmio
import pyemma.msm.analysis as msmana
from bhmm import MLHMM,BHMM

# load observations
o = np.loadtxt('2well_traj_100K.dat', dtype=int)

# hidden states
nstates = 2

# multiple lags
#lags = [300]
#lags = [1,2,5,10,20,50,100,200,300,400,500,600,700,800,900,1000]
lags = [300]
its      = np.zeros((len(lags)))
its_mean = np.zeros((len(lags)))
its_std  = np.zeros((len(lags)))
likelihoods = np.zeros((len(lags)))
for (i,lag) in enumerate(lags):
    # prepare shifted lagged data
    observations = []
    for shift in range(0, lag):
        observations.append(o[shift:][::lag])

    # initial HMM
    em = MLHMM(observations, nstates, kernel='c', output_model_type='discrete')
    em.fit()
    P = em.model.Tij
    its[i] = msmana.timescales(P, tau=lag)[1]
    likelihoods[i] = em.model.likelihood

    # Initialize BHMM, using MLHMM model as initial model.
    print "BHMM for lag ",lag
    bhmm = BHMM(observations, nstates, initial_model=em.model, verbose=True)

    # Sample models.
    nsamples = 10
    models = bhmm.sample(nsamples=nsamples, save_hidden_state_trajectory=False)
    its_sample = np.zeros((nsamples))
    for j in range(nsamples):
        its_sample[j] = msmana.timescales(models[j].Tij, tau=lag)[1]
    print 'ITS samples = ',its_sample
    its_mean[i] = np.mean(its_sample)
    its_std[i] = np.std(its_sample)

print 'Reference:'
P = msmio.read_matrix('2well_P.dat', mode='sparse').toarray()
itsref = msmana.timescales(P, tau=0.1)[1]
print itsref

print 'Resulting timescales:'
for i in range(len(lags)):
    print lags[i], likelihoods[i], its[i], its_mean[i], '+-', its_std[i]

