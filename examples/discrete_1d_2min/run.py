import numpy as np
import pyemma.msm.io as msmio
import pyemma.msm.analysis as msmana
from bhmm import MLHMM

# load observations
o = np.loadtxt('2well_traj_100K.dat', dtype=int)

# multiple lags
#lags = [300]
lags = [1,2,5,10,20,50,100,200,300,400,500,600,700,800,900,1000]
its  = np.zeros((len(lags)))
for (i,lag) in enumerate(lags):
    # prepare shifted lagged data
    observations = []
    for shift in range(0, lag):
        observations.append(o[shift:][::lag])

    # initial HMM
    hmm = MLHMM(observations, 2, kernel='c', output_model_type='discrete')
    hmm.fit()
    P = hmm.model.Tij
    its[i] = msmana.timescales(P, tau=lag)[1]

print 'Reference:'
P = msmio.read_matrix('2well_P.dat', mode='sparse').toarray()
itsref = msmana.timescales(P, tau=0.1)[1]
print itsref

print 'Resulting timescales:'
for i in range(len(lags)):
    print lags[i], its[i]
