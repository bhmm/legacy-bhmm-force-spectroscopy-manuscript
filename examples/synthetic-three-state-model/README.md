# Synthetic three-state model

This example uses synthetic data generated from a three-state system intended to mimic a protein.
The three states are:
* high-compliance, low-force unfolded state
* moderately compliant low-population intermediate at intermediate force
* low-compliance, high-force folded state

The term "compliance" refers to the width of the force or extension distribution characterizing the state.

The model used here is defined in `bhmm.testsystems.force_spectroscopy_model()`.

Plots and a LaTeX table is automatically generated.

## Manifest
* `synthetic-example.py` - command-line executable to run to explore dependence of analysis on number of states and trajectory length
* `generate-figure.py` - source to generate the figure used in paper

## Usage

To use the command-line executable:
```
[LSKI1497:bhmm.jchodera/examples/synthetic-three-state-model] choderaj% python synthetic-example.py --help
usage: synthetic-example.py [-h] [--nstates NSTATES] [--nsamples NSAMPLES]
                            [--verbose]

Maximum-likelihood and Bayesian HMM estimation for synthetic data

optional arguments:
  -h, --help           show this help message and exit
  --nstates NSTATES    number of states
  --nsamples NSAMPLES  number of samples in Bayesian estimator
  --verbose            be loud and noisy
```
