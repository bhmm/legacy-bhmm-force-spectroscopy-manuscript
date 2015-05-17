# Synthetic three-state model

Three-state system intended to mimic a protein, with
* high-compliance, low-force unfolded state
* moderately compliant low-population intermediate at intermediate force
* low-compliance, high-force folded state

The term "compliance" refers to the width of the force or extension distribution characterizing the state.

The model used here is `bhmm.testsystems.force_spectroscopy_model()`.

We generate a trajectory of 100 000 observations, and characterize the BHMM mean parameter estimate and 95\% confidence intervals for a subset of this trajectory of varying lengths.

Plots and a LaTeX table is automatically generated with results.

## Manifest
* `generate-figure.py` - source to generate the figure
* `synthetic-three-state-model-guess-nstates3.pdf` - plot of initial guess for states
* `synthetic-three-state-model-mlhmm-nstates3.pdf` - plot of MLHMM fit
* `synthetic-three-state-model-bhmm-nstates3.pdf` - plot of a typical BHMM sample (from the end of the sampling procedure)
* `synthetic-three-state-model-bhmm-statistics.tex` - automatically-generated LaTeX table summarizing mean parameters and 95% confidence intervals
