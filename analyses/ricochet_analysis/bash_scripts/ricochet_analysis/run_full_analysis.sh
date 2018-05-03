#!/bin/bash

# All runs are done for a signal of 5 events per day, and a
# background of 16, 36, or 50 events per day. 90% certainty level
# sensitivity is extracted by doing runs with a signal of 0 events
# per day, then extracting a signal rate distribution and taking
# the 90\% quantile value as the 90\% certainty level.

# Note that each setup should be run many times with a single chain. This
# is done so that the plots from each chain can be examined,
# and the chain can be discarded if it did not seem to properly converge.
# The results from the remaining chains can then be averaged to obtain
# final results. The number of chains, iterations, and warmup iterations
# can be changed by editing example_ricochet_rate_analyzer.yaml or
# example_ricochet_rate_spec_analyzer.yaml.
# ----------------------------------------------------------------------------------
# Run the 1 dimensional rate only model with just the cns signal
# and a background with flat time dependence. Then fit the fake data
# using two backgrounds, one flat and one with exponential dependence.
# The signal rate is allowed to fluctuate by 5% as discussed in the paper.

./bash_scripts/ricochet_analysis/rate/run_1D_analysis_60pctOn_2back.sh
./bash_scripts/ricochet_analysis/rate/run_1D_analysis_58pctOn_2back.sh
./bash_scripts/ricochet_analysis/rate/run_1D_analysis_80pctOn_2back.sh

./bash_scripts/ricochet_analysis/rate/run_1D_analysis_60pctOn_2back_90pctCL.sh
# ----------------------------------------------------------------------------------
# Now run the 2 dimensional rate+spectrum model. Data is again generated with
# just the cns signal, a flat background, and background with the shape of the
# H3 background in Fig 3 of the paper. That data is then fit with
# four backgrounds: one flat, one with the H3 shape, one with exponential time
# dependence and flat energy dependence, and one with flat time dependence and
# exponential energy dependence. The signal rate is allowed to fluctuate by 5%, and
# the total number of counts is allowed to fluctuate by 15% to account for
# the detector resolution of 15 eV (divided by our bin size of 100 eV).

./bash_scripts/ricochet_analysis/rate_spec/run_2D_analysis_60pctOn_4back.sh

./bash_scripts/ricochet_analysis/rate_spec/run_2D_analysis_60pctOn_4back_90pctCL.sh
# ----------------------------------------------------------------------------------
