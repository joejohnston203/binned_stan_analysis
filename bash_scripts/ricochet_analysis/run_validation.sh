#!/bin/bash

# The following scripts perform various validation runs

#----------------------------------------------------------------------------------------
# Reproduce the initially submitted results, which fit with a
# single flat background and 5% signal variation in each bin.

./bash_scripts/ricochet_analysis/validation/run_1D_analysis_60pctOn_1back.sh
./bash_scripts/ricochet_analysis/validation/run_1D_analysis_58pctOn_1back.sh
./bash_scripts/ricochet_analysis/validation/run_1D_analysis_80pctOn_1back.sh
./bash_scripts/ricochet_analysis/validation/run_2D_analysis_60pctOn_1back_noH3.sh
#----------------------------------------------------------------------------------------
# Do 1 dimensional runs to demonstrate that the fluctuations in the data
# work as expected. Exaggerate the fluctuations to 50% to see their effects.
# We expect the 50% global variation to have the biggest effect, because
# it will cause the uncertainties to saturate around 50%. The 50% bin
# fluctations will have a smaller effect, and this effect will be larger
# for the total bin variation because that will fluctuate the total
# number of counts, while the signal bin variation will just fluctuate
# signal counts.

./bash_scripts/ricochet_analysis/validation/run_1D_analysis_60pctOn_1back_50pctGlob.sh
./bash_scripts/ricochet_analysis/validation/run_1D_analysis_60pctOn_1back_50pctTotBin.sh
./bash_scripts/ricochet_analysis/validation/run_1D_analysis_60pctOn_1back_50pctSigBin.sh
#----------------------------------------------------------------------------------------
# Do 1 dimensional runs with various decay constants, in order to demonstrate that
# the results of the analysis do not depend on the decay rate of the background.

./bash_scripts/ricochet_analysis/validation/run_1D_analysis_60pctOn_2back_exp0.27.sh
./bash_scripts/ricochet_analysis/validation/run_1D_analysis_60pctOn_2back_exp0.027.sh
./bash_scripts/ricochet_analysis/validation/run_1D_analysis_60pctOn_2back_exp0.00027.sh
./bash_scripts/ricochet_analysis/validation/run_1D_analysis_60pctOn_2back_exp0.000027.sh
#----------------------------------------------------------------------------------------
