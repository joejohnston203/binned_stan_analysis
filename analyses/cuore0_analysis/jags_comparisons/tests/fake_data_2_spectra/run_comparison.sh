#!/bin/bash

morpho_venv_activate='/home/joe/MIT_Dropbox/research/stat_code/venv_1.4.1/bin/activate'

# Run stan
. $morpho_venv_activate
python ../../../../../config_builders/cuore_config_builder.py \
       -c stan_inputs/config_builder.yaml
read -p "Stan and Morpho config files generated. Press enter to continue and run Morpho."
morpho -c stan_outputs/cuore0_analysis.yaml
