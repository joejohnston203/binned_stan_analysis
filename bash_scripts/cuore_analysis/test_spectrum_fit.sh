#!/bin/bash

virtualenv_script='./activate_venv.sh'

. $virtualenv_script

yaml_script=./scripts/example_cuore_spec_analyzer.yaml

python python_scripts/pre_morpho_processing.py -c $yaml_script
morpho -c $yaml_script
python ./python_scripts/post_morpho_plots.py -c $yaml_script
