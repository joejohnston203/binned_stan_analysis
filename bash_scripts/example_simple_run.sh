#!/bin/bash

virtualenv_script='../venv/bin/activate'

. $virtualenv_script

#yaml_script=./scripts/example_ricochet_rate_analyzer.yaml
yaml_script=./scripts/example_ricochet_rate_spec_analyzer.yaml

python helper_scripts/shape_fakedata_generator.py -c $yaml_script
morpho -c $yaml_script
python ./helper_scripts/non_morpho_plots.py -c $yaml_script
