#!/bin/bash

# JagsBkgAnalysis must be installed first, ie the following
# executeables must be stored in the path (eg /home/joe/bin):
# back2ground, MergeData-3.0.pl, PrepareData, SetLaunch-3.0.pl,
# jags2root, MergeModel-3.0.pl, PrepareDataMult

morpho_venv_activate='/home/joe/MIT_Dropbox/research/stat_code/venv_1.4.1/bin/activate'
fake_data_dir='../../../data/fake_data/'
fake_data_name='fake_data_2_sim'
simulations_dir='../../../data/simulations/'

while true; do
    read -p "Run Stan analysis? " yn
    case $yn in
        [Yy]* )
	    . $morpho_venv_activate
	    python ../../../../../config_builders/cuore_config_builder.py \
		   -c stan_inputs/config_builder.yaml
	    echo "Stan and Morpho config files generated. Press any key to continue."
	    read -n 1 -s
	    morpho -c stan_outputs/cuore0_analysis.yaml
	    echo "Finished running Stan"
	    echo
	    deactivate
	    break;;
        [Nn]* )
	    echo "Not running Stan"
	    echo
	    break;;
        * ) echo "Please answer yes or no.";;
    esac
done

while true; do
    read -p "Run JAGS analysis? " yn
    case $yn in
        [Yy]* )
	    CWD=`pwd`
	    cd jags_inputs
	    ls
	    PrepareDataMult ../$fake_data_dir$fake_data_name'.root' \
			    ListMC.txt ../$simulations_dir \
			    M1 M2 M2sum
	    echo "JAGS data prepared. Press any key to continue"
	    read -n 1 -s
	    cd JAGS_$fake_data_name
	    jags LaunchJags.cmd
	    echo "JAGS run. Press any key to continue"
	    back2ground ../../$fake_data_dir$fake_data_name'.root' \
			ListMC_opt.txt SpectraM1.root SpectraM2.root \
			SpectraM2sum.root
	    cd ..
	    if [ ! -d ../jags_outputs ]; then
		mkdir ../jags_outputs
	    fi
	    cp -r JAGS_$fake_data_name ../jags_outputs/
	    cp -r PrepareData-log ../jags_outputs/
	    cd $CWD
	    echo "Finished running JAGS"
	    echo
	    break;;
        [Nn]* )
	    echo "Not running JAGS"
	    echo
	    break;;
        * ) echo "Please answer yes or no.";;
    esac
done

while true; do
    read -p "Run Stan vs JAGS comparison? " yn
    case $yn in
        [Yy]* )
	    if [ ! -d comparison_outputs ]; then
		mkdir comparison_outputs
	    fi
	    . $morpho_venv_activate
	    python ../../compare_stan_jags.py \
		   stan_outputs/plots/param_dists/param_distribution_gauss_fits.txt \
		   jags_outputs/JAGS_$fake_data_name/b2g-$fake_data_name'.txt' \
		   comparison_outputs
	    deactivate
	    echo "Finished Stan vs JAGS comparison"
	    echo
	    break;;
        [Nn]* )
	    echo "Not doing Stan vs JAGS comparison"
	    echo
	    break;;
        * ) echo "Please answer yes or no.";;
    esac
done
