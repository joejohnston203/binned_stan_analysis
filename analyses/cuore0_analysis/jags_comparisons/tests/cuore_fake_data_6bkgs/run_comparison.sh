#!/bin/bash

# JagsBkgAnalysis must be installed first, ie the following
# executeables must be stored in the path (eg /home/joe/bin):
# back2ground, MergeData-3.0.pl, PrepareData, SetLaunch-3.0.pl,
# jags2root, MergeModel-3.0.pl, PrepareDataMult

morpho_venv_activate='../../../../../../venv_1.4.1/bin/activate'
fake_data_dir='shared_inputs/'
fake_data_name='fake_data_cuore_6_sim'
simulations_dir='../../../data/simulations_cuore/reduced_sims/combined_trees/'

while true; do
    read -p "Run Stan analysis? " yn
    case $yn in
        [Yy]* )
	    . $morpho_venv_activate
	    python ../../../../../config_builders/cuore_config_builder.py \
		   -c stan_inputs/config_builder.yaml
	    echo "Stan and Morpho config files generated. Press any key to continue."
	    read -n 1 -s
	    morpho -c stan_outputs/cuore_fd_lb_analysis.yaml
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
    read -p "Modify stan_outputs/analysis.yaml and just run morpho? " yn
    case $yn in
        [Yy]* )
	    . $morpho_venv_activate
            sed -i 's/  do_preprocessing: true/  do_preprocessing: false/' stan_outputs/cuore_fd_lb_analysis.yaml
	        morpho -c stan_outputs/cuore_fd_lb_analysis.yaml
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
            mkdir jags_outputs
            cp -r jags_inputs/* jags_outputs/
            cp shared_inputs/Alpha* jags_outputs/
            cp shared_inputs/PeakList* jags_outputs/
	    cd jags_outputs
	    ls
	    PrepareDataMult ../$fake_data_dir$fake_data_name'.root' \
			    ListMC.txt ../$simulations_dir \
			    M1L0 M1L1 M2 M2sum
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
		   jags_outputs/JAGS_$fake_data_name/b2g-$fake_data_name'.txt' \
		   stan_outputs/plots/param_dists/param_distribution_gauss_fits.txt \
		   stan_outputs/plots/param_dists/M1_param_fractions.txt \
		   stan_outputs/plots/param_dists/M2_param_fractions.txt \
		   stan_outputs/plots/param_dists/M2Sum_param_fractions.txt \
	           comparison_outputs \
		   shared_inputs/copied_norm_factors.txt
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
