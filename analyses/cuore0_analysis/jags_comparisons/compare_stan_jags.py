'''
Script to compare the outputs of Stan and JAGS MCMC runs
'''

import sys

import numpy as np

def main():
    """
    Compare Stan and JAGS runs

    Takes two positional command line arguments

    Args:
        sys.argv[1]: Path to the output file of stan
            that contains information about the
            distributions, often with the name 
            param_distribution_gauss_fits.txt
        sys.argv[2]: Path to the text output file of back2ground
            with information about the JAGS run. Default name is 
            b2g-$fake_data_name'.txt'

    Returns:
        None: Creates plots comparing the two
    """
    print("Stan file: %s"%sys.argv[1])
    print("JAGS file: %s"%sys.argv[2])

if __name__=="__main__":
    main()
