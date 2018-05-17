'''
Script to compare the outputs of Stan and JAGS MCMC runs
'''

import sys

import numpy as np
import matplotlib.pyplot as plt

def main():
    """
    Compare Stan and JAGS runs

    Takes three positional command line arguments

    Args:
        sys.argv[1]: Path to the output file of stan
            that contains information about the
            distributions, often with the name 
            param_distribution_gauss_fits.txt
        sys.argv[2]: Path to the text output file of back2ground
            with information about the JAGS run. Default name is 
            b2g-$fake_data_name'.txt'
        sys.argv[3]: Output path to save data

    Returns:
        None: Creates plots comparing the two
    """
    if len(sys.argv)<4:
        print("Must provide stan fit, jags fit, "+
              "and output directory in command line")
    print("Stan file: %s"%sys.argv[1])
    print("JAGS file: %s"%sys.argv[2])
    print("Output Directory: %s"%sys.argv[3])

    stan_results = np.genfromtxt(sys.argv[1], dtype=None)

    i_line = 0
    i_start = -1
    i_end = -1
    with open(sys.argv[2], 'r') as jfile:
        for line in jfile:
            if 'GAUS FIT RESULTS' in line:
                print("setting i_start")
                i_start = i_line+2
            if i_start>=0 and not line.strip():
                i_end = i_line
                break
            i_line += 1
    jags_results = np.genfromtxt(sys.argv[2], dtype=None,
                                 skip_header=i_start, max_rows=i_end-i_start)

    if not len(stan_results)==len(jags_results):
        print("Different numbers of parameters were used in stan and jags")
        print("stan_results:\n%s"%stan_results)
        print("jags_results:\n%s"%jags_results)
        print("Discarding extra parameters")

    results = "#P Num\tJAGS P Name\tJ Gauss Mean\tJ Gauss Sigma\tJ Discovery\tS Gauss Mean\tS Gauss Sigma\tS Discovery\t Stan P Name\n"
    param_names = ""
    normalized_diffs = []
    n_params = min(len(stan_results), len(jags_results))
    for  i_p in range(n_params):
        jags_name = jags_results[i_p][0]
        jags_mean = jags_results[i_p][2]
        jags_sigma = jags_results[i_p][3]
        jags_disc = jags_results[i_p][1]=='gaus_n'
        stan_mean = stan_results[i_p][4]
        stan_sigma = stan_results[i_p][5]
        stan_disc = stan_results[i_p][6]
        stan_name = stan_results[i_p][1]
        results += "%i\t%s\t%.3e\t%.3e\t%s\t%.3e\t%.3e\t%s\t%s\n"%\
                   (i_p+1,
                    jags_name,
                    jags_mean, jags_sigma,
                    jags_disc,
                    stan_mean, stan_sigma,
                    stan_disc,
                    stan_name)
        param_names += jags_name+", "+stan_name+"\n"
        normalized_diffs.append((jags_mean-stan_mean)/
                                np.sqrt(jags_sigma**2+stan_sigma**2))

    res_file = open(sys.argv[3]+"/comparison.txt", 'w')
    res_file.write(results)
    res_file.close()

    param_file = open(sys.argv[3]+"/params.txt", 'w')
    param_file.write(param_names)
    param_file.close()

    fig = plt.figure()
    plt.plot(range(1, n_params+1), normalized_diffs,
             linestyle="None", marker="o")
    plt.title("JAGS vs Stan, Marginal Distribution Gaussian Fits Comparison")
    plt.xlabel("Parameter Number")
    plt.ylabel("(mu_J-mu_S)/sqrt(sig_J**2+sig_S**2)")
    plt.savefig(sys.argv[3]+"/comparison.pdf")

if __name__=="__main__":
    main()
