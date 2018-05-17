'''
Script to compare the outputs of Stan and JAGS MCMC runs
'''

import sys

import numpy as np
import matplotlib.pyplot as plt

def main():
    """
    Compare Stan and JAGS runs

    Takes three positional command line arguments, and
    an optional fourth argument

    Args:
        sys.argv[1]: Path to the output file of stan
            that contains information about the
            distributions, often with the name 
            param_distribution_gauss_fits.txt
        sys.argv[2]: Path to the text output file of back2ground
            with information about the JAGS run. Default name is 
            b2g-$fake_data_name'.txt'
        sys.argv[3]: Output path to save data
        sys.argv[4]: Magnitudes used to generate fake data
            (optional)

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
    try:
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
    except Exception as e:
        print("Got exception when opening stan results:\n%s")
        print("i_line%i, i_start=%i, i_end=%i"%
              (i_line, i_start, i_end))
        print("exiting")
        return

    if not len(stan_results)==len(jags_results):
        print("Different numbers of parameters were used in stan and jags")
        print("stan_results:\n%s"%stan_results)
        print("jags_results:\n%s"%jags_results)
        print("Discarding extra parameters")

    p_name_col_width = 40
    results = "#P Num  "+\
              "JAGS P Name".ljust(p_name_col_width)+"\t"+\
              "J Gauss Mean\t"+\
              "J Gauss Sigma\t"+\
              "J Disc\t"+\
              "S Gauss Mean\t"+\
              "S Gauss Sigma\t"+\
              "S Disc\t"+\
              "Stan P Name\n"
    param_names = ""
    jags_means = []
    jags_sigmas = []
    stan_means = []
    stan_sigmas = []
    normalized_diffs = []
    n_params = min(len(stan_results), len(jags_results))
    for  i_p in range(n_params):
        jags_name = jags_results[i_p][0][:p_name_col_width-2].ljust(p_name_col_width, '.')
        jags_means.append(jags_results[i_p][2])
        jags_sigmas.append(jags_results[i_p][3])
        jags_disc = jags_results[i_p][1]=='gaus_n'
        stan_means.append(stan_results[i_p][4])
        stan_sigmas.append(stan_results[i_p][5])
        stan_disc = stan_results[i_p][6]
        stan_name = stan_results[i_p][1]

        results += "%i\t%s\t%.3e\t%.3e\t%s\t%.3e\t%.3e\t%s\t%s\n"%\
                   (i_p+1,
                    jags_name,
                    jags_means[-1], jags_sigmas[-1],
                    jags_disc,
                    stan_means[-1], stan_sigmas[-1],
                    stan_disc,
                    stan_name)
        param_names += jags_name+", "+stan_name+"\n"
        normalized_diffs.append((jags_means[-1]-stan_means[-1])/
                                np.sqrt(jags_sigmas[-1]**2+stan_sigmas[-1]**2))

    jags_means = np.array(jags_means)
    jags_sigmas = np.array(jags_sigmas)
    stan_means = np.array(stan_means)
    stan_sigmas = np.array(stan_sigmas)

    res_file = open(sys.argv[3]+"/comparison.txt", 'w')
    res_file.write(results)
    res_file.close()

    param_file = open(sys.argv[3]+"/params.txt", 'w')
    param_file.write(param_names)
    param_file.close()

    fig = plt.figure()
    param_indices = range(1, n_params+1)
    plt.plot(param_indices, normalized_diffs,
             linestyle="None", marker="o")
    plt.title("JAGS vs Stan, Marginal Distribution Gaussian Fits Comparison")
    plt.xlabel("Parameter Number")
    plt.ylabel("(mu_J-mu_S)/sqrt(sig_J**2+sig_S**2)")
    plt.savefig(sys.argv[3]+"/comparison.pdf")

    if len(sys.argv)>=5:
        true_mags = np.loadtxt(sys.argv[4])
    else:
        true_mags = None

    fig = plt.figure()
    if not true_mags is None:
        plt.scatter(param_indices, true_mags, color='black',
                    label="Real")
    plt.errorbar(param_indices, jags_means, yerr=jags_sigmas,
                 marker='o', color='green',
                 label="JAGS", linestyle="None")
    plt.errorbar(param_indices, stan_means, yerr=stan_sigmas,
                 marker='^', color='blue',
                 label="Stan", linestyle="None")
    plt.legend(loc=0)
    plt.title("Jags and Stan Param Extractions")
    plt.xlabel("Parameter Number")
    plt.ylabel("Normalization Coefficient")
    plt.savefig(sys.argv[3]+"/extracted_normalizations.png")

    if not true_mags is None:
        fig = plt.figure()
        plt.errorbar(param_indices,
                     (np.array(jags_means)-np.array(true_mags))/np.array(jags_sigmas),
                     marker='o', color='green',
                     label="JAGS", linestyle="None")
        plt.errorbar(param_indices,
                     (np.array(stan_means)-np.array(true_mags))/np.array(stan_sigmas),
                     marker='^', color='blue',
                     label="Stan", linestyle="None")
        plt.legend(loc=0)
        plt.title("Jags and Stan Normalized Norm Extractions")
        plt.xlabel("Parameter Number")
        plt.ylabel("(mu_extracted-mu_true)/sigma_extracted")
        plt.savefig(sys.argv[3]+"/extracted_normalizations_normalized.png")

if __name__=="__main__":
    main()
