'''
Script to compare the outputs of Stan and JAGS MCMC runs
'''

import sys

import numpy as np
import matplotlib.pyplot as plt

from morpho.utilities.list_plotter import plot_curves

def main():
    """
    Compare Stan and JAGS runs

    Takes six positional command line arguments, and
    an optional seventh argument

    Args:
        sys.argv[1]: Path to the text output file of back2ground
            with information about the JAGS run. Default name is
            b2g-$fake_data_name'.txt'
        sys.argv[2]: Path to the output file of stan
            that contains information about the
            distributions, often with the name 
            param_distribution_gauss_fits.txt
        sys.argv[3 to 5]: Path to output files with fractions each
            parameter accounts for in the reconstruction, for
            M1, M2, and M2sum respectively
        sys.argv[6]: Output path to save data
        sys.argv[7]: Magnitudes used to generate fake data
            (optional)

    Returns:
        None: Creates plots comparing the two
    """
    if len(sys.argv)<7:
        print("Must provide jags fit, stan fit,"+
              "three files with stan param fractions, "+
              "and output directory in command line")
    print("JAGS file: %s"%sys.argv[1])
    print("Stan file: %s"%sys.argv[2])
    print("Stan M1 file: %s"%sys.argv[3])
    print("Stan M2 file: %s"%sys.argv[4])
    print("Stan M2sum file: %s"%sys.argv[5])
    print("Output Directory: %s"%sys.argv[6])

    i_line = 0
    i_gaus_start = -1
    i_gaus_end = -1
    i_m1_start = -1
    i_m1_end = -1
    i_m2_start = -1
    i_m2_end = -1
    i_m2sum_start = -1
    i_m2sum_end = -1
    try:
        with open(sys.argv[1], 'r') as jfile:
            for line in jfile:
                if 'GAUS FIT RESULTS' in line:
                    i_gaus_start = i_line+2
                    for i_temp in range(2):
                        next(jfile)
                        i_line += 1
                if (i_gaus_start>=0 and i_gaus_end<0
                    and not line.strip()):
                    i_gaus_end = i_line

                if '/////////     Multiplicity M1     /////////' in line:
                    i_m1_start = i_line+6
                    for i_temp in range(6):
                        next(jfile)
                        i_line += 1
                if (i_m1_start>=0 and i_m1_end<0
                    and not line.strip()):
                    i_m1_end = i_line

                if '/////////     Multiplicity M2     /////////' in line:
                    i_m2_start = i_line+6
                    for i_temp in range(6):
                        next(jfile)
                        i_line += 1
                if (i_m2_start>=0 and i_m2_end<0
                    and not line.strip()):
                    i_m2_end = i_line

                if '/////////     Multiplicity M2sum     /////////' in line:
                    i_m2sum_start = i_line+6
                    for i_temp in range(6):
                        next(jfile)
                        i_line += 1
                if (i_m2sum_start>=0 and i_m2sum_end<0
                    and not line.strip()):
                    i_m2sum_end = i_line
                i_line += 1
            if i_m2sum_end == -1:
                i_m2sum_end = i_line

        jags_gaus_fit = np.genfromtxt(sys.argv[1], dtype=None,
                                      skip_header=i_gaus_start,
                                      max_rows=i_gaus_end-i_gaus_start)
        jags_m1_fracts = np.genfromtxt(sys.argv[1], dtype=None,
                                       skip_header=i_m1_start,
                                       max_rows=i_m1_end-i_m1_start)
        jags_m2_fracts = np.genfromtxt(sys.argv[1], dtype=None,
                                       skip_header=i_m2_start,
                                       max_rows=i_m2_end-i_m2_start)
        jags_m2sum_fracts = np.genfromtxt(sys.argv[1], dtype=None,
                                       skip_header=i_m2sum_start,
                                       max_rows=i_m2sum_end-i_m2sum_start)
    except Exception as e:
        print("Got exception when opening JAGS results:\n%s"%e)
        print("i_gaus_start=%i, i_gaus_end=%i"%(i_gaus_start, i_gaus_end))
        print("i_m1_start=%i, i_m1_end=%i"%(i_m1_start, i_m1_end))
        print("i_m2_start=%i, i_m2_end=%i"%(i_m2_start, i_m2_end))
        print("i_m2sum_start=%i, i_m2sum_end=%i"%(i_m2sum_start, i_m2sum_end))
        print("exiting")
        return

    stan_gaus_fit = np.genfromtxt(sys.argv[2], dtype=None)
    stan_m1_fracts = np.genfromtxt(sys.argv[3], dtype=None)
    stan_m2_fracts = np.genfromtxt(sys.argv[4], dtype=None)
    stan_m2sum_fracts = np.genfromtxt(sys.argv[5], dtype=None)

    if not len(stan_gaus_fit)==len(jags_gaus_fit):
        print("Different numbers of parameters were used in stan and jags")
        print("stan_gaus_fit:\n%s"%stan_gaus_fit)
        print("jags_gaus_fit:\n%s"%jags_gaus_fit)
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
    param_xvals = [[],[],[]]
    normalized_diffs = [[],[],[]]
    normalized_diffs_labels = ["mu-3sig>0, J and S",
                               "mu-3sig>0, J or S",
                               "mu-3sig<0, J and S"]
    normalized_diffs_flat = []
    n_params = min(len(stan_gaus_fit), len(jags_gaus_fit))
    for  i_p in range(n_params):
        jags_name = jags_gaus_fit[i_p][0][:p_name_col_width-2].ljust(p_name_col_width, '.')
        jags_means.append(jags_gaus_fit[i_p][2])
        jags_sigmas.append(jags_gaus_fit[i_p][3])
        jags_disc = jags_gaus_fit[i_p][1]=='gaus_n'
        stan_means.append(stan_gaus_fit[i_p][4])
        stan_sigmas.append(stan_gaus_fit[i_p][5])
        stan_disc = stan_gaus_fit[i_p][6]
        stan_name = stan_gaus_fit[i_p][1]

        results += "%i\t%s\t%.3e\t%.3e\t%s\t%.3e\t%.3e\t%s\t%s\n"%\
                   (i_p+1,
                    jags_name,
                    jags_means[-1], jags_sigmas[-1],
                    jags_disc,
                    stan_means[-1], stan_sigmas[-1],
                    stan_disc,
                    stan_name)
        param_names += jags_name+", "+stan_name+"\n"
        if jags_disc and stan_disc:
            temp_idx = 0
        elif jags_disc ^ stan_disc:
            temp_idx = 1
        elif (not jags_disc) and (not stan_disc):
            temp_idx = 2
        param_xvals[temp_idx].append(i_p+1)
        normalized_diffs[temp_idx].append((jags_means[-1]-stan_means[-1])/
                                          np.sqrt(jags_sigmas[-1]**2+stan_sigmas[-1]**2))
        normalized_diffs_flat.append((jags_means[-1]-stan_means[-1])/
                                     np.sqrt(jags_sigmas[-1]**2+stan_sigmas[-1]**2))
        

    jags_means = np.array(jags_means)
    jags_sigmas = np.array(jags_sigmas)
    stan_means = np.array(stan_means)
    stan_sigmas = np.array(stan_sigmas)

    res_file = open(sys.argv[6]+"/comparison.txt", 'w')
    res_file.write(results)
    res_file.close()

    #param_file = open(sys.argv[6]+"/params.txt", 'w')
    #param_file.write(param_names)
    #param_file.close()

    fig = plt.figure()
    param_indices = range(1, n_params+1)
    plt.fill_between(param_indices, -3., 3.,
                     color="red")
    plt.fill_between(param_indices, -2., 2.,
                     color="yellow")
    plt.fill_between(param_indices, -1., 1.,
                     color="cyan")
    plt.plot(param_xvals[0], normalized_diffs[0],
             label=normalized_diffs_labels[0],
             color="black",
             linestyle="None", marker="o", markersize=5)
    plt.plot(param_xvals[1], normalized_diffs[1],
             label=normalized_diffs_labels[1],
             color="black",
             linestyle="None", marker="^", markersize=5)
    plt.plot(param_xvals[2], normalized_diffs[2],
             label=normalized_diffs_labels[2],
             color="black",
             linestyle="None", marker="*", markersize=5)
    plt.legend(loc=0)
    plt.title("JAGS vs Stan, Marginal Distribution Gaussian Fits Comparison")
    plt.xlabel("Parameter Number")
    plt.ylabel("(mu_J-mu_S)/sqrt(sig_J**2+sig_S**2)")
    plt.savefig(sys.argv[6]+"/comparison.pdf")


    (pulls_hist, pulls_hist_edges) = \
        np.histogram(normalized_diffs_flat, bins=20, range=(-4, 4))
    curves = []
    curves.append((pulls_hist_edges, pulls_hist, "histo_error",
                   {"yerr":np.sqrt(pulls_hist), "color":"black"}))
    def norm_gauss(x):
        return np.exp(-x**2/2.)/np.sqrt(2.*np.pi)
    norm_gauss = np.vectorize(norm_gauss)
    x_pts = np.linspace(pulls_hist_edges[0], pulls_hist_edges[-1])
    y_pts = norm_gauss(x_pts)*len(normalized_diffs_flat)*\
            (pulls_hist_edges[1]-pulls_hist_edges[0])
    curves.append((x_pts, y_pts, "default",
                   {"color":"red"}))
    plot_curves(curves,
                sys.argv[6]+"/comparison_pulls.pdf",
                plotter="matplotlib",
                xlabel="N Sigma", ylabel="N Counts",
                title="Pull Distribution",
                xlog=False, ylog=False)

    if len(sys.argv)>=8:
        true_mags = np.loadtxt(sys.argv[7])
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
    plt.savefig(sys.argv[6]+"/extracted_normalizations.pdf")

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
        plt.savefig(sys.argv[6]+"/extracted_normalizations_normalized.pdf")

    # Save and plot things related to the fractions
    jags_fractions = [jags_m1_fracts, jags_m2_fracts, jags_m2sum_fracts]
    stan_fractions = [stan_m2_fracts, stan_m2_fracts, stan_m2sum_fracts]
    labels = ["M1", "M2", "M2sum"]

    for i in range(3):
        lab = labels[i]
        jf = jags_fractions[i]
        sf = stan_fractions[i]

        jf_pts = []
        jf_err = []
        sf_pts = []
        sf_err = []

        p_name_col_width = 40
        results = "#P Num  "+\
                  "JAGS P Name".ljust(p_name_col_width)+"\t"+\
                  "J Fraction  \t"+\
                  "J Fract Err \t"+\
                  "S Fraction  \t"+\
                  "S Fract Err\n"
        for j in range(len(jf)):
            jags_name = jags_gaus_fit[j][0][:p_name_col_width-2].ljust(p_name_col_width, '.')
            jf_pts.append(jf[j][4])
            jf_err.append(jf[j][5])
            sf_pts.append(sf[j][4])
            sf_err.append(sf[j][5])
            results += "%i\t%s\t%.3e\t%.3e\t%.3e\t%.3e\n"%\
                       (j+1, jags_name,
                        jf_pts[-1], jf_err[-1],
                        sf_pts[-1], sf_err[-1])

        res_file = open(sys.argv[6]+"/fracions_"+lab+".txt", 'w')
        res_file.write(results)
        res_file.close()


        param_indices = range(1, n_params+1)
        fig = plt.figure()
        plt.errorbar(param_indices, jf_pts, yerr=jf_err,
                     marker='o', color='green',
                     label='JAGS', linestyle='None')
        plt.errorbar(param_indices, sf_pts, yerr=sf_err,
                     marker='^', color='blue',
                     label='Stan', linestyle="None")
        plt.legend(loc=0)
        plt.title("Fraction of Total, JAGS and Stan")
        plt.xlabel("Parameter Number")
        plt.ylabel("Fraction")
        plt.xlim(0, n_params+1)
        plt.savefig(sys.argv[6]+"/fractions_"+lab+".pdf")

if __name__=="__main__":
    main()
