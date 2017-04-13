/*
* 2D Binned Data - Analyzer
* -----------------------------------------------------
* Author: Joe Johnston <jpj13@mit.edu>
*
* Date: 10 April 2017
*
* Purpose:
*
* Constructs count rate model in Stan for binned 2D data
* Relies on shape files, etc from shape_fake_data_generator.py
*
*/


functions{
    // Load libraries
    //    include=ricochet;
}

data{
    // Info about each independent variable:
    //real lb_0;             // x (first dimension) lower bound
    //real<lower=lb_0> ub_0; // x upper bound
    int nBins_0;           // x number of bins

    //real lb_1;             // y (second dimension) lower bound
    //real<lower=lb_1> ub_1; // y upper bound
    int nBins_1;           // y number of bins
                                             
    // Binned signal and background shapes, normalized so 
    // the expected number of counts in bin [i,j] is
    // sum_k(signal_rate*x_sig_k[i]*y_sig_k[j])
    //   + sum_l(back_rate*x_back_l[i]*y_back_l[j])
    vector[nBins_0] x_sig_1;
    vector[nBins_1] y_sig_1;
    vector[nBins_0] x_back_1;
    vector[nBins_1] y_back_1;

    // Array with number of data events in each bin
    int fake_data[nBins_0,nBins_1];


    // Bounds on analysis parameters that will be extracted
    real signal_lb;                    // Signal lower bound
    real<lower=signal_lb> signal_ub;   // Signal upper bound
    real background_lb;                      // Background event rate lower bound
    real<lower=background_lb> background_ub; // Background event rate upper bound

}

parameters{
    real<lower=signal_lb,upper=signal_ub> signal_rate;   // Signal rate
    real<lower=signal_lb,upper=signal_ub> signal_rate_global;
    vector[nBins_0] log_signal_rate_bins;          // Signal in each bin, gaussian distributed around actual
    real<lower=background_lb,upper=background_ub> background_rate;   // Background event rate
}

transformed parameters {
    real rate;
    real n_counts_recon[nBins_0,nBins_1];
    vector[nBins_0] signal_rate_bins_0;

    // Calculate and store signal rate in the first bin for plotting purposes
    real signal_rate_smeared;
    signal_rate_smeared = exp(log_signal_rate_bins[1]);

    for(i in 1:nBins_0){
    	signal_rate_bins_0[i] = exp(log_signal_rate_bins[i]);
        for(j in 1:nBins_1) {
	      n_counts_recon[i][j] = signal_rate_bins_0[i]*x_sig_1[i]*y_sig_1[j] +
	      			     background_rate*x_back_1[i]*y_back_1[j];
              if (n_counts_recon[i,j]<0.){
      	      	 print(n_counts_recon[i,j]," Bin ", i,"  ", j);
  	      }
	}
    }
}


model {
  increment_log_prob(sum(log_signal_rate_bins));  // Jacobian 
  signal_rate_global ~ normal(signal_rate,0.05*signal_rate);
  signal_rate_bins_0 ~ normal(signal_rate_global,0.15*signal_rate_global);

// Fit data to poisson distribution of inputted generated data points

for (i in 1:nBins_0){
    for(j in 1:nBins_1){
      target += poisson_lpmf(fake_data[i,j] | n_counts_recon[i,j]);
    }
}

}
