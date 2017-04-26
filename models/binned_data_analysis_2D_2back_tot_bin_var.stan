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
  int nBins_0;           // x number of bins
  int nBins_1;           // y number of bins
                                            
  // Binned signal and background shapes, normalized so 
  // the expected number of counts in bin [i,j] is
  // sum_k(signal_rate*x_sig_k[i]*y_sig_k[j])
  //   + sum_l(back_rate*x_back_l[i]*y_back_l[j])
  vector[nBins_0] x_sig_1;
  vector[nBins_1] y_sig_1;
  vector[nBins_0] x_back_1;
  vector[nBins_1] y_back_1;
  vector[nBins_0] x_back_2;
  vector[nBins_1] y_back_2;
  // Array with number of data events in each bin
  int fake_data[nBins_0,nBins_1];

  // Bounds on analysis parameters that will be extracted
  real signal_lb;                    // Signal lower bound
  real<lower=signal_lb> signal_ub;   // Signal upper bound
  real background_lb;                      // Background event rate lower bound
  real<lower=background_lb> background_ub; // Background event rate upper bound

  // Info about gaussian fluctuations in data
  real<lower=0.0> sig_x_gauss_bin_frac; // Fraction the signal can vary in each bin
  real<lower=0.0> sig_x_gauss_global_frac; // Fraction the overall signal can vary
  real<lower=0.0> total_bin_variation; // Fraction the total counts is varied in each bin

}

parameters{
  real<lower=signal_lb,upper=signal_ub> signal_rate;   // Signal rate
  real<lower=signal_lb,upper=signal_ub> signal_rate_global_gauss;
  vector[nBins_0] log_signal_rate_bins_gauss;          // Signal in each bin, gaussian distributed around actual
  real log_n_counts_recon_envelope[nBins_0,nBins_1]; // gaussian envelope in each bin for n_counts_recon
  real<lower=background_lb,upper=background_ub> background_rate_1;
  real<lower=background_lb,upper=background_ub> background_rate_2;
}

transformed parameters {
  real signal_rate_global;
  real n_counts_recon[nBins_0,nBins_1];
  vector[nBins_0] signal_rate_bins_gauss;
  vector[nBins_0] signal_rate_bins;
  real n_counts_recon_envelope[nBins_0,nBins_1];
  real background_rate;

  // Calculate and store signal rate in the first bin for plotting purposes
  real signal_rate_smeared;

  if(sig_x_gauss_global_frac==0.0){
    signal_rate_global = signal_rate;
  }else{
    signal_rate_global = signal_rate_global_gauss;
  }

  background_rate = background_rate_1+background_rate_2;

  for(i in 1:nBins_0){
    signal_rate_bins_gauss[i] = exp(log_signal_rate_bins_gauss[i]);
    if(sig_x_gauss_bin_frac==0.0){
      signal_rate_bins[i] = signal_rate_global;
    }else{
      signal_rate_bins[i] = signal_rate_bins_gauss[i];
    }
    for(j in 1:nBins_1) {
      n_counts_recon_envelope[i][j] = exp(log_n_counts_recon_envelope[i][j]);
      n_counts_recon[i][j] = signal_rate_bins[i]*x_sig_1[i]*y_sig_1[j] +
	      		       background_rate_1*x_back_1[i]*y_back_1[j]+                              background_rate_2*x_back_2[i]*y_back_2[j];
      if (n_counts_recon[i,j]<0.){
       	 print(n_counts_recon[i,j]," Bin ", i,"  ", j);
      }
      if(total_bin_variation!=0.0){
        n_counts_recon[i][j] = n_counts_recon[i][j]*n_counts_recon_envelope[i][j];
      }
    }
  }
  signal_rate_smeared = exp(log_signal_rate_bins_gauss[1]);
}


model {
  if(sig_x_gauss_global_frac!=0.0){
    signal_rate_global_gauss ~ normal(signal_rate,sig_x_gauss_global_frac*signal_rate);
  }
  if(sig_x_gauss_bin_frac!=0.0){
    increment_log_prob(sum(log_signal_rate_bins_gauss));  // Jacobian
    signal_rate_bins_gauss ~ normal(signal_rate_global,sig_x_gauss_bin_frac*signal_rate_global);
  }
  if(total_bin_variation!=0.0){
    for(i in 1:nBins_0){
      increment_log_prob(sum(log_n_counts_recon_envelope[i])); // Jacobian
      n_counts_recon_envelope[i] ~ normal(1.0,total_bin_variation);
    }
  }

  // Fit data to poisson distribution of inputted generated data points
  for(i in 1:nBins_0){
    for(j in 1:nBins_1){
      target += poisson_lpmf(fake_data[i,j] | n_counts_recon[i,j]);
    }
  }
}
