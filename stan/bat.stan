


functions {
  // Conversion from array to vector takes place by row.
  // Nested elements are appended in sequence. 
  // That is, a N (rows) by M (columns) matrix X becomes 
  // a vector [x11, x12, ..., x1M, x21, ..., xNM].

  
  // Convert an array to a vector based on a binary matrix 
  // indicating non-missing data
  
  vector ragged_vec(vector[] x, int[,] bin) {
    vector[num_elements(x)] out;
    int ind;
    
    ind = 1;
    for (i in 1:size(x)) {
      for (t in 1:num_elements(x[1])) {
        if (bin[i, t] == 1) {
          out[ind] = x[i, t];
          ind += 1;
        }
      }
    }
    // print(out);
    return(out[1:(ind - 1)]);
  }
  
  // Repeat elements of a "row" vector to match with 2-D array vectorization
  vector ragged_row(vector x, int[,] bin) {
    vector[num_elements(bin)] out;
    int ind;
    
    ind = 0;
    for (i in 1:size(bin)) {
      for (t in 1:num_elements(bin[1])) {
        if (bin[i, t] == 1) {
          ind += 1;
          out[ind] = x[t];
        }
      }
    }
    return(out[1:ind]);
  }
  
  // Repeat elements of a "column" vector to match with 2-D array vectorization
  vector ragged_col(vector x, int[,] bin) {
    vector[num_elements(bin)] out;
    int ind;
    
    ind = 0;
    for (i in 1:size(bin)) {
      for (t in 1:num_elements(bin[1])) {
        if (bin[i, t] == 1) {
          ind += 1;
          out[ind] = x[i];
        }
      }
    }
    return(out[1:ind]);
  }
  
  // Mean along a given index of an already-vectorized array. 
  vector ragged_mean(vector x, int[,] indarray, int[] numinds) {
    
    vector[num_elements(x)] meanvec;
    real meant;
    
    // print(num_elements(numinds))

    // initialize counters
    for (t in 1:num_elements(numinds)) {
      // print(indarray[t][1:numinds[t]])
      meant = mean(x[indarray[t][1:numinds[t]]]);
      for (s in 1:numinds[t]) {
        meanvec[indarray[t][s]] = meant;
      }
      
    }
    
    return(meanvec);
  }
  
  // indices of vectorized bin1 that are also in vectorized bin2
  int[] commoninds(int[,] bin1, int[,] bin2) {
    int out[num_elements(bin1)];
    int vecinds[size(bin1), num_elements(bin1[1])];
    int ctr;
    int ind;

    ctr = 0;
    for (i in 1:size(bin1)) {
      for (t in 1:num_elements(bin1[1])) {
        if (bin1[i, t] == 1) {
          ctr += 1;
          vecinds[i, t] = ctr;
        }
      }
    }

    ind = 0;
    for (i in 1:size(vecinds)) {
      for (t in 1:size(vecinds[1])) {
        if (bin2[i, t] == 1) {
          ind += 1;
          out[ind] = vecinds[i, t];
        }
      }
    }
    return(out[1:ind]);
  }
}

data {
  
  // Options
  int<lower=0, upper=1> meas_err; // 0=no; 1=yes;
  
  // Dimensions
  int<lower=0> nx; // number of reaches
  int<lower=0> nt; // number of times
  int<lower=0> ntot_man; // total number of non-missing Manning observations

  // Missing data
  int<lower=0,upper=1> hasdat[nx, nt]; // matrix of 0 (missing), 1 (not missing)
  int<lower=0> indmat[nt, nx];
  int<lower=1> nindvec[nt];

  // *Actual* data
  vector[nt] Wobs[nx]; // measured widths, including placeholders for missing
  vector[nt] Sobs[nx]; // measured slopes
  vector[nt] dAobs[nx]; // measured partial area
  vector[nx] dA_shift; // adjustment from min to median

  real<lower=0> Werr_sd;
  real<lower=0> Serr_sd;
  real<lower=0> dAerr_sd;
 
 
  // Hard bounds on parameters
  real lowerbound_A0; // These must be scalars, unfortunately. 
  real upperbound_A0;

  // *Known* likelihood parameters
  vector<lower=0>[nt] sigma_man[nx]; // Manning error standard deviation

  // Hyperparameters
  real logA0_hat[nx];
  real<lower=0> logA0_sd;
}



transformed data {
  // Transformed data are *vectors*, not arrays. This to allow ragged structure

  vector[nt] dApos_array[nx];
  
  vector[ntot_man] Wobsvec;
  vector[ntot_man] Sobs_vec;
  
  vector[ntot_man] logWobs;
  vector[ntot_man] logSobs;
  vector[ntot_man] dApos_obs;
  vector[nt] sigma_man_adj[nx]; // Adjusted Manning stdev, accounting for subtraction of mean
  vector[ntot_man] sigmavec_man;

  for (i in 1:nx) {
    dApos_array[i] = dAobs[i] - min(dAobs[i]); // make all dA positive
  }
  
  // convert pseudo-ragged arrays to vectors
  Wobsvec = ragged_vec(Wobs, hasdat);
  Sobs_vec = ragged_vec(Sobs, hasdat);
  dApos_obs = ragged_vec(dApos_array, hasdat);
  
  logWobs = log(Wobsvec);
  logSobs = log(Sobs_vec);
  
  // Adjust standard deviation to account for averaging operation. 
  for (i in 1:nx) {
    for (j in 1:nt) {
      sigma_man_adj[i][j] = sigma_man[i,j] * (nindvec[j] + 1.0) / (nindvec[j]);
    }
  }
  
  sigmavec_man = ragged_vec(sigma_man_adj, hasdat);
  // print(sigmavec_man)
}

parameters {
  vector<lower=lowerbound_A0,upper=upperbound_A0>[nx] A0;
  
  vector<lower=0>[ntot_man] Wact[meas_err];
  vector<lower=0>[ntot_man] Sact[meas_err];
  vector[ntot_man] dApos_act[meas_err];
}


transformed parameters {

  vector[ntot_man] logA; // log area for Manning's equation
  vector[ntot_man] logA_ctr;
  vector[ntot_man] man_lhs; // centered slope-width transformed
  vector[ntot_man] man_lhs_ctr; // Mean across space of man_lhs
  

  // Manning params
  if (meas_err) {
    logA = 10. * log(ragged_col(A0, hasdat) + dApos_act[1]);
    man_lhs = 4. * log(Wact[1]) - 3. * log(Sact[1]);
  }
  else {
    logA = 10. * log(ragged_col(A0, hasdat) + dApos_obs);
    man_lhs = 4. * logWobs - 3. * logSobs;
  }
  
  logA_ctr = logA - ragged_mean(logA, indmat, nindvec);
  man_lhs_ctr = man_lhs - ragged_mean(man_lhs, indmat, nindvec);
}

model {
  
  // Prior on A0
  A0 + dA_shift ~ lognormal(logA0_hat, logA0_sd);

  // Likelihood and observation error
    
  // Manning likelihood
  man_lhs_ctr ~ normal(logA_ctr, 6. * sigmavec_man);

  // Latent vars for measurement error
  if (meas_err) {
    Wact[1] ~ normal(Wobsvec, Werr_sd); // W meas err
    Sact[1] ~ normal(Sobs_vec, Serr_sd); // S meas err
    dApos_act[1] ~ normal(dApos_obs, dAerr_sd); // dA meas err
    target += -log(Wact[1]); // Jacobian adjustments
    target += -log(Sact[1]);

  }
}
