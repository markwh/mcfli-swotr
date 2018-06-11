
data {
  
  // Dimensions
  int<lower=0> nx; // number of cross-sections
  int<lower=0> nt; // number of observation times
  
  
  // *Actual* data
  
  vector[nt] dAobs[nx]; // measured area difference from base area
  real<lower=0> dA_shift[nx]; // median(dA) - min(dA) for each location
  vector[nt] xobs[nx];
  
  // *Known* likelihood parameters
  vector<lower=0>[nt] sigma_err[nx]; // This is now a hyperparameter
  
  
  // Hyperparameters
  // vector[nt] logQ_hat;
  real z_hat;
  real logA0_hat[nx];
  
  real<lower=0> muz_sd; # error in predicting mean z
  real<lower=0> logA0_sd;
}

transformed data {
  vector[nt] dA_pos[nx];
  vector[nt] logx[nx];

  for (i in 1:nx) {
    dA_pos[i] = dAobs[i] - min(dAobs[i]); // make all dA positive
    logx[i] = log(xobs[i]);
  }
}

parameters {
  vector[nt] z;
  real<lower=0> sigma_z;
  
  real<lower=0> A0[nx];

  real mu_z;
}

transformed parameters {
  vector[nt] lhs[nx];
  vector[nt] logA[nx]; // log area for Manning's equation
  real A0_med[nx];
  
  for (i in 1:nx) {
    A0_med[i] = A0[i] + dA_shift[i];
    for (t in 1:nt) {
      logA[i, t] = log(A0[i] + dA_pos[i, t]);
    }
    
    lhs[i] = ((5. / 3. * logA[i]) + logx[i] - z) ./ sigma_err[i];
  }
  
  // print(logA[1,1])
  // print(logA[1,1] + logx[1,1] - z[1])
  // print(dA_pos[1,1])
  // print(logx[1, 1])
  // print(sigma_err[1,1])
  // print(lhs[1])
  // 
}

model {
  // Likelihood and observation error
  for (i in 1:nx) {
    lhs[i] ~ normal(0, 1); //already scaled by sigma_man
    
    A0_med[i] ~ lognormal(logA0_hat, logA0_sd);
    target += -(log(A0[i] + dA_pos[i]));
  }
  
  
  // Priors
  // logQ ~ normal(logQ_hat, logQ_sd);
  z ~ normal(mu_z, sigma_z);
  mu_z ~ normal(z_hat, muz_sd);
  sigma_z ~ normal(0, 1);

}
