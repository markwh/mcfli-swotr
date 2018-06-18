
data {
  
  // Dimensions
  int<lower=0> ns;
  int<lower=0> nt;
  
  // *Actual* data
  vector[nt] x[ns];
  vector<lower=0>[nt] dA[ns];
  real<lower=1> dA_shift[ns]; // median(dA) - min(dA) for each location

  // *Known* likelihood parameters
  real<lower=0> sigma_err;
  
  // Hyperparameters
  real mu_hat;
  real<lower=0> mu_sd; 
  real logA0_hat[ns];
  real<lower=0> logA0_sd;
}

transformed data {
  vector[ns] xbari;
  vector[ns] xbar_dev;
  
  for (i in 1:ns) {
    xbari[i] = mean(x[i]);
  }
  xbar_dev = xbari - mean(xbari);
}


parameters {
  vector[nt] y;
  real<lower=0> sigma_y;
  real mu;
  
  real
  
  real<lower=0> logA0bar;
  real<lower=0> logA0_dev[ns];
  
  real<lower=0> truesigma_err;
  
}

transformed parameters {
  
  vector[nt] z[ns];
  real logA0_med[ns];
  
  for (i in 1:ns) {
    logA0_med[i] = log(A0 + dA_shift);
    
    A0[i] = A0bar + A0dev;
    a[i] = (5. / 3. * log(A0[i] + dA[i]))
    z[i] = y - a[i];
  }
}

model {
  // Likelihood
  for (i in 1:ns) {
    x[i] ~ normal(z[i], truesigma_err); //already scaled by sigma_err
    
    logA0_dev[i] ~ normal(xbar_dev[i])
    
    // prior on A0
    logA0bar ~ normal(mean(logA0_hat), logA0_sd);
    // A0[i] + dA_shift[i] ~ lognormal(logA0_hat[i], logA0_sd);
  }
  
  // Priors
  truesigma_err ~ normal(0, sigma_err);
  y ~ normal(mu, sigma_y);
  sigma_y ~ normal(0.96, 0.4);
  mu ~ normal(mu_hat, mu_sd);
}
