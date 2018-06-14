
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

parameters {
  vector[nt] y;
  real<lower=0> sigma_y;
  real mu;
  real<lower=0> A0[ns];
  real<lower=0> truesigma_err;
}

transformed parameters {
  vector[nt] z[ns];

  for (i in 1:ns) {
    z[i] = y - (5. / 3. * log(A0[i] + dA[i]));
  }
}

model {
  // Likelihood
  for (i in 1:ns) {
    x[i] ~ normal(z[i], truesigma_err); //already scaled by sigma_err
    // x[i] ~ normal(z[i], 0.25); //already scaled by sigma_err
        
    // prior on A0
    A0[i] + dA_shift[i] ~ lognormal(logA0_hat, logA0_sd);
  }
  
  // Priors
  truesigma_err ~ normal(0, sigma_err);
  y ~ normal(mu, sigma_y);
  sigma_y ~ normal(0.96, 0.4);
  mu ~ normal(mu_hat, mu_sd);
}
