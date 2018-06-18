
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
  real logA0hat_bar;
  real logA0hat_sd;
  
  logA0hat_bar = mean(logA0_hat);
  logA0hat_sd = sd(logA0_hat);
}

parameters {
  vector[nt] y;
  real<lower=0> sigma_y;
  real mu;
  real<lower=0> A0[ns];
  real<lower=0> truesigma_err;
  real mu_logA0;
}

transformed parameters {
  vector[nt] z[ns];
  real A0_med[ns];
  
  for (i in 1:ns) {
    z[i] = y - (5. / 3. * log(A0[i] + dA[i]));
    A0_med[i] = A0[i] + dA_shift[i];
  }
}

model {
  // Likelihood
  for (i in 1:ns) {
    x[i] ~ normal(z[i], truesigma_err); //already scaled by sigma_err
    // x[i] ~ normal(z[i], 0.25); //already scaled by sigma_err
        
    // prior on A0
    A0_med[i] ~ lognormal(logA0_hat[i], logA0_sd);
    A0_med[i] ~ lognormal(mu_logA0, sigma_logA0);
  }
  
  // Priors
  truesigma_err ~ normal(0, sigma_err);
  y ~ normal(mu, sigma_y);
  sigma_y ~ normal(0.96, 0.4);
  mu ~ normal(mu_hat, mu_sd);
}
