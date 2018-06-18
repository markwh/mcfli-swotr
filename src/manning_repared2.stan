
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
  real<upper=0.1> mindA[ns]; // Just to verify that dA is appropriately shifted.
  vector[nt] dA_med[ns]; // median-shifted dA
  real c[ns];

  for (i in 1:ns) {
    mindA[i] = min(dA[i]);
    dA_med[i] = dA[i] + dA_shift[i];
    c[i] = - 1. / min(dA_med[i]);
  }
}

parameters {
  vector[nt] y;
  real<lower=0> sigma_y;
  real mu_y;
  
  real<lower=1> cA0[ns];
  
  real<lower=0> truesigma_err;
  real mu_ai[ns];
  real mu_a;
  
}

transformed parameters {
  vector[nt] z[ns];
  real A0_med[ns];
  
  for (i in 1:ns) {
    A0_med[i] = cA0[i] / c[i];
    z[i] = y - (5. / 3. * log(A0_med[i] + dA_med[i]));
  }
}

model {
  // Likelihood
  for (i in 1:ns) {
    x[i] ~ normal(z[i], truesigma_err); //already scaled by sigma_err
    // x[i] ~ normal(z[i], 0.25); //already scaled by sigma_err
        
    // prior on A0
    A0_med[i] ~ lognormal(mu_a, sd(logA0_hat));
  }
  
  // Priors
  truesigma_err ~ normal(0, sigma_err);
  y ~ normal(mu_y, sigma_y);
  sigma_y ~ normal(0.96, 0.4);
  mu_y ~ normal(mu_hat, mu_sd);
  mu_a ~ normal(mean(logA0_hat), logA0_sd);
}
