// Simplified version using matrix multiplication
// Only error is in dA measurement.
data {

  // Dimensions
  int<lower=2> nrow; // number of rows of model matrix
  int<lower=2> nx; // number of locations
  
  vector[nrow] rhs;
  matrix[nrow, nx] modmat;
  
  real logA0_hat;
  real sigma_logA0;
  real sigma_err;
}

parameters {
  vector<lower=0>[nx] A0;
  real<lower=0> truesigma_err;
}


model {
  modmat * A0 ~  normal(rhs, truesigma_err);

  A0 ~ lognormal(logA0_hat, sigma_logA0);
  truesigma_err ~ normal(0, sigma_err);
  
}
