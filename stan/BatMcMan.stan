
data {

  // Dimensions
  int<lower=2> nt; // number of observation times
  int<lower=2> nx; // number of locations
  vector[nt] dA[nx];
  vector[nt] W[nx];
  vector[nt] S[nx];
  real logA0_hat;
  real logA0_sd;
  real<lower=0> sigmalogA_hat;
  real<lower=0> sigmalogA_sd;
}

transformed data {
  real<lower = 0> minA0;
  
  minA0 = 1 - min(dA);
}

parameters {
  vector[nt] manningN[nx];
  real<lower=0> A0[nx];
}

transformed parameters {
  vector[nt] A[nx];
  // real logA0;
  
  A = (A0 + dA);
  logA0 = log(A0);
}



model {
  
  logA0 ~ normal(logA0_hat, logA0_sd);
  sigma_logA ~ normal(sigmalogA_hat, sigmalogA_sd);
  logA ~ normal(log(A0), sigma_logA);
  
  target += -logA;
  target += -logA0;

}
