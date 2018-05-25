// This version doesn't use matrix multiplication or pairwise equalitites
// Instead, it infers Qn. and centers LHS on that. 
// First, assume closure term is zero. Only error is from dA measurement.
data {

  // Dimensions
  int<lower=2> nt; // number of observation times
  int<lower=2> nx; // number of locations
  vector[nt] dA[nx];
  vector[nt] W[nx];
  vector[nt] S[nx];

  // real<lower= 0> sigma_dA; // measurement error sd on dA
  real<lower=0> logA0_hat;
  real<lower=0> sigma_logA0;
  
  real logQ_hat;
  real sigma_logQ;
  real logn_hat;
  real sigma_logn;
  
}

transformed data {
  vector<lower=0>[nt] Mobs[nx];
  vector[nt] dA_pos[nx];
  
  for (i in 1:nx) {
    for (t in 1:nt) {
      Mobs[i, t] = W[i, t] ^ (-2. / 5.) * S[i, t] ^ (3. / 10.);
    }
    dA_pos[i] = dA[i] - min(dA[i]); // make all dA positive
  }
  print(Mobs);
}

parameters {
  vector<lower=0>[nx] A0;
  vector<lower=0>[nt] Qn;
  real<lower=0> sigma_dA;
}

transformed parameters {
  vector<lower=0>[nt] A[nx];

  for (i in 1:nx) {
    A[i] = dA_pos[i] + A0[i];
  }
}

model {
  
  for (i in 1:nx) {
    A[i] ~ normal(Qn ./ Mobs[i], sigma_dA);
  }
  
  A0 ~ lognormal(logA0_hat, sigma_logA0);
  Qn ~ lognormal(logQ_hat + logn_hat, sqrt(sigma_logQ^2 + sigma_logn^2));
  sigma_dA ~ normal(0, 1);

}




