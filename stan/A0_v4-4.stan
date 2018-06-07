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
  real sigma_logS; 
  real sigma_logW;
  
}

transformed data {
  // vector<lower=0>[nt] Mobs[nx];
  vector[nt] Mobs[nx];
  vector[nt] dA_pos[nx];
  
  for (i in 1:nx) {
    for (t in 1:nt) {
      Mobs[i, t] = (-2. / 3.) * log(W[i, t]) + 0.5 * log(S[i, t]);
    }
    dA_pos[i] = dA[i] - min(dA[i]); // make all dA positive
  }
  // print(Mobs);
}

parameters {
  vector[nt] Mtrue[nx];
  vector<lower=0>[nx] A0;
  vector<lower=0>[nt] Qn;
  // real<lower=0> sigma_dA;
  real<lower=0> sigma_logA;
}

transformed parameters {
  // vector<lower=0>[nt] A[nx];
  // vector[nt] logA[nx];
  vector[nt] lhs[nx];
  for (i in 1:nx) {
    // logA[i] = 5./3. * log(dA_pos[i] + A0[i]);
    lhs[i] = log(Qn) - 5./3. * log(dA_pos[i] + A0[i]);
  }
}

model {
  
  for (i in 1:nx) {
    // A[i] ~ lognormal((log(Qn) - log(Mobs[i]), sigma_dA);
    // logA[i] ~ normal(log(Qn) - Mobs[i], sigma_logA);
    lhs[i] ~ normal(Mtrue[i], sigma_logA);
    target += -log(A0[i] + dA_pos[i]);
    // target += -log(Qn);
    target += log(5. / 3.);
    
    Mtrue[i] ~ normal(Mobs[i], sqrt(4./9. * sigma_logW^2 + 0.25 * sigma_logS^2));
    
  }
  target += -log(Qn);
  
  A0 ~ lognormal(logA0_hat, sigma_logA0);
  Qn ~ lognormal(logQ_hat + logn_hat, sqrt(sigma_logQ^2 + sigma_logn^2));
  // sigma_dA ~ normal(0, 1);
  sigma_logA ~ normal(0, 1);
}




