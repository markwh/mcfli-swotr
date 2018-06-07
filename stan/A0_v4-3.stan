// This version doesn't use matrix multiplication or pairwise equalitites
// Instead, it infers Qn. and centers LHS on that. 
// First, assume closure term is zero. Only error is from dA measurement.
// Now (v4-3) treat closure as broken into time-variant and space-variant parts. 
data {

  // Dimensions
  int<lower=2> nt; // number of observation times
  int<lower=2> nx; // number of locations
  vector[nt] dA[nx];
  vector[nt] W[nx];
  vector[nt] S[nx];
  vector[nx] x;

  // real<lower= 0> sigma_dA; // measurement error sd on dA
  real<lower=0> logA0_hat;
  real<lower=0> sigma_logA0;
  
  real logQ_hat;
  real sigma_logQ;
  real logn_hat;
  real sigma_logn;
  
}

transformed data {
  vector[nt] Mobs[nx];
  vector[nt] dA_pos[nx];
  vector[nx] x_dev;
  
  for (i in 1:nx) {
    for (t in 1:nt) {
      Mobs[i, t] = (-2. / 3.) * log(W[i, t]) + 0.5 * log(S[i, t]);
    }
    dA_pos[i] = dA[i] - min(dA[i]); // make all dA positive
  }
  
  x_dev = (x - mean(x)) / 1000.; // make this have units of km. 
}

parameters {
  vector<lower=0>[nx] A0;
  vector<lower=0>[nt] Qn;
  real<lower=0> sigma_logA;
  
  vector[nt] dgdx;
  vector[nx] alpha;
  
  real<lower=0> sigma_dgdx;
  real<lower=0> sigma_alpha;
}

transformed parameters {
  // vector<lower=0>[nt] A[nx];
  vector[nt] M[nx];
  vector[nt] lhs[nx];
  for (i in 1:nx) {
    // logA[i] = 5./3. * log(dA_pos[i] + A0[i]);
    lhs[i] = log(Qn) - 5./3. * log(dA_pos[i] + A0[i]);
    M[i] = Mobs[i] + dgdx * x_dev[i] + alpha[i];
  }
}

model {
  
  for (i in 1:nx) {
    lhs[i] ~ normal(M[i], sigma_logA);
    target += -log(A0[i] + dA_pos[i]);
    target += log(5. / 3.);
  }
  target += -log(Qn);
  
  A0 ~ lognormal(logA0_hat, sigma_logA0);
  Qn ~ lognormal(logQ_hat + logn_hat, sqrt(sigma_logQ^2 + sigma_logn^2));
  dgdx ~ normal(0, sigma_dgdx);
  alpha ~ normal(0, sigma_alpha);
  
  sigma_logA ~ normal(0, 1);
  
  sigma_dgdx ~ normal(0, 0.05);
  sigma_alpha ~ normal(0, 0.0001);

}




