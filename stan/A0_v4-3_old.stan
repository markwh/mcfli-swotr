// This version doesn't use matrix multiplication or pairwise equalitites
// Instead, it infers Qn. and centers LHS on that. 
// In v4-1 I assumed closure term is zero. Only error is from dA measurement.
// In v4-2 I treated closure as random noise. 
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
  vector<lower=0>[nt] Mobs[nx];
  vector[nt] dA_pos[nx];
  vector[nx] x_dev;
  for (i in 1:nx) {
    for (t in 1:nt) {
      Mobs[i, t] = W[i, t] ^ (-2. / 5.) * S[i, t] ^ (3. / 10.);
    }
    dA_pos[i] = dA[i] - min(dA[i]); // make all dA positive
  }
  
  x_dev = (x - mean(x)) / 1000.; // Make this have units of km, for numerical ease. 
}

parameters {
  vector<lower=0>[nx] A0;
  vector<lower=0>[nt] Qn;
  real<lower=0> sigma_dA;
  
  // vector<lower=0>[nt] M[nx];
  vector[nt] dgdx;
  vector[nx] alpha;
  
  real<lower=0> sigma_dgdx;
  real<lower=0> sigma_alpha;
  // real<lower=0> sigma_err;
}

transformed parameters {
  vector<lower=0>[nt] A[nx];
  vector<lower=0>[nt] M[nx];
  
  for (i in 1:nx) {
    A[i] = dA_pos[i] + A0[i];
    M[i] = Mobs[i] .* exp(dgdx * x_dev[i] + alpha[i]);
  }
}

model {
  // print(sigma_dgdx);
  // print(dgdx);
  for (i in 1:nx) {
    A[i] ~ normal(Qn ./ M[i], sigma_dA);
    
    // M[i] ./ Mobs[i] ~ lognormal(0, sigma_err);
  }
  
  A0 ~ lognormal(logA0_hat, sigma_logA0);
  Qn ~ lognormal(logQ_hat + logn_hat, sqrt(sigma_logQ^2 + sigma_logn^2));
  dgdx ~ normal(0, sigma_dgdx);
  alpha ~ normal(0, sigma_alpha);
  
  sigma_dA ~ normal(0, 1.);
  // sigma_err ~ normal(0, 0.1);
  sigma_dgdx ~ normal(0, 0.05); // When multiplied by a characteristic reach distance (km)
                                // scale term will yield a characteristic flow change as fraction.
  sigma_alpha ~ normal(0, 0.1);
}




