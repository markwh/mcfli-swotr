
// A version that does not parse out the closure term. 
data {

  // Dimensions
  int<lower=2> nt; // number of observation times
  int<lower=2> nx; // number of locations
  vector[nt] dA[nx];
  vector[nt] W[nx];
  vector[nt] S[nx];
  vector[nx] x;
  
  real<lower= 0> sigma_dA; // measurement error sd on dA
  // real<lower= 0> sigma_dgdx;
  // real<lower= 0> sigma_nu;
  real<lower= 0> sigma_err;
  
  real<lower=0> logA0_hat;
  real<lower=0> logA0_sd;
  
}

transformed data {
  vector[nx] x_dev;
  vector<lower=0>[nt] Mobs[nx];
  vector[nt] dA_pos[nx];
  
  for (i in 1:nx) {
    for (t in 1:nt) {
      Mobs[i, t] = W[i, t] ^ (-2. / 5.) * S[i, t] ^ (3. / 10.);
    }
    dA_pos[i] = dA[i] - min(dA[i]); // make all dA positive
  }
  
  x_dev = x - mean(x);
}

parameters {
  vector<lower=10>[nx] A0;
  // vector<lower=-1. / max(x),upper=1./max(x)>[nt] dgdx;
  // real nu[nx];
  // vector[nt] epsilon[nx];
  // vector<lower=0>[nt] M[nx];
}

transformed parameters {
  vector<lower=10>[nt] A[nx];
  // vector<lower=0>[nt] M[nx];
  // vector[nt] gamma[nx];
  
  for (i in 1:nx) {
    // gamma[i] = dgdx * x_dev[i];
    // M[i] = Mobs[i] .* exp(gamma[i] + nu[i] + epsilon[i]);
    // M[i] = Mobs[i] .* exp(epsilon[i]);
    A[i] = dA_pos[i] + A0[i];
  }

}

model {
  
  for (i in 1:(nx - 1)) {
    for (j in (i + 1):nx) {
      for (t in 1:nt) {
        // print("M", M[i, t]);
        print("A", A[i, t])
        // print(exp(gamma[i, t] + nu[i] + epsilon[i, t]))
        // print(gamma[i, t]);
        
        // M[i, t] * A[i] - M[j, t] * A[j] ~ normal(0, 2. * sigma_dA);
        Mobs[i, t] * A[i] - Mobs[j, t] * A[j] ~ normal(0, 2. * sigma_dA);
        // target += -log(M[i, t]);
        // target += -log(M[j, t]);
      }
    }
  }
  
  A0 ~ lognormal(logA0_hat, logA0_sd);
  // for (i in 1:nx) {
  //   M[i] ~ lognormal(log(Mobs[i]), sigma_err);
  //   // epsilon[i] ~ normal(0, sigma_err);
  // }
  // dgdx ~ normal(0, sigma_dgdx);
  // nu ~ normal(0, sigma_nu);
}




