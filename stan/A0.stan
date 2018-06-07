
data {

  // Dimensions
  int<lower=2> nt; // number of observation times
  int<lower=2> nx; // number of locations
  vector[nt] dA[nx];
  vector[nt] W[nx];
  vector[nt] S[nx];
  vector[nx] x;
  
  // real dA[nt, nx];
  // real W[nt, nx];
  // real S[nt, nx];
  
  // matrix[nt, nx] dA;
  // matrix[nt, nx] W;
  // matrix[nt, nx] S;
  
  real logA0_hat;
  real sigma_logA0;

  real logQ_hat;
  real sigma_logQ;

  real logn_hat;
  real sigma_logn;

  real sigma_dgdx;
  real sigma_alpha;
  real sigma_err;
}

transformed data {
  vector[nt] logW[nx];
  vector[nt] logS[nx];
  vector[nt] logM[nx];
  vector[nx] x_dev;
  
  // real logW[nt, nx];
  // real logS[nt, nx];
  // real logM[nt, nx];
  
  // matrix[nt, nx] logW;
  // matrix[nt, nx] logS;
  // matrix[nt, nx] logM;
  
  vector[nt] dA_pos[nx];
  // vector dA_mins[nx];

  for (i in 1:nx) {
    dA_pos[i] = dA[i] - min(dA[i]); // make all dA positive
  }
  
  logW = log(W);
  logS = log(S);
  
  for (i in 1:nx) {
    logM[i] = -(2. / 5.) * logW[i] + (3. / 10.) * logS[i];
  }

  
  x_dev = x - mean(x);
}

parameters {
  real<lower=0> A0[nx];
  vector[nt] logQn;
  vector[nt] dgdx;
  real alpha[nx];
}

transformed parameters {
  // matrix[nt, nx] logA;
  // matrix[nt, nx] gamma;
  // matrix[nt, nx] lhs;
  // matrix[nt, nx] rhs;

  vector[nt] logA[nx];
  vector[nt] gamma[nx];
  vector[nt] lhs[nx];
  vector[nt] rhs[nx];
  
  for (i in 1:nx) {
    for (t in 1:nt) {
      logA[i, t] = log(A0[i] + dA_pos[i, t]);
      gamma[i, t] = dgdx[t] * x_dev[i];
    }
    
    lhs[i] = logA[i] + logM[i];
    rhs[i] = logQn[i] + gamma[i] + alpha[i];
  }

}

model {
  for (i in 1:nx) {
    lhs[i] ~ normal(rhs[i], sigma_err);
    
    target += -logA[i];
  }
  
  A0 ~ lognormal(logA0_hat, sigma_logA0);
  logQn ~ normal(logQ_hat + logn_hat, sqrt(sigma_logQ^2 + sigma_logn^2));
  dgdx ~ normal(0, sigma_dgdx);
  alpha ~ normal(0, sigma_alpha);
}




