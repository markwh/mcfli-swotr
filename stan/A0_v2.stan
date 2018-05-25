// Simplified version using matrix multiplication
// Only error is in dA measurement.
data {

  // Dimensions
  int<lower=2> nt; // number of observation times
  int<lower=2> nx; // number of locations

  vector[nx] x;
  vector[nx * (nx - 1) / 2] omega[nx];

  matrix[nt, nx] dA;
  matrix[nt, nx] W;
  matrix[nt, nx] S;
  
  real logA0_hat;
  real sigma_logA0;

  // real logQ_hat;
  // real sigma_logQ;
  // 
  // real logn_hat;
  // real sigma_logn;

  // real sigma_dgdx;
  // real sigma_alpha;
  real sigma_err;
}

transformed data {

  vector[nx] x_dev;
  matrix[nt, nx] Mobs;
  
  matrix[nt, nx] modmats[(nx * (nx - 1) / 2)];
  matrix[nt * (nx * (nx - 1) / 2), nx] modmat;
  
  matrix[nx, nx] Omega[(nx * (nx - 1) / 2)];
  matrix[nt, nx] dA_pos;
  
  vector[nt] rhss[(nx * (nx - 1) / 2)];
  vector[nt * (nx * (nx - 1) / 2)] rhs;

  x_dev = x - mean(x);
  for (i in 1:nx) {
    for (t in 1:nt) {
      Mobs[t, i] = W[t, i] ^ (-2. / 5.) * S[t, i] ^ (3. / 10.);
    }
    dA_pos[, i] = dA[, i] - min(dA[, i]); // make all dA positive
  }
  
  modmat = rep_matrix(0, nt * (nx * (nx - 1) / 2), nx);
  rhs = rep_vector(0, nt * (nx * (nx - 1) / 2));
  for (i in 1:(nx * (nx - 1) / 2)) {
    Omega[i] = diag_matrix(omega[i]);
    modmats[i] = Mobs * Omega[i];
    rhss[i] = -(Mobs .* dA_pos) * omega[i]; 
    
    for (t in 1:nt) {
      modmat[(i - 1) * nx + t, ] = modmats[i][t, ];
      rhs[(i - 1) * nx + t] = rhss[i][t];
    }
  }
  
  print(modmat);
}

parameters {
  vector<lower=0>[nx] A0;
  // vector[nt] logQn;
  // vector[nt] dgdx;
  // real alpha[nx];
}

transformed parameters {
  // matrix[nt, nx] logA;
  // matrix[nt, nx] gamma;
  // matrix[nt, nx] lhs;
  // matrix[nt, nx] rhs;

  // vector[nt] logA[nx];
  // vector[nt] gamma[nx];
  // vector[nt] lhs[nx];
  // vector[nt] rhs[nx];
  
  // for (i in 1:nx) {
  //   for (t in 1:nt) {
  //     logA[i, t] = log(A0[i] + dA_pos[i, t]);
  //     gamma[i, t] = dgdx[t] * x_dev[i];
  //   }
  // }

}

model {
  // for (i in 1:(nx * (nx - 1) / 2)) {
    modmat * A0 ~  normal(rhs, sigma_err);
  //   lhs[i] ~ normal(rhs[i], sigma_err);
  //   
  //   target += -logA[i];
  // }
  
  A0 ~ lognormal(logA0_hat, sigma_logA0);
  // logQn ~ normal(logQ_hat + logn_hat, sqrt(sigma_logQ^2 + sigma_logn^2));
  // dgdx ~ normal(0, sigma_dgdx);
  // alpha ~ normal(0, sigma_alpha);
}




