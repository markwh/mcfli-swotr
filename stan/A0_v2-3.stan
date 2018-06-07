// Simplified version using matrix multiplication
// Now adding in closure uncertainty as independent errors. 
// This version puts a prior on sigma_err.
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

  real<lower= 0> sigma_dA; // measurement error sd on dA
  // real<lower= 0> sigma_dgdx;
  // real<lower= 0> sigma_nu;
  // real<lower= 0> sigma_err;
}

transformed data {
  vector[nx] x_dev;
  matrix[nt, nx] Mobs;
  
  matrix[nx, nx] Omega[(nx * (nx - 1) / 2)];
  matrix[nt, nx] dA_pos;


  x_dev = x - mean(x);
  for (i in 1:nx) {
    for (t in 1:nt) {
      Mobs[t, i] = W[t, i] ^ (-2. / 5.) * S[t, i] ^ (3. / 10.);
    }
    dA_pos[, i] = dA[, i] - min(dA[, i]); // make all dA positive
  }
  
  for (i in 1:(nx * (nx - 1) / 2)) {
    Omega[i] = diag_matrix(omega[i]);
  }
}

parameters {
  vector<lower=0>[nx] A0;
  matrix<lower = 0>[nt, nx] M;
  real<lower=0> sigma_err;
  // vector[nt] dgdx;
  // real alpha[nx];
}

transformed parameters {
  matrix[nt, nx] modmats[(nx * (nx - 1) / 2)];
  matrix[nt * (nx * (nx - 1) / 2), nx] modmat;
  vector[nt] rhss[(nx * (nx - 1) / 2)];
  vector[nt * (nx * (nx - 1) / 2)] rhs;
  
  // initialize model matrix and response to be filled in via loop
  modmat = rep_matrix(0, nt * (nx * (nx - 1) / 2), nx);
  rhs = rep_vector(0, nt * (nx * (nx - 1) / 2));
  
  // fill in model matrix and response vector
  for (i in 1:(nx * (nx - 1) / 2)) {
    modmats[i] = M * Omega[i];
    rhss[i] = -(M .* dA_pos) * omega[i]; 
    
    for (t in 1:nt) {
      modmat[(i - 1) * nx + t, ] = modmats[i][t, ];
      rhs[(i - 1) * nx + t] = rhss[i][t];
    }
  }
}

model {
  
  for (i in 1:nx) {
    M[, i] ./ Mobs[, i] ~ lognormal(0, sigma_err);
  }

  // for (i in 1:(nx * (nx - 1) / 2)) {
    modmat * A0 ~  normal(rhs, 2 * sigma_dA);
  //   lhs[i] ~ normal(rhs[i], sigma_err);
  //   
  //   target += -logA[i];
  // }
  
  A0 ~ lognormal(logA0_hat, sigma_logA0);
  // dgdx ~ normal(0, sigma_dgdx);
  // alpha ~ normal(0, sigma_alpha);
  sigma_err ~ normal(0, 0.01);
}




