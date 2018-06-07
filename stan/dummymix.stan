// Testing (possibly breaking) my sanity

data {
  vector[100] zobs;
  real xhat;
  real<lower=0> x_sd;
  real yhat;
  real<lower=0> y_sd;
  real<lower=0> z_sd;
}

parameters {
  real x;
  real y;
  // real z;
  // real<lower=0> z_sd;
}

transformed parameters {
  real z;
  z = x - y;
}

model {
  zobs ~ normal(z, z_sd);
  // z_sd ~ normal(0, 1);
  
  // z ~ normal(xhat + y, x_sd);
  x ~ normal(xhat, x_sd);
  y ~ normal(yhat, y_sd);
}

