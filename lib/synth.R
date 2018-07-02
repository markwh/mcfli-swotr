# Synthetic data

synth_new <- function(nx = 20, nt = 100, mu = runif(1, 0, 10),
                      noise_ar = "none", sigma_err = runif(1, 0, 0.3)) {
  
  sigma_z <- runif(1, 0.5, 1.5)
  
  if (noise_ar == "none") {
    z <- rnorm(nt, mu, sigma_z)
  } else {
    z <- arima.sim(model = list(ar = noise_ar), n = nt)
    z <- (z - mean(z)) * sigma_z / sd(z) + mu
  }
  
  errmat <- matrix(rnorm(nx * nt, 0, sigma_err), nrow = nx)
  lhs <- swot_vec2mat(z, errmat) + errmat
  
  mu_logA53 <- runif(nx, 5, 8)
  logA53 <- lhs + matrix(rnorm(nx * nt, rep(mu_logA53, nt), sigma_z / 2),
                         nrow = nx, byrow = FALSE)
  logx <- lhs - logA53
  
  neglogW23 <- logx / 2
  logS12 <- logx - neglogW23
  W <- exp(-3/2 * neglogW23)
  S <- exp(2 * logS12)
  
  mu_logA <- mu_logA53 * 3 / 5
  logA <- 3 / 5 * logA53
  
  A0vec <- apply(exp(logA), 1, min)
  dA <- exp(logA) %>% - swot_vec2mat(A0vec, logA)
  dA_shift <- apply(dA, 1, function(x) median(x) - min(x))
  
  out <- list(
    nx = nx,
    nt = nt,
    dAobs = dA,
    dA_shift = dA_shift,
    Wobs = W,
    Sobs = S,
    sigma_man = swot_vec2mat(rep(sigma_err, nx), dA),
    logQ_hat = 5,
    logQ_sd = 2,
    logn_hat = -3.5,
    logn_sd = 0.25,
    logA0_hat = rnorm(nx, log(A0vec), 0.25),
    logA0_sd = 0.25,
    lowerbound_logQ = 1,
    upperbound_logQ = 20,
    lowerbound_logn = -10,
    upperbound_logn = 10,
    lowerbound_A0 = 0,
    upperbound_A0 = 1e5
  )
  
  out
}