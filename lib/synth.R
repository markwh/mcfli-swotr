# Synthetic data


#' Produces a new swotlist, which can then be put into a bamr or pared model. 
synth_new <- function(nx = 20, nt = 100, mu_q = runif(1, 0, 10),
                      sigma_q = runif(1, 0.5, 1.5),
                      noise_ar = "none", sigma_err = runif(1, 0, 0.3)) {
  
  if (noise_ar == "none") {
    logq <- rnorm(nt, mu_q, sigma_q)
  } else {
    logq <- arima.sim(model = list(ar = noise_ar), n = nt)
    logq <- (logq - mean(logq)) * sigma_q / sd(logq) + mu_q
  }
  
  
  errmat <- matrix(rnorm(nx * nt, 0, sigma_err), nrow = nx)
  lhs <- swot_vec2mat(logq, errmat) + errmat
  
  mu_logA53 <- runif(nx, 5, 8)
  
  # logA should be strongly correlated with logQ.
  
  logA53 <- lhs + matrix(rnorm(nx * nt, rep(mu_logA53, nt), sigma_q / 2),
                         nrow = nx, byrow = FALSE)
  logx <- lhs - logA53
  
  neglogW23 <- logx / 2
  logS12 <- logx - neglogW23
  W <- exp(-3/2 * neglogW23)
  S <- exp(2 * logS12)
  
  Q <- swot_vec2mat(exp(logq), W)
  
  mu_logA <- mu_logA53 * 3 / 5
  logA <- 3 / 5 * logA53
  
  A0vec <- apply(exp(logA), 1, min)
  dA <- exp(logA) %>% - swot_vec2mat(A0vec, logA)
  dA_shift <- apply(dA, 1, function(x) median(x) - min(x))
  
  out <- list(
    dA = dA,
    A = exp(logA),
    W = W,
    S = S,
    Q = Q
  )
  
  attr(out, "QWBM") <- rlnorm(1, mu_q, 1)
  
  out
}
