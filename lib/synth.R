# Synthetic data

#' Simulate y given x, where y has known mean, sd, and cor with x
sim_cor <- function(x, mu_y, sd_y, rho_xy) {
  x_scl <- as.numeric(scale(x))
  sig_er <- sqrt(rho_xy^(-2) - 1)
  y_scl <- x_scl + rnorm(length(x), 0, sig_er)
  out <- y_scl * sd_y + mu_y
  out
}


#' Produces a new swotlist, which can then be put into a bamr or pared model. 
synth_new <- function(nr = 20, nt = 100, 
                      logq_ar1 = "none") {
  
  moms <- synth_moments(nr = nr)
  
  if (logq_ar1 == "none") {
    logq <- rnorm(nt, moms$mu_q, moms$sigma_q)
  } else {
    logq <- arima.sim(model = list(ar = logq_ar1), n = nt)
    logq <- as.numeric((logq - mean(logq)) * moms$sigma_q / sd(logq) + 
                         moms$mu_q)
  }
  
  errmat <- matrix(rnorm(nr * nt, 0, moms$sigma_err), nrow = nr)
  lhs <- (swot_vec2mat(logq, errmat) + errmat) %>% 
    t() %>% 
    as.data.frame()
  
  logA <- pmap(list(x = lhs, mu = moms$mu_a, 
                    sd = moms$sigma_a, rho = moms$rho_qa),
               function(x, mu, sd, rho) sim_cor(x = x, mu_y = mu, sd_y = sd, 
                                             rho_xy = rho))
  
  logW <- pmap(list(x = logA, mu = moms$mu_w, 
                    sd = moms$sigma_w, rho = moms$rho_aw),
               function(x, mu, sd, rho) sim_cor(x = x, mu_y = mu, 
                                                sd_y = sd, rho_xy = rho))
  
  logn <- rnorm(1, -3.5, 0.25)
  
  logS <- pmap(list(q = lhs, a = logA, w = logW),
               function(q, a, w) 2 * q  - 10 / 3 * a + 4 / 3 * w + 2 * logn)
  
  tomat <- function(lst) {
    reduce(map(lst, ~t(.)), rbind)
  }
  
  W <- exp(tomat(logW))
  S <- exp(tomat(logS))
  A <- exp(tomat(logA))
  
  Q <- swot_vec2mat(exp(logq), W)
  
  A0vec <- apply(A, 1, min)
  dA <- A - swot_vec2mat(A0vec, W)

  swotlist <- list(
    dA = dA,
    A = A,
    W = W,
    S = S,
    Q = Q
  )
  
  attr(swotlist, "QWBM") <- rlnorm(1, moms$mu_q, 1)
  
  params <- c(moms, list(logn = logn))
  
  out <- list(swotlist = swotlist, params = params)
  
  out
}


synth_moments <- function(nr = 1) {
  mu_q <- rnorm(1, 6, 1) # Made this larger
  mu_a <- 1.85 + 0.753 * mu_q + rnorm(nr, 0, 0.65)
  mu_w <- 1.31 + 0.587 * mu_a + rnorm(nr, 0, 0.295)
  
  sigma_q <- abs(rnorm(1, 0.8, 0.2)) # Tightened this from 0.4, and lowered the mean from 1.0
  sigma_a <- abs(0.02 + 0.50 * sigma_q + rnorm(nr, 0, 0.2)) # Tightened from 0.26
  sigma_w <- abs(-0.02 + 0.53 * sigma_a + rnorm(nr, 0, 0.13))
  sigma_err <- runif(1, 0, 0.3)
  
  rhotrans_qa <- rnorm(nr, -3, 0.7) # Fudged this one to make it behave
  rhotrans_aw <- rnorm(nr, -1.67, 0.7) # Fudged just a bit (tightened)
  
  rho_qa <- 1 - exp(rhotrans_qa)
  rho_aw <- 1 - exp(rhotrans_aw)
  
  rho_qa[rho_qa < 0] <- runif(sum(rho_qa < 0), 0, 1)
  rho_aw[rho_aw < 0] <- runif(sum(rho_aw < 0), 0, 1)
  
  out <- list(mu_q = mu_q, 
              mu_a = mu_a,
              mu_w = mu_w,
              
              sigma_q = sigma_q,
              sigma_a = sigma_a,
              sigma_w = sigma_w,
              sigma_err = sigma_err,
              
              rho_qa = rho_qa,
              rho_aw = rho_aw)
  
  out
}

