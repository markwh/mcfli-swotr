# functions I use but that need refining.

# from 04/07 notebook

omegaProduct <- function(ws35mat, symmetric = FALSE) {
  omegas <- make_omegas(nrow(ws35mat), symmetric = symmetric)
  out <- map(omegas, function(x) t(ws35mat) %*% diag(x)) %>% 
    Reduce(rbind, .)
  out
}

omegaProduct_y <- function(ws35mat, dAmat, symmetric = FALSE) {
  omegas <- make_omegas(nrow(ws35mat), symmetric = symmetric)
  
  out <- map(omegas, function(x) t(-ws35mat * dAmat) %*% x) %>% 
    unlist()
  out
}


# Symmetric finite difference functions that preserve dimension
# from 0417
findif_x <- function(dawgmat) {
  difmat <- apply(dawgmat, 2, diff)
  out1 <- rbind(difmat[1, ], difmat)
  out2 <- rbind(difmat, difmat[nrow(difmat), ])
  out <- (out1 + out2) / 2
  out
}

findif_t <- function(dawgmat) {
  difmat <- t(apply(dawgmat, 1, diff))
  out1 <- cbind(difmat[, 1], difmat)
  out2 <- cbind(difmat, difmat[, ncol(difmat)])
  out <- (out1 + out2) / 2
  out
}

# data prep for Stan A0 inference
# from 0420

bat_data <- function(swotlist, logA0_hat = 7, sigma_logA0 = 7, 
                     logQ_hat = 7, sigma_logQ = 7, 
                     logn_hat = -5, sigma_logn = 10,
                     sigma_dgdx = 1e-5, sigma_alpha = 1e-1, sigma_err = 1e-1) {
  W <- swotlist$W
  dA <- swotlist$dA
  S <- swotlist$S
  x <- swotlist$x[, 1]
  
  out <- list(nt = ncol(W), nx = nrow(W), 
              dA = dA, W = W, S = S, x = x, 
              logA0_hat = logA0_hat, sigma_logA0 = sigma_logA0,
              logQ_hat = logQ_hat, sigma_logQ = sigma_logQ,
              logn_hat = logn_hat, sigma_logn = sigma_logn,
              sigma_dgdx = sigma_dgdx, sigma_alpha = sigma_alpha, 
              sigma_err = sigma_err)
  out
}
