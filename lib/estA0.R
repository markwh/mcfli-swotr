# A0 estimation using Durand 2010 method

# These functions from 20180321. 
make_omega <- function(n, p1, p2) {
  stopifnot(p1 < p2 && p2 <= n && p1 >= 1)
  out <- rep(0L, n)
  out[p1] <- 1L
  out[p2] <- -1L
  out
}

make_omegas <- function(n) {
  p1s <- as.list(1:(n - 1))
  p2s <- map(p1s, function(x) (x + 1):n)
  out <- map2(p1s, p2s, function(x, y) map(y, ~make_omega(n, x, .)))
  unlist(out, recursive = FALSE)
}


estA0 <- function(wmat, smat, damat, zero = c("first", "minimum", "median"),
                  random_omega = FALSE,
                  wexp = -2/3, sexp = 1/2, aexp = 5/3) {
  
  stopifnot(dim(wmat) == dim(smat) && dim(smat) == dim(damat))
  
  zero <- match.arg(zero)
  
  # Assume these inputs come as DAWG-style time-across, space-down matrices. 
  # Transpose.
  wmat <- t(wmat)
  smat <- t(smat)
  damat <- t(damat)
  xmat <- (wmat ^ (wexp) * smat ^ (sexp)) ^ (1 / aexp)
  
  nx <- ncol(xmat)
  
  
  if (random_omega) {
    omegas <- lapply(1:100, function(x) rnorm(nx))
    omegas <- lapply(omegas, function(x) x - mean(x))
  } else {
    omegas <- make_omegas(nx)
  }

  
  modmats <- lapply(omegas, function(x) xmat %*% diag(x))
  
  modmat <- Reduce(rbind, modmats)
  
  rsps <- lapply(omegas, function(x) (xmat * -damat) %*% x) 
  rsp <- unlist(rsps)
  
  lmdf <- setNames(as.data.frame(cbind(rsp, modmat)), c("y", paste0("x", 1:nx)))
  
  lmout <- lm(y ~ 0 + ., data = lmdf, y = TRUE)
  # out <- coef(lmout)
  # out
  lmout
}


# Modified from notebook20180312.Rmd
estA0_list <- function(obslist, keeplocs = "all", keeptimes = "all") {
  if (length(keeplocs) == 1 && keeplocs == "all")
    keeplocs <- 1:nrow(obslist$W)
  if (length(keeptimes) == 1 && keeptimes == "all")
    keeptimes <- 1:nrow(obslist$W)
  
  with(obslist, estA0(wmat = W[keeplocs, ], 
                        smat = S[keeplocs, ], 
                        damat = dA[keeplocs, ]))
}


# Fix calcdA_mat function to give different zero reference points. 

calcdA_mat <- function (w, h, zero = c("first", "minimum", "median")) {
  zero <- match.arg(zero)
  stopifnot(all(dim(w) == dim(h)))
  dA <- w
  for (i in 1:nrow(dA)) {
    dA[i, ] <- calcdA_vec(w[i, ], h[i, ], zero = zero)
  }
  dA
}

calcdA_vec <- function(w, h, zero = c("first", "minimum", "median")) {
  zero <- match.arg(zero)
  words <- order(w)
  warr <- w[words]
  harr <- h[words]
  delh <- c(0, diff(harr))
  delA <- cumsum(warr * delh)
  dA <- 1:length(w)
  dA[words] <- delA
  
  if (zero == "minimum")
    dA <- dA - min(dA, na.rm = TRUE)
  
  if (zero == "median") 
    dA <- dA - median(dA, na.rm = TRUE)
  
  dA
}