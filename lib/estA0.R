# A0 estimation using Durand 2010 method

# These functions from 20180321. 
make_omega <- function(n, p1, p2) {
  stopifnot(p1 < p2 && p2 <= n && p1 >= 1)
  out <- rep(0L, n)
  out[p1] <- 1L
  out[p2] <- -1L
  out
}

make_omegas <- function(n, symmetric = FALSE) {
  p1s <- as.list(1:(n - 1))
  p2s <- map(p1s, function(x) (x + 1):n)
  out1 <- map2(p1s, p2s, function(x, y) map(y, ~make_omega(n, x, .)))
  out <- unlist(out1, recursive = FALSE)
  
  if (symmetric) {
    out <- c(out, lapply(out, `-`))
  }
  out
}


estA0_lm <- function(swotlist, zero = c("same", "first", "minimum", "median"),
                     weight = TRUE, 
                     ndot = 1, random_omega = FALSE, intercept = FALSE, 
                     symmetric = FALSE) {
  
  zero <- match.arg(zero)
  
  ws35 <- manning_ws35(swotlist = swotlist, ndot = ndot)
  dA <- swotlist$dA
  
  
  if (is.null(swotlist$x)) {
    message("No distance info present. Using unweighted regression")
    xvec <- NULL
    weight <- FALSE
  } else{
    xvec <- as.vector(scale(swotlist$x[, 1]))
  }

  lmout <- estA0_lm_ws35(ws35mat = ws35, dAmat = dA, xvec = xvec, 
                         weight = weight, zero = zero, 
                         random_omega = random_omega, intercept = intercept, 
                         symmetric = symmetric)
  
  lmout
}

estA0 <- function(swotlist, zero = c("same", "first", "minimum", "median"),
                  weight = TRUE,
                  ndot = 1, random_omega = FALSE, symmetric = FALSE, intercept = FALSE) {
  
  zero <- match.arg(zero)
  outlm <- estA0_lm(swotlist = swotlist, zero = zero, ndot = ndot, 
                    weight = weight,
                    random_omega = random_omega, symmetric = symmetric,
                    intercept = intercept)
  out <- coef(outlm)
  out
}


# Stepwise location selection ---------------------------------------------

estA0_prune_condno <- function(swotlist) {
  
  wmat <- swotlist$W
  smat <- swotlist$S
  amat <- swotlist$A
  dAmat <- swotlist$dA
  
  ws35 <- t(wmat ^ (-2/3) * smat ^ (1/2)) ^ (3/5)
  recres <- recondition(ws35)
  
  rmlocs <- Reduce(c, recres$elim[-1], accumulate = TRUE)
  keeplocs <- c(list(0), map(rmlocs, `-`))
  
  sl_ss <- map(keeplocs, ~swot_sset(swotlist = swotlist, keeplocs = .))
  A0lms <- map(sl_ss, ~estA0_lm(swotlist = ., zero = "minimum"))
  
  A0r2 <- map_dbl(A0lms, A0_R2)
  
  allocs <- 1:nrow(wmat)
  coeflocs <- c(list(allocs), map(rmlocs, ~allocs[-.]))
  A0ests <- map2(A0lms, coeflocs, function(x, y) setNames(coef(x), y)) %>% 
    map(~data.frame(loc = names(.), A0_est = .)) %>% 
    setNames(0:(length(.) - 1)) %>% 
    bind_rows(.id = "n_elim") %>% 
    mutate(n_elim = as.integer(n_elim), 
           loc = as.integer(loc))
  
  out <- left_join(A0ests, recres, by = "n_elim") %>% 
    left_join(A0r2, by = "n_elim") %>% 
    select(loc, condno, R2, elim, n_elim, A0_est)
  
  if (!is.null(amat)) {
    A0reals <- data.frame(loc = 1:nrow(wmat), A0_real = (amat - dAmat)[, 1])
    out <- out %>% 
      left_join(A0reals, by = "loc")
  }
  
  out
}

#' Iteratively estimate A0 and ndot.
estA0_ndot <- function(swotlist, iters = 10L, trace = FALSE, count = FALSE) {
  
  # initialize ndot, taking all n equal
  ndot <- 1
  
  qdots <- list()
  ndots <- list()
  qbarns <- list()
  A0hats <- list()
  
  for (i in 1:iters) {
    # estimate A0
    A0hat <- estA0(swotlist, ndot = ndot)
    minW <- apply(swotlist$W, 1, min)
    A0hat[A0hat < 0] <- minW
    swotlist$Ahat <- swotlist$dA + swot_vec2mat(A0hat, swotlist$dA)
    # estimate ndot
    ndot <- manning_ndot(swotlist, Avar = "Ahat", log = FALSE, meanspace = "linA")
    
    
    qbarn <- manning_qbarn(swotlist, Avar = "Ahat", log = FALSE)
    ndots[[i]] <- log(ndot)
    qdots[[i]] <- manning_qdot(swotlist = swotlist, Avar = "Ahat", log = TRUE)
    qbarns[[i]] <- qbarn
    A0hats[[i]] <- A0hat
    
    if (trace) {
      cat(sprintf("iter: \t%s \nQn: \t%s \nA0: \t%s \nndot: \t%s \n\n", 
                  i, 
                  paste(round(qbarn, digits = 3)),
                  paste(round(A0hat), collapse = "\t"), 
                  paste(round(log(ndot), digits = 3), collapse = "\t")))
    } else if (count) {
      cat(ifelse(i %% 10 == 0, i, "."))
    }
  }
  
  out <- list(A0hats = A0hats, Qdots = qdots, ndots = ndots, Qbarn = qbarns)
  out
  
}


estA0_lm_ws35 <- function(ws35mat, dAmat, xvec = NULL, weight = !is.null(xvec),
                          zero = c("same", "first", "minimum", "median"),
                          random_omega = FALSE, intercept = FALSE,
                          symmetric = FALSE) {
  
  zero <- match.arg(zero)
  
  if (!zero == "same") {
    dAmat <- rezero_dA(dAmat, zero = zero)
  } 
  
  lmdf <- estA0_moddf(ws35mat = ws35mat, dAmat = dAmat, xvec = xvec,
                      weight = weight,
                      random_omega = random_omega, symmetric = symmetric)
  
  if (intercept) {
    lmout <- lm(y ~ ., data = lmdf, y = TRUE)
  } else {
    lmout <- lm(y ~ 0 + ., data = lmdf, y = TRUE)
  }
  
  lmout
}

estA0_moddf <- function(ws35mat, dAmat, xvec = NULL, weight = !is.null(xvec),
                        random_omega = FALSE, symmetric = FALSE) {
  
  xmat <- t(ws35mat)
  dAmat <- t(dAmat)
  nx <- ncol(xmat)
  
  if (random_omega) {
    omegas <- lapply(1:100, function(x) rnorm(nx))
    omegas <- lapply(omegas, function(x) x - mean(x))
  } else {
    omegas <- make_omegas(nx, symmetric = symmetric)
    
    if (weight) {
      dists <- abs(map_dbl(omegas, ~(xvec %*% .)))
      omegas <- map2(omegas, dists, ~(.x / .y))
    }
  }
  
  modmats <- lapply(omegas, function(x) xmat %*% diag(x))
  modmat <- Reduce(rbind, modmats)
  moddf <- setNames(as.data.frame(modmat), 
                    paste0("x", 1:ncol(modmat)))
  
  rsps <- lapply(omegas, function(x) (xmat * -dAmat) %*% x) 
  moddf$y <- unlist(rsps)
  
  moddf
}


## Moved to swotr package
# 
# realA0 <- function(swotlist, 
#                    rezero = c("none", "first", "minimum", "median")) {
#   rezero = match.arg(rezero)
#   
#   dA <- swotlist$dA
#   if (rezero != "none") {
#     dA <- rezero_dA(dA, zero = rezero)
#   }
#   
#   out <- (swotlist$A - dA)[, 1]
#   out
# }
