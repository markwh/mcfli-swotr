# utils.R
# 3/19/2018
# Mark Hagemann
# Various utility functions

#' @param swotlist A named list of matrices 
swot_tidy <- function(swotlist) {
  nt <- ncol(swotlist[[1]])
  nx <- nrow(swotlist[[1]])
  
  outvecs <- map(swotlist, ~as.vector(.))
  
  times <- rep(1:nt, each = nx)
  locs <- rep(1:nx, nt)
  
  out <- as.data.frame(outvecs)
  out$time <- times
  out$loc <- locs
  out
}

swot_untidy <- function(swotdf) {
  matnames <- setdiff(names(swotdf), c("time", "loc"))
  # browser()
  times <- unique(swotdf$time)
  newtimes <- order(times)
  
  locs <- unique(swotdf$loc)
  newlocs <- order(locs)
  
  nc <- max(newtimes)
  nr <- max(newlocs)
  
  
  swotdf <- arrange(swotdf, time, loc)
  out <- map(swotdf[matnames], 
             ~matrix(., nrow = nr, ncol = nc, byrow = FALSE))
  out
}


swot_timelag <- function(swotlist, lags) {
  # browser()
  mats <- vapply(swotlist, is.matrix, logical(1))
  
  swotdf <- swot_tidy(swotlist[mats])
  
  lagdf <- data.frame(loc = 1:max(swotdf$loc), lag = lags)
  
  swotdf <- left_join(swotdf, lagdf, by = "loc") %>% 
    mutate(time = time + lag) %>% 
    filter(time >= min(time) + max(abs(lag)), 
           time <= max(time) - max(abs(lag)))
  out <- swot_untidy(swotdf)
  # browser()
  out
}

ccf_lag <- function(swotlist, Aest = TRUE) {
  
  # initialize dl, the datalist that will be modified
  dl <- swotlist
  
  for(i in 1:10) {
    # cat(i, "\n")
    ntimes <- ncol(dl$W)
    tomat <- function(x) {
      matrix(rep(x, ntimes), ncol = ntimes)
    }
    
    W <- dl$W
    S <- dl$S
    dA <- rezero_dA(dl$dA, "minimum")
    
    # initialize A0
    delta <- 10
    A0 <- coef(estA0(W, S, dA + delta))
    cat(A0, "\n")
    
    A0[A0 <= 0] <- min(W[W>0])
    A <- dA + tomat(A0)
    
    logmanlist <- (W^(-2/3) * S^(1/2) * A^(5/3)) %>% 
      t() %>% 
      log() %>% 
      as.data.frame()
    
    ccs <- map(2:length(logmanlist), 
               ~ccf(logmanlist$V1, logmanlist[[.]], plot = FALSE, na.action = na.pass))
    bestlags <- c(0L, map_int(ccs, ~as.integer(.$lag)[which.max(.$acf)]))
    cat(bestlags, "\n")
    
    
    if (sum(!bestlags == 0) == 0)
      break
    
    dl <- swot_timelag(dl, bestlags)
  }
  
  out <- dl
  
  if (Aest)
    out$A <- A
  
  out
}



rezero_dA <- function(dAmat, zero = c("first", "minimum", "median")) {
  zero <- match.arg(zero)
  if (zero == "first") {
    shifts <- dAmat[, 1]
  } else if (zero == "minimum") {
    shifts <- apply(dAmat, 1, min, na.rm = TRUE)
  } else if (zero == "median") {
    shifts <- apply(dAmat, 1, median, na.rm = TRUE)
  }
  shiftmat <- matrix(rep(shifts, ncol(dAmat)), ncol = ncol(dAmat))
  out <- dAmat - shiftmat
  out
}

swot_A <- function(dAmat, A0vec) {
  A0mat <- matrix(rep(A0vec, ncol(dAmat)), ncol = ncol(dAmat))
  Amat <- A0mat + dAmat
  Amat
}

manning_qdot <- function(Wmat, Smat, Amat, log = FALSE) {
  outmat <- Amat^(5/3) * Wmat^(-2/3) * Smat^(1/2)
  if (log) 
    outmat <- log(outmat)
  
  outmat
}


# Estimate flow imbalance as percent 

# swot_dqdx <- function(swotlist) {
#   rhsmat <- with(swotlist, manning_qdot(W, S, A, log = FALSE))
#   diffs <- apply(rhsmat, 2, function(x) abs(diff(x)))
#   dxerr <- mean(diffs / rhsmat[-1, ]) / 2 * 100 # times 100 makes it a percent
#   dxerr
# }

# This version uses log.
swot_dqdx <- function(swotlist) {
  rhsmat <- with(swotlist, manning_qdot(W, S, A, log = TRUE))
  sds <- apply(rhsmat, 2, sd, na.rm = TRUE)
  dxerr <- mean(sds) * 100 # times 100 makes it a percent
  dxerr
}

true_dqdx <- function(qmat) {
  sds <- apply(log(qmat), 2, sd, na.rm = TRUE)
  dxerr <- mean(sds) * 100 # times 100 makes it a percent
  dxerr
}

# An R2 function that works as expected for no-intercept models
A0_R2 <- function(mod) {
  true <- mod$model$y
  pred <- predict(mod)
  out <- markstats::R2(true, pred)
  out
}

# Condition number in 3/5 space
condNo <- function(mat) {
  eigvals <- eigen(mat)[["values"]]
  out <- abs(eigvals[1]) / abs(eigvals[length(eigvals)])
  out
}

A0_condno <- function(swotlist) {
  wm35 <- with(swotlist, t(W ^ (-2/3) * S ^ (1/2)) ^ (3/5))
  out <- condNo(cor(wm35))
  out
}
