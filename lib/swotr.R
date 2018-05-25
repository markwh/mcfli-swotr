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


#' Subset dimensions of swot-like data.
#' 
#' @param swotlist a list of swot-like matrices
#' @param keeptimes Indices of times (columns) to keep. Default is 0, which 
#'   keeps all times. Negative indices are allowed. 
#' @param keeplocs Indices of locations (rows) to keep. Default is 0, which 
#'   keeps all locations. Negative indices are allowed. 
swot_sset <- function(swotlist, keeptimes = 0L, keeplocs = 0L) {
  
  nr <- nrow(swotlist$W)
  nc <- ncol(swotlist$W)
  
  if (length(keeptimes) == 1 && keeptimes == 0)
    keeptimes <- 1:nc
  if (length(keeplocs == 1) && keeplocs == 0L)
    keeplocs <- 1:nr
  
  out <- lapply(swotlist, `[`, keeplocs, keeptimes)
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

#' Convert a vector to a suitably-dimensioned matrix
#' 
#' @param vec
#' @param pattern a matrix with the desired dimensions
swot_vec2mat <- function(vec, pattern) {
  nr <- nrow(pattern)
  nc <- ncol(pattern)
  
  if (nr == nc) 
    stop("Doesn't work for square matrices")
  
  if (length(vec) == nr) {
    repvec <- rep(vec, nc)
  } else if (length(vec) == nc) {
    repvec <- rep(vec, each = nr)
  } else {
    stop(paste("vec length must be equal to either number",
         "of rows or columns of pattern."))
  }
  
  out <- matrix(repvec, nrow = nr, ncol = nc, byrow = FALSE)
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
           time <= max(time) - max(abs(lag))) %>% 
    select(-lag)
  out <- swot_untidy(swotdf)
  # browser()
  out
}

#' Returns a lagged swotlist, possibly with a new Ahat piece
ccf_lag <- function(swotlist, Ahat = TRUE, verbose = FALSE) {
  
  # initialize dl, the datalist that will be modified
  dl <- swotlist
  # browser()
  for(i in 1:10) {
    # cat(i, "\n")
    ntimes <- ncol(dl$W)
    
    W <- dl$W
    S <- dl$S
    dA <- rezero_dA(dl$dA, "minimum")
    
    # initialize A0
    delta <- 10
    A0 <- estA0(swotlist = dl, zero = "minimum", random_omega = FALSE, ndot = 1)
    if (verbose) cat(A0, "\n")
    
    A0[A0 <= 0] <- min(W[W>0])
    A <- dA + swot_vec2mat(A0, pattern = W)
    
    logmanlist <- (W^(-2/3) * S^(1/2) * A^(5/3)) %>% 
      t() %>% 
      log() %>% 
      as.data.frame()
    
    ccs <- map(2:length(logmanlist), 
               ~ccf(logmanlist$V1, logmanlist[[.]], plot = FALSE, na.action = na.pass))
    bestlags <- c(0L, map_int(ccs, ~as.integer(.$lag)[which.max(.$acf)]))
    if (verbose) cat(bestlags, "\n")
    
    
    if (sum(!bestlags == 0) == 0)
      break
    
    dl <- swot_timelag(dl, bestlags)
  }
  
  out <- dl
  
  if (Ahat)
    out$Ahat <- A
  
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

swot_A <- function(A0vec, dAmat) {
  A0mat <- swot_vec2mat(A0vec, dAmat)
  Amat <- A0mat + dAmat
  Amat
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


# Plot all variables, or a subset thereof, as timeseries

swot_plot <- function(swotlist, vars = "all"){
  if (!(length(vars) == 1 && vars == "all")) {
    swotlist <- swotlist[vars]
  }
  plotdf <- swot_tidy(swotlist) %>% 
    gather(key = variable, value = value, -time, -loc) %>% 
    mutate(loc = as.factor(loc))

  out <- ggplot(plotdf, aes(x = time, y = value, color = loc)) +
    geom_line() +
    geom_point() +
    facet_wrap(~variable, scales = "free_y")
  out
}

plot_DAWG <- function (dawgmat) {
  dawgdf <- as.data.frame(t(dawgmat)) %>% setNames(1:nrow(dawgmat)) %>% 
    dplyr::mutate(time = 1:ncol(dawgmat)) %>% 
    reshape2::melt(id.vars = "time", 
         variable.name = "reach") %>% 
    dplyr::mutate(reach = as.numeric(reach))
  ggplot(dawgdf, aes(x = time, y = value, group = reach)) + 
    geom_line(aes(color = reach, group = reach)) + 
    scale_color_gradient()
}
