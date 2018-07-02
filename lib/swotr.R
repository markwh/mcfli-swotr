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
  
  attr(out, "QWBM") <- attr(swotlist, "QWBM")
  
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
  
  attr(out, "QWBM") <- attr(swotlist, "QWBM")
  
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
  
  attr(out, "QWBM") <- attr(swotdf, "QWBM")
  
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
  
  attr(out, "QWBM") <- attr(swotlist, "QWBM")
  
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
  
  attr(out, "QWBM") <- attr(swotlist, "QWBM")
  
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

# Purge all times containing NA's
swot_purge_nas <- function(swotlist, purge = c("times", "locs")) {
  purge = match.arg(purge)
  getnainds <- function(mat) {
    which(is.na(mat), arr.ind = TRUE)
  }
  inddf <- purrr::map(swotlist, getnainds) %>% 
    map(as.data.frame) %>% 
    dplyr::bind_rows() %>% 
    unique()
  
  if (nrow(inddf) == 0L) {
    return(swotlist)
  }
  
  if (purge == "times") {
    out <- swot_sset(swotlist, keeptimes = -unique(inddf[[2]]))
  } else {
    out <- swot_sset(swotlist, keeplocs = -unique(inddf[[1]]))
  }
  
  attr(out, "QWBM") <- attr(swotlist, "QWBM")
  
  out
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

#' Read a netcdf file to a list.
#'
#' @param file string providing the location of a netcdf file.
#' @export
nc_list <- function(file) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("The ncdf4 package is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  nc <- ncdf4::nc_open(file)
  on.exit(ncdf4::nc_close(nc))
  
  vars <- names(nc$var)
  
  out <- setNames(lapply(vars, ncdf4::ncvar_get, nc = nc), make.names(vars))
  out
}

nc_reach <- function (file, good_only = FALSE) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("The ncdf4 package is needed for this function to work. Please install it.", 
         call. = FALSE)
  }
  nclist <- nc_list(file)
  t <- as.Date(nclist$Reach_Timeseries.t - 1, origin = "0000-01-01")
  good_reaches <- nclist$River_Info.gdrch
  W <- nclist$Reach_Timeseries.W
  H <- nclist$Reach_Timeseries.H
  S <- nclist$Reach_Timeseries.S
  dA <- calcdA_mat(w = W, h = H)
  A <- nclist$Reach_Timeseries.A
  QWBM <- nclist$River_Info.QWBM[1]
  inbounds <- 1:nrow(W)
  if (good_only) 
    inbounds <- good_reaches
  Q <- nclist$Reach_Timeseries.Q
  # Qobs <- apply(Q, 2, median)
  
  rbnd <- as.vector(nclist$River_Info.rch_bnd)
  nreach <- length(rbnd) - 1
  upbnd = rbnd[1:nreach][inbounds]
  dnbnd = rbnd[2:(nreach + 1)][inbounds]
  x <- (upbnd + dnbnd) / 2
  
  ptrn <- W[inbounds, ]
  
  out <- list(W = W[inbounds, ], S = S[inbounds, ], dA = dA[inbounds, ], 
              H = H[inbounds, ], 
              t = swot_vec2mat(t, ptrn), x = swot_vec2mat(x, ptrn),
              reachid = swot_vec2mat(inbounds, ptrn), Q = Q[inbounds, ], 
              A = A[inbounds, ])
  attr(out, "QWBM") <- QWBM
  out
}


#' Create bamdata object from a swotlist. 
#' 
#' @param swotlist a list of SWOT observables
#' @param Qhat Prior guess of mean discharge. If NULL, will attempt to get 
#'   from \code{attr(swotlist, "QWBM")}
#' @export
swot_bamdata <- function(swotlist, Qhat = NULL, ...) {
  if (!requireNamespace("bamr", quietly = TRUE)) {
    stop("The bamr package is needed for this function to work. Please install it.", 
         call. = FALSE)
  } 
  
  if (is.null(Qhat)) {
    Qhat <- attr(swotlist, "QWBM")
    if (is.null(Qhat)) {
      stop("QWBM must be supplied if no QWBM attribute is present in swotlist.\n")
    }
  }
  
  bd <- bamr::bam_data(w = swotlist$W, s = swotlist$S, dA = swotlist$dA, Qhat = Qhat, ...)
  
}



