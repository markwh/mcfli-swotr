# McMan closure term functions. 


# Calculates closure term from Manning equation. 
# Naively this could be considered Manning's n, but it varies in time and space
# Also includes model error and (by default) flow imbalance error. 
manning_closure <- function(swotlist, log = FALSE, mc = TRUE, mcfun = mean) {
  
  W <- swotlist$W
  S <- swotlist$S
  A <- swotlist$A
  Q <- swotlist$Q
  
  if (mc) {
    Qvec <- apply(Q, 2, mcfun)
    Q <- swot_vec2mat(Qvec, Q)
  }
  
  out <- A ^ (5/3) * W ^ (-2/3) * S ^ (1/2) / Q
  
  if (log) 
    out <- log(out)
  
  out
}


# Produces a matrix that comports with McMan in linear A space. 
# Obsmat should be a WS35 matrix, and by default will be calculated as ws35 for 
# the supplied swotlist. The output will be coerced to have the same row-by-row 
# mean as obsmat. 

manning_linA_closed <- function(swotlist, obsmat = NULL, adjust = TRUE) {
  if (is.null(obsmat)) {
    obsmat <- manning_ws35(swotlist = swotlist)
  }
  clos <- manning_closure(swotlist)
  xmat <- obsmat * clos ^ (-3/5)
  
  rmeans_W <- apply(obsmat, 2, mean)
  rmeans_X <- apply(xmat, 2, mean)
  
  if (adjust) {
    out <- xmat * swot_vec2mat(rmeans_W / rmeans_X, xmat)
  } else {
    out <- xmat
  }
  
  out
}

# Produces the zero-mean error matrix in linear A space. 
manning_linA_closure <- function(swotlist, obsmat = NULL) {
  if (is.null(obsmat)) {
    obsmat <- manning_ws35(swotlist = swotlist)
  }
  
  closedmat <- manning_linA_closed(swotlist = swotlist, obsmat = obsmat)
  out <- closedmat - obsmat
  out
}



# Beta functions below--need testing and documenting ----------------------


# from 4/11 notebook
# Update 4/17-4/18: now log-transformed by default, name changed to reflect future packaging

#' gamma matrix from McMan error decomposition
manning_gamma <- function(swotlist, log = TRUE) {
  Qmat <- swotlist$Q
  gmeans <- apply(Qmat, 2, geomMean)
  
  if (sum(apply(Qmat, 2, function(x) length(unique(x)) == 1)) == ncol(Qmat))
    warning("swotlist doesn't seem to have space-varying Q\n")
  
  out <- Qmat / swot_vec2mat(gmeans, Qmat)
  
  if (log) {
    out <- log(out)
  }
  
  out
}

#' nu matrix from McMan error decomposition
manning_nu <- function(swotlist, log = TRUE) {
  nmat <- manning_closure(swotlist, log = TRUE, mc = FALSE)
  out <- nmat - mean(nmat)
  
  if (!log) {
    out <- exp(out)
  }
  
  out
}

# from 4/17 notebook
#' Returns 2 matrices, gammahat and gammaerr, and one vector, 
#' dgdx, from McMan gamma decomposition.
#' 
decomp_gamma <- function (gammamat, xmat) {
  
  dgdx <- findif_x(gammamat) %>% 
    `/`(findif_x(xmat)) %>% 
    apply(2, mean) %>% 
    swot_vec2mat(pattern = xmat)
  
  deltax <- apply(xmat, 2, function(x) x - mean(x))
  gammahat <- deltax * dgdx
  
  gammaerr <- gammamat - gammahat
  
  out <- list(gammahat = gammahat, gammaerr = gammaerr, dgdx = dgdx[1, ])
  out
}

#' Returns 2 matrices, nuhat and nuerr, and one vector, nubar, from McMan nu decomposition
decomp_nu <- function(numat) {
  nubars <- apply(numat, 1, mean)
  nubarmat <- swot_vec2mat(nubars, numat)
  residmat <- numat - nubarmat
  
  out <- list(nuhat = nubarmat, nuerr = residmat, nubar = nubars)
  out
}


decomp_gamma_var <- function (gammamat, xmat) {
  dg1 <- decomp_gamma(gammamat, xmat)
  # out <- c(
  #   gamma = var(as.vector(gammamat)), 
  #   gammahat = var(as.vector(dg1$gammahat)),
  #   resid = var(as.vector(dg1$gammaerr))
  # )
  
  out <- var(as.vector(dg1$gammaerr)) / var(as.vector(gammamat))
  out
}



# Characterize closure term, as parameterized in 4/19 notebook
characterize_closure <- function(swotlist) {
  gma <- manning_gamma(swotlist)
  gd <- decomp_gamma(gma, swotlist$x)
  
  sig_dgdx <- sd(gd$dgdx)
  
  nu <- manning_nu(swotlist)
  nd <- decomp_nu(nu)
  
  sig_alpha <- sd(nd$nuhat[, 1])
  sig_err <- sd(gd$gammaerr + nd$nuerr)
  dx <- sd(swotlist$x[, 1])
  
  out <- data.frame(dgdx = sig_dgdx, alpha = sig_alpha, err = sig_err, 
                    dx = dx, dQ_pct = sig_dgdx * dx * 100)
  out
}


closure_lm <- function(swotlist) {
  clos <- manning_closure(swotlist, log = TRUE, mc = TRUE)
  modmat <- swot_tidy(list(clos = clos, x = swotlist$x)) %>% 
    group_by(time) %>% 
    mutate(xbar = mean(x), xdev = x - xbar) %>% 
    ungroup() %>% 
    mutate(time = as.factor(time), loc = as.factor(loc))
  contrasts(modmat$loc) <- contr.sum
  
  mod <- lm(clos ~ time : xdev + loc, data = modmat)
  mod
} 

closure_coefs <- function(swotlist) {
  mod <- closure_lm(swotlist)
  coefs <- coef(mod)
  nubar0 <- coefs[2:(nrow(swotlist$W))]
  logn <- coefs[1]
  nubar <- c(nubar0, -sum(nubar0))
  dgdx <- coefs[-1:-nrow(swotlist$W)]
  
  out <- list(logn = logn, nubar = nubar, dgdx = dgdx)
  out
}

closure_anova <- function(swotlist) {
  mod <- closure_lm(swotlist)
  out <- broom::tidy(mod)
  out
}
