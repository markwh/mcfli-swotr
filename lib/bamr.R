# bamr.R
# Mark Hagemann
# 04/03/2018
# Functions to put into the bamr package, when I'm done here


#' convert a list of bamr inputs to a list suitable for pared model input. 
pare_baminps <- function(baminplist) {
  out <- with(baminplist, list(
    ns = nx,
    nt = nt,
    x = 1/2 * log(Sobs) - 2/3 * log(Wobs),
    dA = rezero_dA(dAobs, "minimum"),
    dA_shift = dA_shift,
    sigma_err = median(sigma_man),
    mu_hat = logQ_hat + logn_hat,
    mu_sd = sqrt(logQ_sd^2 + logn_sd^2),
    logA0_hat = logA0_hat,
    logA0_sd = logA0_sd
  ))
  out
}

#' Add closure characterization parameters to input list
#' 
#' @param inplist List of Stan inputs
#' @param swotlist a swotlist of matrices
#' @param tenkm Represent distances in units of 10 km? This assumes swotlist uses 1-m units. 
#' @param method Passed to characterize_closure
#' 
add_closure_char <- function(inplist, swotlist, tenkm = TRUE, 
                             method = c("decomp", "anova")) {
  method <- match.arg(method)
  if (tenkm) {
    swotlist$x <- swotlist$x / 10000
  }
  
  char <- characterize_closure(swotlist = swotlist, method = method)
  
  out <- inplist
  out$dist_km <- swotlist$x[, 1]
  out$sigma_dgdx <- char$dgdx
  out$sigma_nubar <- char$nuhat
  out$sigma_err <- char$err
  
  out
}

#' Augment a swotcase with posterior mean from BAM
#' 
#' @param swotlist a list of SWOT-like observations
#' @param stanfit a stanfit object, as produced by \code{bam_estimate}
#' 
meancase <- function(swotlist, stanfit) {
  out <- swotlist
  A0mat <- get_posterior_mean(stanfit, pars = "A0")
  A0 <- A0mat[, ncol(A0mat)]
  out$A <- swot_A(A0, out$dA)
  
  qmat <- get_posterior_mean(stanfit, pars = "logQ")
  qvec <- bam_qpred(stanfit)$mean
  
  out$Q <- swot_vec2mat(qvec, out$A)
  
  out
}