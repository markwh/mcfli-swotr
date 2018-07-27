# Gibbs sampler functions for McMan

gibbs_inputs <- function(swotlist, Qhat = NULL) {
  # Specify known quantities
  bdat <- swot_bamdata(swotlist, Qhat = Qhat)

  bi <- bamr:::compose_bam_inputs(bdat)
  bi$logQ_hat <- mean(bi$logQ_hat)
  bi$logQ_sd <- mean(bi$logQ_sd)
  bi$dAobs <- rezero_dA(bi$dAobs, zero = "median")
  ws <- -2 / 3 * log(bi$Wobs) + 1 / 2 * log(bi$Sobs)
  A0_min <- -apply(bi$dAobs, 1, min)
  dx <- findif_x(swotlist$x)
  deltax <- swotlist$x - mean(swotlist$x)
  
  out <- c(bi, list(ws = ws, A0_min = A0_min, dx = dx, deltax = deltax))
  out
}

gibbs_init_high <- function(inputs) {
  
  # Initialize parameters
  
  out <- list(
    A0 = 20 * (inputs$A0_min + 2 * max(inputs$Wobs)),
    qn = rep(50, inputs$nt),
    qbar = 100,
    sigsq_epsilon = 1e5,
    sigsq_q = 1e5,
    n = 10,
    dgdx = rep(0, inputs$nt),
    sigsq_dgdx = 1 / max(inputs$dx),
    nubar = rep(0, inputs$nx),
    sigsq_nubar = 1e5
  )
  out$logA <- log(swot_A(out$A0, inputs$dAobs))
  out$gamma <- inputs$dx * swot_vec2mat(out$dgdx, inputs$dx)
  out$nu <- swot_vec2mat(out$nubar, inputs$ws)
  
  out
}

gibbs_init_low <- function(inputs) {
  out <- list(
    A0 =  inputs$A0_min + 0.1,
    qn = rep(-20, inputs$nt),
    qbar = -50,
    sigsq_epsilon = 1e-5,
    sigsq_q = 1e-5,
    n = -15,
    dgdx = rep(0, inputs$nt),
    sigsq_dgdx = 1e-5 / max(inputs$dx),
    nubar = rep(0, inputs$nx),
    sigsq_nubar = 1e-5
  )
  out$logA <- log(swot_A(out$A0, inputs$dAobs))
  out$gamma <- inputs$dx * swot_vec2mat(out$dgdx, inputs$dx)
  out$nu <- swot_vec2mat(out$nubar, inputs$ws)
  
  out
}

gibbs_init_rand <- function(inputs, sd = 2) {
  out <- list(
    A0 =  exp(truncnorm::rtruncnorm(inputs$nx, a = log(inputs$A0_min),
                                b = Inf, mean = inputs$logA0_hat, 
                                sd = sd * inputs$logA0_sd)),
    qn = rnorm(inputs$nt, mean = inputs$logQ_hat + inputs$logn_hat,
               sd = sd * sqrt(inputs$logQ_sd^2 + inputs$logn_sd^2)),
    qbar = rnorm(1, mean = inputs$logQ_hat + inputs$logn_hat,
                 sd = sd * sqrt(inputs$logQ_sd^2 + inputs$logn_sd^2)),
    sigsq_epsilon = abs(rnorm(1, 0, 0.1 * sd)),
    sigsq_q = abs(rnorm(1, 0, sd)),
    n = rnorm(1, inputs$logn_hat, sd * inputs$logn_sd),
    dgdx = rnorm(inputs$nt, 0, 0.1 * sd / max(inputs$dx)),
    sigsq_dgdx = abs(rnorm(1, 0, 0.1 * sd)),
    nubar = rnorm(inputs$nx, 0, 0.1 * sd),
    sigsq_nubar = abs(rnorm(1, 0, 0.1 * sd))
  )
  out$logA <- log(swot_A(out$A0, inputs$dAobs))
  out$gamma <- inputs$dx * swot_vec2mat(out$dgdx, inputs$dx)
  out$nu <- swot_vec2mat(out$nubar, inputs$ws)

  out
}


gibbs_inits <- function(inputs, chains, method = "rand") {
  
  method <- match.arg(method, choices = c("rand", "high", "low"), 
                      several.ok = TRUE)
  
  method <- rep(method, length.out = chains)
  
  print(method)
  
  inits <- purrr::map(method, ~get(paste0("gibbs_init_", .))(inputs))
  inits
}



# Sampling functions ------------------------------------------------------


# Sampler utility functions ------------------------
postmu <- function(primu, prisigsq, likmu, liksigsq, n) {
  mu <- (liksigsq * primu + n * prisigsq * likmu) / (n * prisigsq + liksigsq)
  mu
}

postsigsq <- function(prisigsq, liksigsq, n) {
  sigsq <- (liksigsq * prisigsq) / (n * prisigsq + liksigsq)
  sigsq
}

# Formulas from https://arxiv.org/pdf/1605.01019.pdf
invgam_mom <- function(mean, sd) {
  alpha <- mean^2 / (sd^2) + 2
  beta <- mean * (alpha - 1)
  out <- c(alpha = alpha, beta = beta)
  out
}

# Formulas from http://www.stat.cmu.edu/~brian/463-663/week10/Chapter%2004.pdf
# Also https://www.coursera.org/lecture/mcmc-bayesian-statistics/computing-example-with-normal-likelihood-Ilg9Z
postalpha <- function(prialpha, n) {
  alpha <- prialpha + n / 2
  alpha
}

postbeta <- function(pribeta, y, mu) {
  beta <- pribeta + 1 / 2 * sum((y - mu)^2)
  beta
}

sample_normal <- function(priormu, priorsigsq, likmu, liksigsq, likn) {
  mupost <- postmu(primu = priormu, 
                 prisigsq = priorsigsq,
                 likmu = likmu, 
                 liksigsq = liksigsq, n = likn)
  sigsqpost <- postsigsq(prisigsq = priorsigsq,
                       liksigsq = liksigsq, n = likn)
  out <- rnorm(length(mupost), mupost, sqrt(sigsqpost))
  attr(out, "llik") <- sum(dnorm(out, mean = mupost, sd = sqrt(sigsqpost), log = TRUE))
  out
}

# sample_invgamma <- function(sigmahat, sigmasd, obsvec, obsmean, 
#                             likn = length(obsvec)) {
#   alpha <- postalpha(prialpha = 2.17, n = inputs$nt)
#   beta <- postbeta(pribeta = 0.0206, y = state$dgdx, mu = 0)
#   
#   drawgamma <- rgamma(1, shape = alpha, rate = beta)
#   out <- 1 / drawgamma
#   out
# }

# Sampling functions

sample_qn <- function(inputs, state) {
  priormu <- state$qbar + state$n
  priorsigsq <- 1
  
  likmu <- apply(inputs$ws + 5 / 3 * state$logA - state$gamma - state$nu, 
                 2, mean)
  liksigsq <- state$sigsq_epsilon
  likn <- inputs$nx
  
  out <- sample_normal(priormu, priorsigsq, likmu, liksigsq, likn)
  out
}

sample_qbar <- function(inputs, state) {
  priormu <- inputs$logQ_hat
  priorsigsq <- (inputs$logQ_sd)^2
  
  likmu <- mean(state$qn) - state$n
  liksigsq <- var(state$qn)
  # liksigsq <- state$sigsq_q
  likn <- inputs$nt
  
  out <- sample_normal(priormu, priorsigsq, likmu, liksigsq, likn)
  out
}

sample_n <- function(inputs, state) {
  priormu <- inputs$logn_hat
  priorsigsq <- inputs$logn_sd
  
  likmu <- mean(state$qn - state$qbar)
  liksigsq <- var(state$qn)
  # liksigsq <- state$sigsq_q
  likn <- inputs$nt
  
  out <- sample_normal(priormu, priorsigsq, likmu, liksigsq, likn)
  out
}

sample_A0 <- function(inputs, state) {
  
  # Sample the matrix of a = logA
  minmat <- 5 / 3 * log(rezero_dA(inputs$dAobs, zero = "minimum") + 0.1)
  minvec <- as.vector(minmat)
  
  qnmat <- swot_vec2mat(state$qn, inputs$ws)
  muvec <- as.vector(qnmat - inputs$ws + state$gamma + state$nu)
  
  a53vec <- truncnorm::rtruncnorm(inputs$nx * inputs$nt, 
                                  a = minvec, b = Inf,
                                  mean = muvec, 
                                  sd = sqrt(state$sigsq_epsilon))
  amat <- matrix(a53vec * 3 / 5, nrow = inputs$nx)
  A0est <- exp(amat) - inputs$dAobs
  
  priormu <- inputs$logA0_hat
  priorsigsq <- (inputs$logA0_sd)^2
  
  likmu <- apply(log(A0est), 1, mean)
  
  liksigsq <- (3 / 5)^2 * state$sigsq_epsilon
  likn <- inputs$nt
  
  mupost <- postmu(primu = priormu,
                   prisigsq = priorsigsq,
                   likmu = likmu,
                   liksigsq = liksigsq, n = likn)
  sigsqpost <- postsigsq(prisigsq = priorsigsq,
                         liksigsq = liksigsq, n = likn)
  
  out <- exp(rnorm(inputs$nx, mean = mupost, sd = sqrt(sigsqpost)))
  attr(out, "llik") <- sum(dnorm(log(out), mupost, sqrt(sigsqpost), 
                                 log = TRUE))
  out
}

sample_sigsq_epsilon <- function(inputs, state) {
  qnmat <- swot_vec2mat(state$qn, inputs$ws)
  errmat <- qnmat - 5 / 3 * state$logA - inputs$ws + state$gamma + state$nu
  errvec <- as.vector(errmat)
  
  alpha <- postalpha(prialpha = 2, n = inputs$nx * inputs$nt)
  beta <- postbeta(pribeta = 0.01, y = errvec, mu = 0)
  
  drawgamma <- rgamma(1, shape = alpha, rate = beta)
  out <- 1 / drawgamma
  
  attr(out, "llik") <- sum(dgamma(drawgamma, shape = alpha, rate = beta, 
                                  log = TRUE))
  out
}

sample_sigsq_q <- function(inputs, state) {

  alpha <- postalpha(prialpha = 7, n = inputs$nt) # priors based on mu=0.7, sd=0.3
  beta <- postbeta(pribeta = 4.5, y = state$qn, mu = (state$qbar + state$n))
  
  drawgamma <- rgamma(1, shape = alpha, rate = beta)
  out <- 1 / drawgamma
  attr(out, "llik") <- sum(dgamma(drawgamma, shape = alpha, rate = beta, 
                                  log = TRUE))
  out
}

# Extra parameters to explain more of the error
sample_dgdx <- function(inputs, state) {
  gamma <- inputs$ws + 5 / 3 * state$logA -
    state$nu - swot_vec2mat(state$qn, inputs$dx)
  dgmat <- findif_x(gamma)
  
  priormu <- 0
  priorsigsq <- state$sigsq_dgdx
  
  likmu <- apply(dgmat / inputs$dx, 2, mean)
  liksigsq <- state$sigsq_epsilon / (mean(inputs$dx[, 1])^2 * 2)
  likn <- inputs$nx
  
  out <- sample_normal(priormu, priorsigsq, likmu, liksigsq, likn)
  out
}

sample_sigsq_dgdx <- function(inputs, state) {
  # priors based on mu=0.0176, sd=0.0425
  alpha <- postalpha(prialpha = 2.17, n = inputs$nt)
  beta <- postbeta(pribeta = 0.0206, y = state$dgdx, mu = 0)
  
  drawgamma <- rgamma(1, shape = alpha, rate = beta)
  out <- 1 / drawgamma
  attr(out, "llik") <- sum(dgamma(drawgamma, shape = alpha, rate = beta, 
                                  log = TRUE))
  out
}

sample_nubar <- function(inputs, state) {
  nu <- inputs$ws + 5 / 3 * state$logA -
    state$gamma - swot_vec2mat(state$qn, inputs$dx)
  
  priormu <- 0
  priorsigsq <- state$sigsq_nubar
  
  likmu <- apply(nu, 1, mean)
  liksigsq <- state$sigsq_epsilon
  likn <- inputs$nt
  
  out <- sample_normal(priormu, priorsigsq, likmu, liksigsq, likn)
  out
}

sample_sigsq_nubar <- function(inputs, state) {
  # # priors based on mu=0.017, sd=0.021
  # alpha <- postalpha(prialpha = 2.65, n = inputs$nx)
  # beta <- postbeta(pribeta = 0.028, y = state$nubar, mu = 0)
  
  # These priors artificially tightened (mu=0.002, sd=0.001)
  alpha <- postalpha(prialpha = 6.00, n = inputs$nx)
  beta <- postbeta(pribeta = 0.01, y = state$nubar, mu = 0)

  
  drawgamma <- rgamma(1, shape = alpha, rate = beta)
  out <- 1 / drawgamma
  attr(out, "llik") <- sum(dgamma(drawgamma, shape = alpha, rate = beta, 
                                  log = TRUE))
  out
}


# Sample a single chain
gibbs_sample_chain <- function(inputs, inits, iter, thin = 1,
                               gamma = FALSE, nu = FALSE,
                               verbose = TRUE, printint = 100, 
                               chain = 1, chains = 1) {
  if (verbose) {
    cat(sprintf("Sampling chain %s of %s\n", chain, chains))
  }
  
  # Allocate chain
  nx <- inputs$nx
  nt <- inputs$nt
  A0chain <- as.data.frame(map(1:nx, ~numeric(iter)), 
                           col.names = sprintf("A0[%s]", 1:nx),
                           optional = TRUE) %>% 
    as.matrix()
  qchain <- as.data.frame(map(1:nt, ~numeric(iter)),
                          col.names = sprintf("logQ[%s]", 1:nt),
                          optional = TRUE) %>% 
    as.matrix()
  qnchain <- as.data.frame(map(1:nt, ~numeric(iter)),
                           col.names = sprintf("qn[%s]", 1:nt),
                           optional = TRUE) %>% 
    as.matrix()
  qbarchain <- numeric(iter)
  ssqep_chain <- numeric(iter)
  ssqq_chain <- numeric(iter)
  dgdxchain <- as.data.frame(map(1:nt, ~numeric(iter)),
                             col.names = sprintf("dg[%s]", 1:nt),
                             optional = TRUE) %>% 
    as.matrix()
  ssqdgdx_chain <- numeric(iter)
  nubarchain <- as.data.frame(map(1:nx, ~numeric(iter)),
                              col.names = sprintf("nu[%s]", 1:nx),
                              optional = TRUE) %>% 
    as.matrix()
  ssqnubar_chain <- numeric(iter)
  nchain <- numeric(iter)
  lpchain <- numeric(iter)
  
  state <- inits
  
  # Sampling
  printint <- printint * thin
  for (i in 1:(iter * thin)) {

    draw <- (i %% thin == 0)
    if (verbose && (i %% printint == 0)) {
      # browser()
      cat(i / thin, "\n")
    }
    
    state$qbar <- sample_qbar(inputs, state)
    state$qn <- sample_qn(inputs, state)
    state$n <- sample_n(inputs, state)
    state$A0 <- sample_A0(inputs, state)
    state$logA <- log(swot_A(A0vec = state$A0, dAmat = inputs$dAobs))
    state$sigsq_epsilon <- sample_sigsq_epsilon(inputs, state)
    state$sigsq_q <- sample_sigsq_q(inputs, state)

    # Additional parameters (optional)
    if (gamma) {
      state$dgdx <- sample_dgdx(inputs, state)
      state$gamma <- state$dgdx * inputs$deltax
      state$sigsq_dgdx <- sample_sigsq_dgdx(inputs, state)
    } else {
      state$dgdx <- rep(0, inputs$nt)
      state$gamma <- matrix(0, nrow = inputs$nx, ncol = inputs$nt)
      state$sigsq_dgdx <- NA
    }
    
    if (nu) {
      state$nubar <- sample_nubar(inputs, state)

      state$sigsq_nubar <- sample_sigsq_nubar(inputs, state)
    } else {
      state$nubar <- rep(0, inputs$nx)
      attr(state$nubar, "llik") <- 0
      state$nu <- swot_vec2mat(state$nubar, inputs$ws)
      state$sigsq_nubar <- NA
    }
    
    if (draw) {
      ci <- i / thin
      
      qbarchain[ci] <- state$qbar
      qnchain[ci, ] <- state$qn
      nchain[ci] <- state$n
      qchain[ci, ] <- state$qn - state$n
      A0chain[ci, ] <- state$A0
      ssqep_chain[ci] <- state$sigsq_epsilon
      ssqq_chain[ci] <- state$sigsq_q
      
      sampparams <- c("A0", "qn", "qbar", "sigsq_epsilon", "sigsq_q", "n")

      if (gamma) {
        dgdxchain[ci, ] <- state$dgdx
        ssqdgdx_chain[ci] <- state$sigsq_dgdx
        sampparams <- c(sampparams, "dgdx", "sigsq_dgdx")
      }
      if (nu) {
        nubarchain[ci, ] <- state$nubar
        ssqnubar_chain[ci] <- state$sigsq_nubar
        sampparams <- c(sampparams, "nubar", "sigsq_nubar")
      }
      # Calculate log posterior
      totlp <- sum(map_dbl(state[sampparams], ~attr(., "llik")))
      lpchain[ci] <- totlp
    }
  }
  
  out <- list(
    A0 = as.data.frame(A0chain),
    logn = nchain,
    logQ = as.data.frame(qchain),
    logQbar = qbarchain,
    sigsq_epsilon = ssqep_chain,
    sigsq_q = ssqq_chain,
    lp__ = lpchain # log posterior
  )
  
  if (gamma) {
    out <- c(out, list(
      dgdx = as.data.frame(dgdxchain),
      sigsq_dgdx = ssqdgdx_chain   
    ))
  }
  if (nu) {
    out <- c(out, list(
      nubar = as.data.frame(nubarchain),
      sigsq_nubar = ssqnubar_chain
    ))
  }
  out
}


# Sample several chains in parallel
gibbs_sample <- function(inputs, chains, init_method = "rand", 
                         iter, thin = 1, 
                         gamma = FALSE, nu = FALSE, 
                         serial = FALSE, 
                         cores = parallel::detectCores(),
                         printint = 100, verbose = TRUE) {

  inits <- gibbs_inits(inputs = inputs, chains = chains, method = init_method)

  sampfun <- function(i) {
    out <- gibbs_sample_chain(inputs = inputs, inits = inits[[i]],
                              iter = iter, thin = thin,
                              gamma = gamma, nu = nu,
                              printint = printint, 
                              verbose = verbose, chain = i,
                              chains = length(inits))
    out
  }
  
  if (serial) {
    samples <- lapply(1:length(inits), sampfun)
  } else {
    clust <- parallel::makeCluster(min(cores, length(inits)), outfile = "")
    on.exit(parallel::stopCluster(clust))
    exportfuns <- list("gibbs_sample_chain", "postmu", "postsigsq", "postalpha",
                       "postbeta", "sample_normal", "sample_qn", "sample_n",
                       "sample_qbar", "sample_A0", "sample_dgdx", 
                       "sample_nubar", "sample_sigsq_nubar",
                       "sample_sigsq_epsilon", "%>%", "findif_x",
                       "sample_sigsq_q", "sample_sigsq_dgdx", "map", "map_dbl")
    parallel::clusterEvalQ(clust, library("swotr"))
    
    parallel::clusterExport(clust, exportfuns)
    
    samples <- parallel::parLapply(cl = clust, X = 1:length(inits), fun = sampfun)
  }
  attr(samples, "thin") <- thin

  out <- samples
  out
}


gibbs_as_dfs <- function(samps) {
  numers <- sapply(samps[[1]], is.numeric)
  convfun <- function(x, name) {
    if (is.numeric(x)) {
      return(setNames(data.frame(x), name))
    } else {
      return(x)
    }
  }
  dflist <- map(samps, ~map2(., names(.), ~convfun(.x, .y))) %>% 
    map(~Reduce(cbind, .))
  dflist
}

#' Create a stan-like object from gibbs sampler output
#' 
#' This allows gibbs sampler output to be analyzed using the 
#' utilities in swotr and bayesplot packages. 
#' 
#' @param samps output from gibbs sampler
#' @param warmup How many samples to use for warmup? Default is half the total. 
gibbs_stanfit <- function(samps, warmup = NULL) {
  
  thin <- attr(samps, "thin")
  
  sampdfs <- gibbs_as_dfs(samps)
  iter <- nrow(sampdfs[[1]])
  
  if (is.null(warmup)) {
    warmup = ceiling(iter / 2)
  }
  meanpars <- lapply(sampdfs, function(x) apply(x[-1:-warmup, ], 2, mean))
  # meanpardf <- Reduce(cbind, meanpars0)
  # meanpars <- apply(meanpardf, 1, mean)
  
  dimfun <- function(x) {
    if (is.data.frame(x)) {
      return(length(x))
    } else {
      return(numeric(0))
    }
  }
  chains <- length(samps)
  
  sparams <- list(accept_stat__ = runif(iter),
                        stepsize__ = rep(1, iter),
                        treedepth__ = rep(1, iter),
                        n_leapfrog__ = rep(1, iter),
                        divergent__ = rep(0, iter),
                        energy__ = rep(1, iter))
  
  samplists <- map2(sampdfs, meanpars, 
                    ~structure(as.list(.x), mean_pars = .y,
                               sampler_params = sparams))
  
  sim <- list(
    samples = samplists, # May need to add attributes here
    iter = iter,
    thin = 1L,
    warmup = warmup,
    chains = chains,
    n_save = rep(iter, chains), # TODO
    warmup2 = rep(warmup, chains), # TODO
    permutation = lapply(1:chains, function(x) sample.int(iter - warmup)), # TODO
    pars_oi = names(samps[[1]]), # TODO
    dims_oi = lapply(samps[[1]], dimfun), # TODO
    fnames_oi = names(sampdfs[[1]]), # TODO
    n_flatnames = NULL # TODO
  )

  stanargs <- lapply(1:chains, function(x) {
    list(chain_id = x, iter = iter, thin = 1L, seed = NA,
         warmup = warmup, init = NA, algorithm = "gibbs", method = "sampling")
  })
  
  out <- new("stanfit",
             model_name = NA_character_,
             model_pars = names(samps[[1]]),
             par_dims = lapply(samps[[1]], dimfun),
             mode = 0L,
             sim = sim,
             inits = list(NA), # TODO
             stan_args = stanargs, 
             stanmodel = structure(NA, class = "stanmodel"),
             date = date(),
             .MISC = new.env())
  out
}
