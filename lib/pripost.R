# pripost
# Mark Hagemann
# 2/24/2018
# Scoping out priors vs posteriors vs truth.


#' Caclulate true Manning parameters and derivatives
#' 
#' @param nclist A netcdf-derived list obtained from swotData::nc_list()
#' 
calcManParams_nc <- function(nclist) {
  
  goods <- nclist$River_Info.gdrch
  S <- nclist$Reach_Timeseries.S[goods, ]
  W <- nclist$Reach_Timeseries.W[goods, ]
  A <- nclist$Reach_Timeseries.A[goods, ]
  Q <- nclist$Reach_Timeseries.Q[goods, ]
  
  N <- 1 / Q * W^(-2/3) * A^(5/3) * S^(1/2)
  n <- median(N, na.rm = TRUE)
  
  A0 <- apply(A, 1, min, na.rm = TRUE)
  
  qvec <- apply(Q, 2, median)
  logQbar <- mean(log(Q))
  logQdot <- apply(log(Q), 2, mean) - logQbar
  
  out <- list(n = n, A0 = A0, logQbar = logQbar, logQdot = logQdot, Q = qvec)
  out
}


pripost_n <- function(bampriors, stanfit, true_n = NULL) {
  posts <- rstan::extract(stanfit, pars = c("logn"), permuted = FALSE) %>% 
    reshape2::melt()
  
  nchains <- stanfit@sim$chains
  chain <- 1:nchains
  stopifnot(is.numeric(chain))
  pstats <- posts %>% 
    dplyr::mutate(chains = gsub("^chain:", "", chains)) %>% 
    dplyr::filter(chains %in% chain) %>% 
    dplyr::mutate(value = ifelse(grepl("^log", parameters), exp(value), value),
                  parameters = gsub("^log", "", parameters))
  
  # Make the data frame for plotting. Should include stan samples and values of prior density function
  
  priorlen <- 250
  priormin <- exp(bampriors$lowerbound_logn)
  priormax <- exp(bampriors$upperbound_logn)
  priorx <- seq(priormin, priormax, length.out = priorlen)
  priory <- dlnorm(x = priorx, 
                   meanlog = bampriors$logn_hat, sdlog = bampriors$logn_sd)
  priory <- priory / max(priory) 
  priordf <- data.frame(priorx = priorx, priory = priory)
  
  out <- pstats %>% 
    filter(parameters == "n") %>% 
    ggplot(aes(x = value)) + 
    stat_density(aes(y = ..scaled.., color = "darkblue"), geom = "line") + 
    geom_line(data = priordf, aes(x = priorx, y = priory, color = "magenta")) +
    scale_color_manual(name = "", values = c("magenta" = "magenta", "darkblue" = "darkblue", "black" = "black"), 
                       labels = c("magenta" = "prior", "darkblue" = "posterior", "black" = "truth")) +
    # geom_histogram(aes(y = ..ncount..), position = "identity") +
    xlab("Manning's n") + 
    ylab("Scaled Probability") + theme_bw()
  
  if (!is.null(true_n)) {
    truth <- true_n
    out <- out + geom_vline(aes(xintercept = truth, color = "black"))
  }
    
  out
  
}


pripost_A0 <- function(bampriors, stanfit, true_A0 = NULL) {
  
  posts <- rstan::extract(stanfit, pars = c("A0"), permuted = FALSE) %>% 
    reshape2::melt()
  
  nchains <- stanfit@sim$chains
  chain <- 1:nchains
  stopifnot(is.numeric(chain))
  pstats <- posts %>% 
    dplyr::mutate(chains = gsub("^chain:", "", chains)) %>% 
    dplyr::filter(chains %in% chain) %>% 
    dplyr::mutate(value = ifelse(grepl("^log", parameters), exp(value), value),
                  parameters = gsub("^log", "", parameters)) %>% 
    dplyr::mutate(index = gsub("^.+\\[", "", parameters), 
                  index = gsub("\\]$", "", index), 
                  index = as.numeric(index), 
                  parameters = gsub("\\[.+\\]$", "", parameters)) %>%
    dplyr::arrange(index)
  
  # Make the data frame for plotting. Should include stan samples and values of prior density function
  
  priorlen <- 250
  nx <- length(bampriors$logA0_hat)
  vecmax <- function(x, y) mapply(max, x, y)
  vecmin <- function(x, y) mapply(min, x, y)
  priormin <- vecmax(rep(bampriors$lowerbound_A0, nx),
                     exp(bampriors$logA0_hat - 2.5 * bampriors$logA0_sd))
  priormax <- vecmin(bampriors$upperbound_A0, 
                     exp(bampriors$logA0_hat + 2.5 * bampriors$logA0_sd))
  priorx <- seq(min(priormin), max(priormax), length.out = priorlen)
  priory <- mapply(dlnorm, meanlog = bampriors$logA0_hat, sdlog = bampriors$logA0_sd, 
                   MoreArgs = list(x = priorx))
  maxmat <- matrix(rep(apply(priory, 2, max), nrow(priory)), 
                   nrow = nrow(priory), byrow = TRUE)
  priory <- priory / maxmat
  priordf <- as.data.frame(cbind(priorx, priory)) %>% 
    setNames(c("priorx", 1:ncol(priory))) %>% 
    gather(key = "index", value = "priory", - priorx)
  
  out <- pstats %>% 
    filter(parameters == "A0") %>% 
    ggplot(aes(x = value)) + 
    stat_density(aes(y = ..scaled.., color = "darkblue"), geom = "line") + 
    geom_line(data = priordf, aes(x = priorx, y = priory, color = "magenta")) +
    scale_color_manual(name = "", 
                       values = c("magenta" = "magenta", 
                                  "darkblue" = "darkblue", 
                                  "black" = "black"), 
                       labels = c("magenta" = "prior", 
                                  "darkblue" = "posterior", 
                                  "black" = "truth")) +
    facet_wrap(~index) +     
    xlab("A0 (minimum)") + 
    ylab("Scaled Probability") + theme_bw()
  
  
  
  if (!is.null(true_A0)) {
    truth <- true_A0
    realx <- data.frame(realx = truth, index = 1:length(truth))
    out <- out + geom_vline(data = realx, aes(xintercept = realx, color = "black"))
  }
  
  out
  
}

pripost_qdot <- function(bampriors, stanfit, true_Q = NULL, conf.level = 0.95) {
  
  alpha <- (1 - conf.level) / 2
  q1 <- alpha
  q2 <- 1 - alpha
  
  posts <- rstan::extract(stanfit, pars = c("logQ", "logn"), permuted = FALSE) %>%
    reshape2::melt() %>%
    spread(key = parameters, value = value) %>% 
    gather(key = index, value = logQ, -iterations, -chains, -logn) %>% 
    dplyr::mutate(index = gsub("^.+\\[", "", index), 
           index = gsub("\\]$", "", index), 
           index = as.numeric(index), 
           chains = gsub("^chain:", "", chains)) %>% 
    group_by(chains, iterations) %>% 
    dplyr::mutate(logQbar = mean(logQ),
           logQdot = logQ - logQbar,
           alpha = logn + logQbar)
  
  out <- posts %>% 
    group_by(index) %>% 
    dplyr::summarize(mean = mean(logQdot), 
                     lwr = quantile(logQdot, q1), 
                     upr = quantile(logQdot, q2)) %>%  
    ungroup() %>%
    ggplot(aes(x = index, y = mean)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#dddddd") +
    geom_line(linetype = 2) +
    xlab("Time index") + 
    ylab("Centered log-transformed flow") + theme_bw()
  
  if (!is.null(true_Q)) {
    truth <- log(true_Q) - mean(log(true_Q))
    out <- out + geom_line(aes(y = truth), linetype = 1)
  }
  out
}



# Now logQbar
pripost_qbar <- function(bampriors, bamdata, stanfit, true_Q = NULL) {
  
  posts <- rstan::extract(stanfit, pars = c("logQ", "logn"), permuted = FALSE) %>%
    reshape2::melt() %>%
    spread(key = parameters, value = value) %>% 
    gather(key = index, value = logQ, -iterations, -chains, -logn) %>% 
    dplyr::mutate(index = gsub("^.+\\[", "", index), 
           index = gsub("\\]$", "", index), 
           index = as.numeric(index), 
           chains = gsub("^chain:", "", chains)) %>% 
    group_by(chains, iterations) %>% 
    dplyr::mutate(logQbar = mean(logQ),
           logQdot = logQ - logQbar,
           alpha = logn + logQbar)
  
  pstats <- posts %>% 
    group_by(chains, iterations) %>% 
    dplyr::summarize(logQbar = mean(logQbar), logn = mean(logn), alpha = mean(alpha))
  
  lqhat <- mean(bamdata$logQ_hat, na.rm = TRUE)
  priorlen <- 250
  priormin <- exp(max(bampriors$lowerbound_logQ, 
                  lqhat - 2.5 * max(bampriors$logQ_sd)))
  priormax <- exp(min(bampriors$upperbound_logQ, 
                  lqhat + 2.5 * max(bampriors$logQ_sd)))
  
  
  priorx <- seq(priormin, priormax, length.out = priorlen)
  priory <- dlnorm(priorx, meanlog = lqhat, sdlog = max(bampriors$logQ_sd))
  priory <- priory / max(priory)
  priordf <- data.frame(priorx = priorx, priory = priory)
  
  out <- ggplot() + 
    # geom_density(aes(y = ..scaled.., color = "darkblue")) +
    stat_density(data = pstats, aes(x = exp(logQbar), y = ..scaled.., 
                                    color = "darkblue"), geom = "line") + 
    # geom_histogram(aes(y = ..ncount..), position = "identity") +
    geom_line(data = priordf, aes(x = priorx, y = priory, color = "magenta")) +
    scale_color_manual(name = "", values = c("magenta" = "magenta", "darkblue" = "darkblue", "black" = "black"), 
                         labels = c("magenta" = "prior", "darkblue" = "posterior", "black" = "truth")) +
    xlab("Geometric mean flow") + 
    ylab("Scaled Probability") + theme_bw()
  
  if (!is.null(true_Q)) {
    truth <- exp(mean(log(true_Q)))
    out <- out + geom_vline(aes(xintercept = truth, color = "black"))
  }
  out
  
}


# Now for alpha

pripost_alpha <- function(bampriors, bamdata, stanfit, true_n = NULL, true_Q = NULL) {
  
  posts <- rstan::extract(stanfit, pars = c("logQ", "logn"), permuted = FALSE) %>%
    reshape2::melt() %>%
    spread(key = parameters, value = value) %>% 
    gather(key = index, value = logQ, -iterations, -chains, -logn) %>% 
    dplyr::mutate(index = gsub("^.+\\[", "", index), 
           index = gsub("\\]$", "", index), 
           index = as.numeric(index), 
           chains = gsub("^chain:", "", chains)) %>% 
    group_by(chains, iterations) %>% 
    dplyr::mutate(logQbar = mean(logQ),
           logQdot = logQ - logQbar,
           alpha = logn + logQbar) %>% 
    ungroup()
  
  pstats <- posts %>% 
    group_by(chains, iterations) %>% 
    dplyr::summarize(logQbar = mean(logQbar), logn = mean(logn), alpha = mean(alpha)) %>% 
    ungroup()
  
  lqhat <- mean(bamdata$logQ_hat, na.rm = TRUE)
  priorlen <- 250
  priorsd <- sqrt(bampriors$logn_sd^2 + max(bampriors$logQ_sd)^2)
  priormin <- max(bampriors$lowerbound_logQ + bampriors$lowerbound_logn, 
                  lqhat + bampriors$logn_hat - 2.5 * priorsd)
  priormax <- min(bampriors$upperbound_logQ + bampriors$upperbound_logn, 
                  lqhat + bampriors$logn_hat + 2.5 * priorsd)
  
  priorx <- seq(priormin, priormax, length.out = priorlen)
  priory <- dnorm(priorx, mean = lqhat + bampriors$logn_hat, sd = priorsd)
  priory <- priory / max(priory)
  priordf <- data.frame(priorx = priorx, priory = priory)
  
  out <- ggplot() + 
    stat_density(data = pstats, aes(x = alpha, y = ..scaled.., 
                                    color = "darkblue"), geom = "line") + 
    geom_line(data = priordf, aes(x = priorx, y = priory, color = "magenta")) +
    scale_color_manual(name = "", values = c("magenta" = "magenta", "darkblue" = "darkblue", "black" = "black"), 
                       labels = c("magenta" = "prior", "darkblue" = "posterior", "black" = "truth")) +
    
    # geom_histogram(aes(y = ..ncount..), position = "identity") +
    xlab("mean(log flow) + log(Manning's n)") + 
    ylab("Scaled Probability") + theme_bw()
  
  if (!is.null(true_n) && !is.null(true_Q)) {
    truth <- log(true_n) + mean(log(true_Q))
    out <- out + geom_vline(aes(xintercept = truth, color = "black"))
  }
  
  out
}

pripost_q <- function(bampriors, stanfit, true_Q = NULL, conf.level = 0.95) {
  
  alpha <- (1 - conf.level) / 2
  q1 <- alpha
  q2 <- 1 - alpha
  
  posts <- rstan::extract(stanfit, pars = c("logQ", "logn"), permuted = FALSE) %>%
    reshape2::melt() %>%
    spread(key = parameters, value = value) %>% 
    gather(key = index, value = logQ, -iterations, -chains, -logn) %>% 
    dplyr::mutate(index = gsub("^.+\\[", "", index), 
           index = gsub("\\]$", "", index), 
           index = as.numeric(index), 
           chains = gsub("^chain:", "", chains))
  
  out <- posts %>% 
    dplyr::group_by(index) %>% 
    dplyr::summarize(mean = mean(exp(logQ)),
                     lwr = quantile(exp(logQ), q1),
                     upr = quantile(exp(logQ), q2)) %>% 
    ggplot(aes(x = index, y = mean)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#dddddd") +
    geom_line(linetype = 2) +
    xlab("Time index") + 
    ylab("Discharge (cms)") + theme_bw()
  
  if (!is.null(true_Q)) {
    truth <- true_Q
    out <- out + geom_line(aes(y = truth), linetype = 1)
  }
  out
}


pripost_suite <- function(stanfit, bamdata, swotlist, bampriors = NULL, conf.level = 0.95) { 
  if (is.null(bampriors)) {
    bampriors <- bam_priors(bamdata)
  }
  
  swotlist$dA <- rezero_dA(swotlist$dA, "minimum")
  real_A0 <- realA0(swotlist)
  real_logn <- mean(manning_closure(swotlist, log = TRUE, mc = TRUE), na.rm = TRUE)
  real_n <- exp(real_logn)
  real_Q <- swotlist$Q
  realQvec <- apply(real_Q, 2, geomMean)
  pp_A0 <- pripost_A0(bampriors = bampriors, stanfit = stanfit, true_A0 = real_A0)
  pp_alpha <- pripost_alpha(bampriors = bampriors, bamdata = bamdata, stanfit = stanfit, 
                            true_n = real_n, true_Q = real_Q)
  pp_n <- pripost_n(bampriors = bampriors, stanfit = stanfit, true_n = real_n)
  pp_q <- pripost_q(bampriors = bampriors, stanfit = stanfit, true_Q = realQvec)
  pp_qbar <- pripost_qbar(bampriors = bampriors, bamdata = bamdata, 
                          stanfit = stanfit, true_Q = realQvec)
  pp_qdot <- pripost_qdot(bampriors = bampriors, stanfit = stanfit, true_Q = realQvec)
  
  out <- list(Q = pp_q, Qbar = pp_qbar, Qdot = pp_qdot, 
              logQn = pp_alpha, n = pp_n, A0 = pp_A0)
  out
}