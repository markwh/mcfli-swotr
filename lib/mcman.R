# mcman.R
# Mark Hagemann
# 3/27/2018
# Migrated from other files, eventually to be packaged up.

geomMean <- function(x, na.rm = FALSE) {
  exp(mean(log(x), na.rm = na.rm))
}

manning_ws35 <- function(swotlist, ndot = 1) {
  
  wmat <- swotlist$W^(-2/3)
  smat <- swotlist$S^(1/2)
  
  if (length(ndot) > 1 && !is.matrix(ndot)) {
    ndot <- swot_vec2mat(ndot, wmat)
  }
  
  out <- (wmat * smat / ndot) ^ (3/5)
  out
}



# Condition number in 3/5 space
condno <- function(mat, center = TRUE, scale = FALSE) {
  scmat <- scale(mat, center = center, scale = scale)
  sqmat <- t(scmat) %*% scmat
  eigvals <- eigen(sqmat)[["values"]]
  out <- sqrt(eigvals[1]) / sqrt(eigvals[length(eigvals)])
  out
}

ws35_condno <- function(swotlist, center = TRUE, scale = FALSE) {
  ws35 <- t(manning_ws35(swotlist))
  out <- condno(ws35, center = center, scale = scale)
  out
}

A0_condno <- function(mod, center = FALSE, scale = FALSE) {
  modmat <- mod$model[, -1]
  out <- condno(modmat, center = center, scale = scale)
  out
}




# Manning eval functions --------------------------------------------------


# An R2 function that works as expected for no-intercept models
A0_R2 <- function(mod) {
  true <- mod$model[[1]]
  pred <- predict(mod)
  out <- 1 - sum((true - pred)^2)/sum((true - mean(true))^2)
  out
}


#' Manning equation fit
#' 
#' @param swotlist A DAWG-formatted list of swot variables
#' @param mc Impose steady-state mass conservation?

val_linA_df <- function(swotlist, mc = TRUE) {
  W <- swotlist$W
  S <- swotlist$S
  A <- swotlist$A
  Q <- swotlist$Q
  
  if (mc) {
    Qvec <- apply(Q, 2, median, na.rm = TRUE)
    Q <- swot_vec2mat(Qvec, Q)
  }
  
  rhs <- as.data.frame(t((W^(-2/3) * S^(1/2))^(3/5))) %>% 
    setNames(1:ncol(.)) %>% 
    mutate(var = "WS35", time = 1:nrow(.))
  
  lhs <- as.data.frame(t(Q^(3/5) / A)) %>% 
    setNames(1:ncol(.)) %>% 
    mutate(var = "Q35.A", time = 1:nrow(.))
  
  out <- rbind(rhs, lhs) %>% 
    gather(key = loc, value = value, -var, -time) %>% 
    mutate(loc = as.numeric(loc)) %>% 
    spread(key = var, value = value)
  out
}

val_linA_plot <- function(swotlist, mc = TRUE) {
  plotdf <- val_linA_df(swotlist = swotlist, mc = mc) %>% 
    mutate(loc = as.factor(loc))
  
  out <- ggplot(plotdf, aes(x = WS35, y = Q35.A)) +
    geom_point(aes(color = loc)) +
    xlab(expression(paste("(W"^"-2/3", "S"^"1/2", ")"^"3/5"))) + 
    ylab(expression(paste("Q"^"3/5", "A"^"-1")))
  out
}

val_linA_lm <- function(swotlist, mc = TRUE) {
  
  df1 <- val_linA_df(swotlist = swotlist, mc = mc)
  dfs <- split(df1, f = df1$loc)
  lms <- map(dfs, ~lm(Q35.A ~ 0 + WS35, data = ., y = TRUE))
  lms
}

val_linA_R2 <- function(swotlist, mc = TRUE) {
  
  lms <- val_linA_lm(swotlist, mc = mc)
  
  out <- map_dbl(lms, A0_R2)
  out
}

val_log_lm <- function(swotlist, mc = TRUE) {
  
  W <- swotlist$W
  S <- swotlist$S
  A <- swotlist$A
  Q <- swotlist$Q
  
  if (mc) {
    Qvec <- apply(Q, 2, median, na.rm = TRUE)
    Q <- matrix(rep(Qvec, nrow(Q)), nrow = nrow(Q), byrow = TRUE)
  }
  
  rhs <- as.data.frame(t(manning_qdot(swotlist, Avar = "A", log = TRUE)))
  lhs <- as.data.frame(t(log(Q)))
  
  dfs <- map2(rhs, lhs, function(x, y) data.frame(x = x, y = y))
  lms <- map(dfs, ~lm(y ~ x, data = ., y = TRUE))
}

val_log_R2 <- function(swotlist, mc = TRUE) {
  
  lms <- val_log_lm(swotlist, mc = mc)
  
  out <- map_dbl(lms, A0_R2)
  out
}

val_plot <- function(swotlist, mc = TRUE, plot = TRUE, log = TRUE) {
  lms <- val_log_lm(swotlist, mc = mc)
  names(lms) <- 1:length(lms)
  valdfs <- map(lms, ~data.frame(true = .$y, pred = predict(.)))
  valdf <- bind_rows(valdfs, .id = "loc")
  valdf$loc <- as.factor(as.integer(valdf$loc))
  
  if (!log) {
    valdf$true <- exp(valdf$true)
    valdf$pred <- exp(valdf$pred)
  }
  
  if (!plot) 
    return(valdf)
  
  out <- ggplot(valdf, aes(x = pred, y = true)) +
    geom_point(aes(color = loc)) +
    geom_abline(aes(slope = 1, intercept = 0))
  out
}

#' Generate termplots for Manning variables
#'
#' Modified from a similar function in 20180312 notebook.
val_log_termplot <- function(swotlist, mc = TRUE, plot = TRUE, scales = "free") {
  W <- swotlist$W
  S <- swotlist$S
  A <- swotlist$A
  Q <- swotlist$Q
  
  if (mc) {
    Qvec <- apply(Q, 2, median, na.rm = TRUE)
    Q <- matrix(rep(Qvec, nrow(Q)), nrow = nrow(Q), byrow = TRUE)
    swotlist$Q <- Q
  }
  
  fullresids <- log(Q) - manning_qdot(swotlist = swotlist, Avar = "A", 
                                      log = TRUE)
  
  wresids <- fullresids - 2/3 * log(W)
  sresids <- fullresids + 1/2 * log(S)
  aresids <- fullresids + 5/3 * log(A)
  residlist <- list(W = wresids, S = sresids, A = aresids)
  
  datadf <- swotlist[c("W", "S", "A")] %>% 
    map(log) %>% 
    swot_tidy() %>% 
    gather(key = "variable", value = "logval", -time, -loc)
  residdf <- swot_tidy(residlist) %>% 
    gather(key = "variable", value = "p_resid", -time, -loc) %>% 
    left_join(datadf, by = c("variable", "time", "loc")) %>% 
    group_by(variable) %>% 
    mutate(p_resid = p_resid - mean(p_resid, na.rm = TRUE)) %>% 
    ungroup()
  
  if (!plot)
    return(residdf)
  
  plotdf <- residdf %>% 
    group_by(variable) %>% 
    mutate(meanval = mean(logval, na.rm = TRUE),
           # meanresid = mean(p_resid, na.rm = TRUE),
           slope = ifelse(variable == "W", -2/3, 
                          ifelse(variable == "S", 1/2, 
                                 ifelse(variable == "A", 5/3, 
                                        stop("variable not recognized")))),
           intercept = - slope * meanval,
           loc = as.factor(as.integer(loc))) 
  
  out <- ggplot(plotdf, aes(x = logval, y = p_resid, color = loc)) +
    geom_point() + 
    geom_abline(aes(slope = slope, intercept = intercept)) +
    facet_wrap(~variable, scales = scales)
  
  out
}

# 2 versions: one as log-space sd, the other using mean relative diff. 
manning_dqdx <- function(swotlist, Avar = c("Ahat", "A"), 
                         method = c("logsd", "reldiff")) {
  
  method <- match.arg(method)
  Avar <- match.arg(Avar)

  
  rhsmat <- manning_qdot(swotlist = swotlist, Avar = Avar, log = TRUE)  
  if (method == "logsd") {
    sds <- apply(rhsmat, 2, sd, na.rm = TRUE)
    dxerr <- mean(sds) * 100 # times 100 makes it a percent
  } else if (method == "reldiff") {
    rhsmat <- exp(rhsmat)
    diffs <- apply(rhsmat, 2, function(x) abs(diff(x)))
    dxerr <- mean(diffs / rhsmat[-1, ]) / 2 * 100 # times 100 makes it a percent
    dxerr
  }
  
  dxerr
}

val_dqdx <- function(qmat, method = c("logsd", "reldiff")) {
  
  method <- match.arg(method)
  
  if (method == "logsd") {
    sds <- apply(log(qmat), 2, sd, na.rm = TRUE)
    dxerr <- mean(sds) * 100 # times 100 makes it a percent
  } else if (method == "reldiff") {
    diffs <- apply(qmat, 2, function(x) abs(diff(x)))
    dxerr <- mean(diffs / qmat[-1, ]) / 2 * 100
  }
  
  dxerr
}


# Manning ANOVA -----------------------------------------------------------

# RHS of manning equation
manning_qdot <- function(swotlist, Avar = c("Ahat", "A"), 
                         ndot = 1, log = FALSE) {
  
  Avar <- match.arg(Avar)
  
  if (is.null(swotlist[[Avar]])) {
    if (Avar != "Ahat") stop("Missing A data from swotlist")
    
    A0vec <- estA0(swotlist)
    swotlist$Ahat <- swot_A(dAmat = swotlist$dA, A0vec = A0vec)
  }
  
  Wmat <- swotlist$W
  Smat <- swotlist$S
  Amat <- swotlist[[Avar]]
  
  if (length(ndot) > 1 && !is.matrix(ndot)) {
    ndot <- swot_vec2mat(ndot, Wmat)
  }
  
  outmat <- Amat^(5/3) * Wmat^(-2/3) * Smat^(1/2) / ndot
  if (log) 
    outmat <- log(outmat)
  
  outmat
}

# geometric mean (in time and space) of Q and n
manning_qbarn <- function(swotlist, Avar = c("Ahat", "A"), log = FALSE, 
                          meanspace = c("log", "linQ", "linA")) {
  Avar = match.arg(Avar)
  meanspace = match.arg(meanspace)
  
  qdot <- manning_qdot(swotlist = swotlist, Avar = Avar, log = FALSE)
  
  if (meanspace == "log") {
    out <- mean(log(qdot))
  } else if (meanspace == "linQ") {
    out <- log(mean(qdot))
  } else if (meanspace == "linA") {
    out1 <- mean(qdot^(3/5)) ^ (5/3)
    out <- log(out1)
  }

  if (!log) {
    out <- exp(out)
  }
  
  out
}

# Difference in geometric mean n across locations, relative to global geometric mean
manning_ndot <- function(swotlist, Avar = c("Ahat", "A"), log = FALSE, 
                         meanspace = c("log", "linQ", "linA")) {
  Avar = match.arg(Avar)
  meanspace = match.arg(meanspace)
  
  qdot <- manning_qdot(swotlist = swotlist, Avar = Avar, log = FALSE)
  
  if (meanspace == "log") {
    means <- apply(log(qdot), 1, mean, na.rm = TRUE)
  } else if (meanspace == "linQ") {
    means <- log(apply(qdot, 1, mean, na.rm = TRUE))
  } else if (meanspace == "linA") {
    means1 <- apply(qdot, 1, function(x) mean(x^(3/5))) ^ (5/3)
    means <- log(means1)
  }
  
  globmean <- manning_qbarn(swotlist = swotlist, Avar = Avar, log = TRUE,
                            meanspace = meanspace)
  ndots <- means - globmean
  
  if (!log)
    ndots <- exp(ndots)
  
  ndots
}


# Location subset via backwards regression --------------------------------

# Modified from notebook20180315.Rmd.
recondition <- function(modmat) {
  nlocs <- ncol(modmat)
  
  elims <- numeric(0)
  remains <- 1:nlocs
  condnos <- condno(cov(modmat))
  for (i in 1:nlocs) {
    nlocs_remain <- nlocs - i + 1
    if (nlocs_remain <= 2) break
    condnos_i <- map_dbl(1:nlocs_remain, ~condno(cov(modmat[, -.])))
    elim_i <- remains[which.min(condnos_i)]
    
    remains <- remains[-which.min(condnos_i)]
    modmat <- modmat[, -which.min(condnos_i)]
    
    elims <- c(elims, elim_i)
    condnos <- c(condnos, condno(cov(modmat)))
  }
  
  out <- data.frame(elim = c(NA_integer_, elims),
                    n_elim = 0:length(elims),
                    condno = condnos)
  out
}



# Stan stuff --------------------------------------------------------------

prepdata <- function(swotlist) {
  list0 <- with(swotlist, list(
    nt = ncol(W),
    nx = nrow(W),
    dA = dA,
    W = W, 
    S = S, 
    logA0_hat = log(max(dA) - min(dA)),
    logn_hat = -3.5,
    sigma_logn = 0.5))
  
  list1 <- with(c(swotlist, list0), 
                list(sigma_logA0 = 1, 
                     logQ_hat = mean(5/3 * logA0_hat - 2/3 * log(W) + 1/2 * log(S) - logn_hat), 
                     sigma_logQ = 1))
  out <- c(list0, list1)
  out
}

check_rhat <- function(stanfit, ..., plot = TRUE) {
  rhats <- stan_rhat(stanfit, ...)$data
  
  if (plot) {
    plot(rhats$stat, ylab = "Rhat statistic", xlab = "Parameter")
  }
  invisible(rhats)
}

max_rhat <- function(stanfit) {
  rhats <- stan_rhat(stanfit)$data
  max(rhats$stat)
}


# Peek at parameters ------------------------------------------------------

mcman_peek <- function(swotlist, method = c("decomp", "anova"),
                       par = c("mu_q", "sigma_q", "sigma_gprime", 
                               "sigma_nubar", "sigma_man")) {
  
  method <- match.arg(method)
  par <- match.arg(par, several.ok = TRUE)
  meanQ <- apply(swotlist$Q, 2, function(x) mean(log(x)))
  mu_q <- mean(meanQ)
  sigma_q <- sd(meanQ)
  
  charclose <- characterize_closure(swotlist, method = method)
  sigma_gprime <- charclose$dgdx
  sigma_nubar <- charclose$nuhat
  sigma_man <- charclose$err
  
  outlist <- list(mu_q = mu_q, sigma_q = sigma_q, sigma_gprime = sigma_gprime,
                    sigma_nubar = sigma_nubar, sigma_man = sigma_man)
  out <- as.data.frame(outlist[par])
  out
}

mcman_sigma <- function(swotlist, mc = TRUE, na.rm = FALSE) {
  clos <- manning_closure(swotlist, log = TRUE, mc = mc)
  sigma <- mean(apply(clos, 1, sd, na.rm = na.rm), na.rm = na.rm)
  sigma
}
