# loadCases.R
# Mark Hagemann
# 3/20/2018
# Modified from various notebooks, including 20180316.Rmd

load("../SWOT/cache/nc_r.RData")
pep1cases <- names(nc_r)
pep1ncs <- sprintf("../swotData/data/NC_files/%s.nc", pep1cases)

nclists <- map(pep1ncs, nc_list) %>% 
  setNames(pep1cases)

nclists$Platte$River_Info.gdrch <- 1:12 # get rid of 13-14.

# Convert x and t to matrices. 
xtToMat <- function(lst) {
  x0 <- lst$x
  x <- (x0[-1] + x0[-length(x0)]) / 2
  lst$x <- swot_vec2mat(x ,lst$W)
  lst$t <- swot_vec2mat(lst$t, lst$W)
  lst
}

gdrchs <- map(nclists, ~.$River_Info.gdrch)

keepvars <- c("W", "S", "A", "Q", "t")

reachparts <- paste0("Reach_Timeseries.", keepvars) %>% 
  c("River_Info.rch_bnd")

zeroToNA <- function(mat) {
  zeros <- mat <= 0
  mat[zeros] <- NA
  mat
}

reachdata0 <- nclists %>% 
  map(~.[reachparts]) %>% 
  map(setNames, c(keepvars, "x")) %>% 
  map(xtToMat) %>% 
  map2(gdrchs, function(x, y) map(x, function(z) z[y, ])) %>% 
  map(~map(., zeroToNA))

# Get rid of problematic times (spinup issue?)
reachdata0$Severn <- swot_sset(reachdata0$Severn, keeptimes = -1)
reachdata0$Platte <- swot_sset(reachdata0$Platte, keeptimes = -1)

reachncols <- map(reachdata0, ~ncol(.$A))

reachA0 <- reachdata0 %>% 
  map(~.$A) %>% 
  map(apply, 1, median) %>% 
  map2(reachncols, function(x, y) matrix(rep(x, y), ncol = y, byrow = FALSE))

reachdA <- reachdata0 %>% 
  map(~.$A) %>% 
  map2(reachA0, function(x, y) x - y)

reachdata1 <- map2(reachdata0, reachdA, function(x, y) {x$dA <- y; x})

Qmats <- nclists %>% 
  map(~.$Reach_Timeseries.Q[.$River_Info.gdrch, ])
Qhats <- nclists %>% 
  map(~.$River_Info.QWBM[1])

reachdata <- map2(reachdata1, Qhats, 
                  function(x, y) {attr(x, "QWBM") <- y; x})

cache("reachdata")
cache("Qmats")
cache("Qhats")


# Cross-section data ------------------------------------------------------

xsparts <- c("XS_Timeseries.t", "XS_Timeseries.Z", 
"XS_Timeseries.xs_rch", "XS_Timeseries.X", "XS_Timeseries.W", 
"XS_Timeseries.Q", "XS_Timeseries.H", "XS_Timeseries.A", "XS_Timeseries.P", 
"XS_Timeseries.n")

xsnames <- gsub("XS_Timeseries.", "", xsparts)

xsdata <- nclists %>% 
  map(~.[xsparts]) %>% 
  map(~setNames(., xsnames))


# lisflood cases

is_ss <- function(...) {
  arglist <- list(...)
  ssvars <- map(arglist, function(x) c(1, diff(x))) %>% 
    map(function(x) x == 0) %>% 
    Reduce(`*`, x = .)
  out <- as.logical(ssvars)
  out
}

# sscase <- lis_profiles("lisflood/toy_1/results_const_v2/") %>% 
#   lis_reaches(slope_method = "s_median", agfun = median) %>% 
#   filter(loc != 4) %>%
#   arrange(loc, time) %>% 
#   group_by(loc) %>% 
#   mutate(steady = is_ss(W, S, Q, H, D)) %>% 
#   group_by(time) %>% 
#   filter(sum(steady) == 3) %>% 
#   group_by(Q) %>% 
#   filter(time == time[1]) %>% 
#   mutate(A = W * D) %>% 
#   select(-steady) %>% 
#   swot_untidy()
attr(sscase, "QWBM") <- 120 # Fake QWBM to ease bam etc estimation

# uscase <- lis_profiles("lisflood/toy_1/results_simple_v2/") %>% 
#   lis_reaches(slope_method = "s_mean", agfun = mean) %>% 
#   filter(loc != 4) %>%
#   mutate(A = W * D) %>% 
#   swot_untidy() %>% 
#   swot_sset(keeptimes = -1:-20)
attr(uscase, "QWBM") <- 120 # Fake QWBM to ease bam etc estimation

cache("sscase")
cache("uscase")