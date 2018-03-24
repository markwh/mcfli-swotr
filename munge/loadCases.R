# loadCases.R
# Mark Hagemann
# 3/20/2018
# Modified from various notebooks, including 20180316.Rmd

pep1cases <- names(Pepsi_v2)
pep1ncs <- sprintf("../swotData/data/NC_files/%s.nc", pep1cases)

nclists <- map(pep1ncs, nc_list) %>% 
  setNames(pep1cases)

nclists$Platte$River_Info.gdrch <- 1:12 # get rid of 13-14.

gdrchs <- map(nclists, ~.$River_Info.gdrch)

keepvars <- c("W", "S", "A", "Q")
reachparts <- paste0("Reach_Timeseries.", keepvars)

zeroToNA <- function(mat) {
  zeros <- mat <= 0
  mat[zeros] <- NA
  mat
}

reachdata0 <- nclists %>% 
  map(~.[reachparts]) %>% 
  map(setNames, keepvars) %>% 
  map2(gdrchs, function(x, y) map(x, function(z) z[y, ])) %>% 
  map(~map(., zeroToNA))

# Get rid of problematic times (spinup issue?)
reachdata0$Severn$W <- reachdata0$Severn$W[, -1]
reachdata0$Severn$A <- reachdata0$Severn$A[, -1]
reachdata0$Severn$S <- reachdata0$Severn$S[, -1]
reachdata0$Severn$Q <- reachdata0$Severn$Q[, -1]

reachdata0$Platte$W <- reachdata0$Platte$W[, -1]
reachdata0$Platte$A <- reachdata0$Platte$A[, -1]
reachdata0$Platte$S <- reachdata0$Platte$S[, -1]
reachdata0$Platte$Q <- reachdata0$Platte$Q[, -1]

reachncols <- map(reachdata0, ~ncol(.$A))

reachA0 <- reachdata0 %>% 
  map(~.$A) %>% 
  map(apply, 1, median) %>% 
  map2(reachncols, function(x, y) matrix(rep(x, y), ncol = y, byrow = FALSE))

reachdA <- reachdata0 %>% 
  map(~.$A) %>% 
  map2(reachA0, function(x, y) x - y)

reachdata <- map2(reachdata0, reachdA, function(x, y) {x$dA <- y; x})

Qmats <- nclists %>% 
  map(~.$Reach_Timeseries.Q[.$River_Info.gdrch, ])

cache("reachdata")
cache("Qmats")