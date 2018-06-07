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
  
  nc <- max(swotdf$time)
  nr <- max(swotdf$loc)
  
  
  swotdf <- arrange(swotdf, time, loc)
  out <- map(swotdf[matnames], 
             ~matrix(., nrow = nr, ncol = nc, byrow = FALSE))
  out
}


