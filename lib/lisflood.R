# lisflood.R
# Mark Hagemann
# 3/21/2018
# functions to assist reading and writing lisflood inputs/outputs.

#' @importFrom readr read_fwf
read_lisProfile <- function(file) {
  infile <- readr::read_fwf(file, 
                     fwf_empty(file, 
                               c("ChanX", "ChanY", "Chainage", "Width", 
                                 "Mannings", "Slope", "BankZ", "BedElev", 
                                 "WaterElev", "WaterDepth", "Flow"), skip = 2), 
                     skip = 2, col_types = "ddddddddddd")
  infile
}

# function to read lisflood output into swot-like data format
lis_profiles <- function(resultdir) {
  pnums <- list.files(resultdir, 
                      pattern = "\\.profile$", full.names = FALSE) %>% 
    stringr::str_extract("[0-9]{4}")
  
  profiles <- list.files(resultdir, 
                         pattern = "\\.profile$", full.names = TRUE) %>% 
    map(read_lisProfile) %>% 
    setNames(pnums) %>% 
    bind_rows(.id = "profile") %>% 
    mutate(time = as.numeric(profile)) %>% 
    arrange(time, ChanX) %>% 
    group_by(time) %>% 
    mutate(Slope = - c(NA, diff(WaterElev)) / c(NA, diff(ChanX))) %>% 
    ungroup()
  profiles
}

lis_reaches <- function(resultdir, reachlen_m = 10000) {
  profiles <- lis_profiles(resultdir)
  reaches <- profiles %>% 
    mutate(loc = floor(ChanX / reachlen_m) + 1) %>% 
    group_by(loc, time) %>% 
    summarize(W = median(Width), 
              S = median(Slope, na.rm = TRUE), 
              H = median(WaterElev),
              Q = median(Flow)) %>% 
    group_by(loc) %>% 
    mutate(dA = calcdA_vec(w = W, h = H, zero = "minimum")) %>% 
    ungroup()
  reaches
}
