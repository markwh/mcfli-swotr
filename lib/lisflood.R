# lisflood.R
# Mark Hagemann
# 3/21/2018
# functions to assist reading and writing lisflood inputs/outputs.

#' @importFrom readr read_fwf
read_lisProfile <- function(file) {
  # infile <- readr::read_fwf(file, 
  #                    fwf_empty(file, 
  #                              c("ChanX", "ChanY", "Chainage", "Width", 
  #                                "Mannings", "Slope", "BankZ", "BedElev", 
  #                                "WaterElev", "WaterDepth", "Flow"), skip = 2), 
  #                    skip = 2, col_types = "ddddddddddd")
  
  infile <- read.table(file, header = TRUE, skip = 1, sep = "", stringsAsFactors = FALSE)
  infile
}

read_lispar <- function(file) {
  
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
    mutate(dH = - c(NA, diff(WaterElev)),
           dx = c(NA, diff(ChanX)), 
           Slope = dH / dx) %>% 
    ungroup()
  profiles
}

avgSlope <- function(dz, dx, method = c("median_s", "mean_s", "s_mean", "s_median")) {
  method = match.arg(method)
  
  if (method == "median_s") {
    out <- median(dz / dx, na.rm = TRUE)
  } else if (method == "mean_s") {
    out <- mean(dz / dx, na.rm = TRUE)
  } else if (method == "s_mean") {
    out <- mean(dz, na.rm = TRUE) / mean(dx, na.rm = TRUE)
  } else if (method == "s_median") {
    out <- median(dz, na.rm = TRUE) / median(dx, na.rm = TRUE)
  }
  
  out
}

lis_reaches <- function(lisprofiles, reachlen_m = 10000, 
                        slope_method = c("median_s", "mean_s", "s_mean", "s_median"),
                        agfun = NULL) {
  slope_method <- match.arg(slope_method)
  if (is.null(agfun)) {
    agfun <- median
  }
  
  reaches <- lisprofiles %>% 
    mutate(loc = floor(ChanX / reachlen_m) + 1) %>% 
    group_by(loc, time) %>% 
    summarize(W = agfun(Width), 
              S = avgSlope(dz = dH, dx = dx, method = slope_method),
              # S = median(Slope, na.rm = TRUE),
              H = agfun(WaterElev),
              D = agfun(WaterDepth),
              Q = agfun(Flow),
              x = median(ChanX)) %>% 
    group_by(loc) %>% 
    mutate(dA = calcdA_vec(w = W, h = H, zero = "minimum")) %>% 
    ungroup()
  reaches
}



# Boundary condition generation -------------------------------------------

# Modified from work in 20180313 notebook.

lis_happroxfun <- function(qinfile, houtfile, sep = "\t") {
  qin <- read.csv(qinfile, sep = sep)[["q_m3s"]]
  hout <- read.csv(houtfile, sep = sep)[["h_m"]]
  out <- approxfun(x = qin, y = hout)
  out
}

lis_bdry <- function(q_impose, h_impose, time_sustain, time_change) {
  q_write <- rep(q_impose, each = 2)
  h_write <- rep(h_impose, each = 2)
  t_write <- cumsum(c(0, rep(c(time_sustain, time_change), 
                             length.out = length(h_write)))) %>% 
    `[`(1:length(h_write))
  
  out_q <- data.frame(q_m3s = q_write, t_s = t_write)
  out_h <- data.frame(h_m = h_write, t_s = t_write)
  
  out <- list(q = out_q, h = out_h)
  
  out
}

lis_write_bdry <- function(bdrylist, file, header = "") {
  write_lines(header, file, append = FALSE)
  write_lines("upstream1", file, append = TRUE)
  write_lines(sprintf("%s\t\t seconds", nrow(bdrylist$q)), file, append = TRUE)
  write_tsv(bdrylist$q, file, append = TRUE)
  write_lines("", file, append = TRUE)
  write_lines("downstream1", file, append = TRUE)
  write_lines(sprintf("%s\t\t seconds", nrow(bdrylist$h)), file, append = TRUE)
  write_tsv(bdrylist$h, file, append = TRUE)
}


# Tack metroman function here for now -------------------------------------

read_metroman <- function(file) {
  filechar <- readChar(file, file.info(file)$size)
  filestr <- read_file(file) %>% 
    stringr::str_split(pattern = "\\r\\n") %>% 
    `[[`(1)
  
  txtlines <- stringr::str_which(filestr, "[a-zA-Z]")
  numstarts <- txtlines + 1
  numends <- c(txtlines[-1] - 1, length(filestr))
  
  numvals <- map2(numstarts, numends, function(x, y) filestr[x:y]) %>% 
    map(stringr::str_trim) %>% 
    map(stringr::str_split, pattern = " +") %>% 
    map(~map(., as.numeric)) %>% 
    map(as.data.frame) %>% 
    map(unname) %>% 
    map(as.matrix) %>% 
    map(t)
  
  timeline <- grep("^time", filestr, ignore.case = TRUE)
  tpart <- which(txtlines == timeline)
  wline <- grep("^width", filestr, ignore.case = TRUE)
  wpart <- which(txtlines == wline)
  sline <- grep("^slope", filestr, ignore.case = TRUE)
  spart <- which(txtlines == sline)
  hline <- min(grep("^height", filestr, ignore.case = TRUE))
  hpart <- which(txtlines == hline)
  
  out <- numvals[c(tpart, wpart, spart, hpart)] %>% 
    setNames(c("days", "W", "S", "H"))
  bamda <- calcdA_mat(out$W, out$H)
  out$days <- as.vector(out$days)
  out$S <- out$S * 1e-5 # previously in cm/km. Now unitless. 
  out$dA <- bamda
  
  out
}
