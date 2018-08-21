# main.R
# version-testing bamr runs.

# CHANGE THIS

ref <- "consolidate"
case <- "Cumberland"
seed <- 6346324389

# Install bamr

devtools::install_github("markwh/bamr", ref = ref, local = FALSE)

# Select a case to test

load("reachdata.RData")

swotlist <- reachdata[[case]] %>% 
  swotr::swot_sset(keeptimes = round(seq(1, ncol(swotlist$W), length.out = 15)),
                   keeplocs = 1:4) %>% 
  swotr::swot_purge_nas()

qobs <- apply(swotlist$Q, 2, median)
  
testbd <- swotlist %>% 
  swotr::swot_bamdata()

# Run!

bamest1 <- bamr::bam_estimate(data = testbd, variant = "manning", 
                              meas_error = FALSE, reparam = FALSE,
                              cores = 2, chains = 2, seed = seed)

bamest2 <- bamr::bam_estimate(data = testbd, variant = "manning_amhg",
                              meas_error = FALSE, reparam = FALSE,
                              cores = 2, chains = 2, seed = seed)

bamest3 <- bamr::bam_estimate(data = testbd, variant = "amhg",
                              meas_error = FALSE, reparam = FALSE,
                              cores = 2, chains = 2, seed = seed)