
testcase <- within(reachdata$Seine, {x = x / 10000})
# testcase <- within(uscase, {x = x / 10000; dA = rezero_dA(dA, "median")}) %>%
# swot_sset(keeptimes = 1:25 * 3)
inputs <- testcase %>% gibbs_inputs()
# gibbs_inputs(Qhat = 120)


samps1 <- gibbs_sample(inputs = inputs, chains = 3, iter = 3000,
                       thin = 1, serial = FALSE)

samps2 <- gibbs_sample(inputs = inputs, chains = 3, iter = 3000,
                       gamma = TRUE, nu = FALSE,
                       serial = FALSE)

samps3 <- gibbs_sample(inputs = inputs, chains = 3, iter = 3000,
                       gamma = FALSE, nu = TRUE,
                       serial = FALSE)

samps4 <- gibbs_sample(inputs = inputs, chains = 3, iter = 3000,
                       gamma = TRUE, nu = TRUE,
                       serial = FALSE)

samps5 <- gibbs_sample(inputs = inputs, chains = 3, init_method = "low",
                       iter = 1000,
                       gamma = TRUE, nu = TRUE, thin = 1,
                       serial = FALSE)
bam1 <- bam_estimate(swot_bamdata(testcase, Qhat = 120), variant = "manning")


stan_trace(gibbs_stanfit(samps1), pars = "A0", inc_warmup = TRUE)
stan_trace(gibbs_stanfit(samps2), pars = "A0", inc_warmup = TRUE)
stan_trace(gibbs_stanfit(samps3), pars = "A0", inc_warmup = TRUE)
stan_trace(gibbs_stanfit(samps4), pars = "A0", inc_warmup = TRUE)
