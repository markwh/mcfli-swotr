# anova.R
# Mark Hagemann
# 3/20/2018
# Funcitons to assist ANOVA on Manning's n. 

swot_n_anova_lm <- function(swotlist) {
  anovalm <- swotlist %>% 
    swot_tidy() %>% 
    mutate(D = A / W, N = W^(-2/3) * A^(5/3) * S^(1/2) / Q, 
           locfac = as.factor(loc), logn = log(N), logd = log(D)) %>% 
    lm(data = ., formula = logn ~ locfac * logd)
  anovalm
}

swot_n_anova <- function(swotlist) {
  out <- broom::tidy(anova(swot_n_anova_lm(swotlist)))
  out
}