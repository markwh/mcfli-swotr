# anova.R
# Mark Hagemann
# 3/20/2018
# Funcitons to assist ANOVA on Manning's n. 

val_anova_lm <- function(swotlist) {
  anovalm <- swotlist %>% 
    swot_tidy() %>% 
    mutate(D = A / W, N = W^(-2/3) * A^(5/3) * S^(1/2) / Q, 
           locfac = as.factor(loc), logn = log(N), logd = log(D)) %>% 
    lm(data = ., formula = logn ~ locfac * logd)
  anovalm
}

val_anova <- function(swotlist) {
  out <- broom::tidy(anova(val_anova_lm(swotlist))) %>% 
    mutate(pctTotVar = sumsq / sum(sumsq) * 100,
           term = as.character(term),
           term  = plyr::mapvalues(term, 
               from = c("locfac", "locfac:logd", "logd"),
               to = c("location", "location:depth", "depth")))
  out
}

val_nd_plot <- function(swotlist, plot = TRUE) {
  
  plotdf <- swotlist %>% 
    swot_tidy() %>% 
    mutate(logn = -2/3 * log(W) + 5/3 * log(A) + 1/2 * log(S) - log(Q), 
         logd = log(A / W),
         loc = as.factor(loc))
  
  if (!plot) return(plotdf)
  
  out <- plotdf %>% 
    ggplot(aes(x = logd, y = logn)) +
    geom_point(aes(color = loc), size = 0.2) + 
    stat_smooth(aes(group = loc, color = loc), method = "lm", se = FALSE)
  out
}