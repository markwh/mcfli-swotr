---
title: "Untitled"
author: "Mark Hagemann"
date: "February 15, 2018"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

Make some synthetic data to try and understand A0 calculation. 

If I have a system that varies in time and space governed by Manning's equation. Distributions shouldn't matter. 

Simulate

- Q_t
- A_0i
- delta_it
- x_it


```{r}
nt <- 10000
nx <- 2
qt <- rnorm(nt, 5, 0.5)
qit <- matrix(rep(qt, nx), nrow = nt)
xit <- matrix(rnorm(nt * nx, 3, 0.3), nrow = nt)
ait <- qit / xit
# a0 <- rnorm(nx, 2, 0.2)
a0 <- c(1, 2)
a0it <- matrix(rep(a0, nt), nrow = nt, byrow = TRUE)
dait <- ait - a0

yvec <- (xit * -dait) %*% c(1, -1)

Omega <- matrix(c(1, 0, 0, -1), nr = 2)

desgn <- xit %*% Omega


moddf <- data.frame(x1 = desgn[, 1], x2 = desgn[, 2], y = yvec)

modlm <- lm(y ~ 0 + x1 + x2, moddf)

summary(modlm)
a0

pairs(moddf)

foo <- desgn %*% a0
bar <- desgn %*% c(1.5, 1.5)

head(foo)
head(bar)
head(yvec)

plot(foo, yvec)
plot(bar, yvec)
```

Why isn't my equality holding?

```{r}
mylhs <- xit %*% (a0 * c(1, -1))
myrhs <- (-dait * xit) %*% c(1, -1)

head(mylhs)
head(myrhs)

plot(mylhs, myrhs)
```

OK, try this (working upstream)

```{r}
summary((ait  * xit) %*% c(1, -1))
```

That checks out

