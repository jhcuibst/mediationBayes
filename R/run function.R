
library(tidyverse)
library(brms)

# source("d:/01 dissertation research/Bayesian, Zero-Inflated mediation/02 hurdle model/simulationg_Feb2024/binary_0201/datasim.R")

source("datasim.R")

# prior.m = priors <- c(
#   set_prior("student_t(3, 0, 10)", class = "b"),                      # Coefficients for the count model
#   set_prior("student_t(3, 0, 2.5)", class = "Intercept"),             # Intercept for the count model
#   set_prior("student_t(3, 0, 2.5)", class = "b", dpar = "hu"),        # Coefficients for the hurdle part
#   set_prior("student_t(3, 0, 10)", class = "Intercept", dpar = "hu"), # Intercept for the hurdle part
#   set_prior("gamma(0.01, 0.01)", class = "shape")                     # Shape parameter for the NB distribution
# )
# 
# prior.y = set_prior("student_t(3, 0, 2.5)", class="b") + 
#           set_prior("student_t(3, 0, 10)", class="Intercept")

# fit models
hnb.m <- brm(bf(mnb~ x + age, hu ~ x + age), data = data1, family = hurdle_negbinomial(), chains = 4, iter = 2000) # prior = prior.m,  
hnb.y <- brm(y ~  x + mnb + age + im, data = data1, family = bernoulli(), chains = 4, iter = 2000) # prior = prior.y,

source("brms_mediation_jc.R")

out.edit = medbayes(hnb.m, hnb.y, mediator="mnb", treat="x", outcome = "y", ind_mediator = "im",
                    control.value = 0, treat.value = 1)

rst.rr <- out.edit$effects.rr
rst.rr

