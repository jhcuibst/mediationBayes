# mediaitonBayes: mediation analysis approach leveraging Bayesian framework. 

## Overview

This R package provides a flexible Bayesian mediation analysis framework particularly suitable for zero-inflated data. Both zero-inflated mediators and zero-inflated outcomes are handled by this package, compatible with a wide range of outcome/mediator distributions. This novel technique employs Bayesian models for both the mediator and the outcome, utilizing Markov Chain Monte Carlo (MCMC) algorithms for parameter estimation. This method can be used to conduct mediation analysis with not only zero-inflated mediators but also a wide range of mediator and outcome distributions that can be fitted with GLMs. 

Priors in the R package _brms_ are used by default, please refer to R package _brms_ for more details. Howerver, author self-defined priors can be incorporated. For example, in the mediation analysis with ZI mediator fitted with Hurdle model, prior can be set by:

*prior.m = c(  
  set_prior("student_t(3, 0, 2.5)", class = "b"), # Coefficients for the count model  
  set_prior("student_t(3, 0, 2.5)", class = "b", dpar = "hu") # Coefficients for the hurdle part  
)*

Author: Jinhong Cui jhcui@uab.edu, cuijinhongqk@gmail.com; Nengjun Yi nyi@uab.edu. 

Maintainer: Jinhong Cui jhcui@uab.edu, cuijinhongqk@gmail.com; 

## Installation
devtools::install_github("jhcuibst/mediationBayes")

# Example 

### load data using: data(example_data)

This data (n = 500) is simulated given the 'datasim' function in our package. Users can define the coefficients and distribution and generate data as they want. This existed example dataset only used as a simple example.

\arguments{
     \item{x}{
     exposure
     }

  \item{mnb}{
    Mediator with zero-infalted negative bionomial (ZINB) distribution.
  }

  \item{xm}{
    The interaction between exposure and mediator. User can defined it directly in the model fitting.
  }

  \item{im}{
    The indicator of zerl-inflated mediator, which is used for estimating decomposed Narutal Indirect Effect (NIE).
  }

  \item{y}{
    The outcome variable.
  }

}

\usage{
data(example_data)
# Set up priors
prior.m = c(
  set_prior("student_t(3, 0, 2.5)", class = "b"),
  set_prior("student_t(3, 0, 2.5)", class = "b", dpar = "hu") )

prior.y = set_prior("student_t(3, 0, 2.5)", class="b")

# Fit Bayesian models
hnb.m <- brm(bf(mnb~ x + age, hu ~ x + age), data = example_data, prior = prior.m,
             family = hurdle_negbinomial(), chains = 4, iter = 2000)

hnb.y <- brm(y ~  x + mnb + age + im, data = example_data, prior = prior.y,
            family = bernoulli(), chains = 4, iter = 2000)

# Get outputs of mediation analysis
outmed = medbayes(hnb.m, hnb.y, mediator="mnb", treat="x", outcome = "y", ind_mediator = "im",
          control.value = 0, treat.value = 1)

# Results under Risk Difference (RD) scale
outmed$effects.rd
# Results under Relative Risk (RR) scale
outmed$effects.rr
}
