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