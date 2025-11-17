# mediaitonBayes: mediation analysis approach leveraging Bayesian framework. 

## Overview

This R package provides a flexible Bayesian mediation analysis framework particularly suitable for zero-inflated data. Both zero-inflated mediators and zero-inflated outcomes are handled by this package, compatible with a wide range of outcome/mediator distributions. This novel technique employs Bayesian models for both the mediator and the outcome, utilizing Markov Chain Monte Carlo (MCMC) algorithms for parameter estimation. This method can be used to conduct mediation analysis with not only zero-inflated mediators but also a wide range of mediator and outcome distributions that can be fitted with GLMs. 

Priors in the R package _brms_ are used by default, please refer to R package _brms_ for more details. Howerver, author self-defined priors can be incorporated. For details, please check the example sectio below.

Author: Jinhong Cui jhcui@uab.edu, cuijinhongqk@gmail.com; Nengjun Yi nyi@uab.edu. 

Maintainer: Jinhong Cui jhcui@uab.edu, cuijinhongqk@gmail.com; 

## Installation
```r
devtools::install_github("jhcuibst/mediationBayes")
```

## Example 

### Load data: 

```r
data(example_data)
```

This data (n = 500) is simulated given the 'datasim' function in our package. Users can define the coefficients and distribution and generate data as they want. This existed example dataset only used as a simple example.

|Name | Meanings of variables in the dataset|
|-----|----------------------------|
|x    |Exposure variable.|
|mnb  |Mediator with zero-infalted negative bionomial (ZINB) distribution.|
|xm   |The interaction between exposure and mediator. User can defined it directly in the model fitting using the same way as thos in R package _glm_.|
|im   |The indicator of zero-inflated mediator, which is used for estimating decomposed Narutal Indirect Effect (NIE).|
|y   |The outcome variable.|

### Conduct mediation analysis 

#### Set up priors
```r
prior.m = c(
          set_prior("student_t(3, 0, 2.5)", class = "b"),  
          set_prior("student_t(3, 0, 2.5)", class = "b", dpar = "hu") )

prior.y = set_prior("student_t(3, 0, 2.5)", class="b")
```

#### Fit Bayesian models
```r
hnb.m <- brm(bf(mnb~ x + age, hu ~ x + age), data = example_data, prior = prior.m,  
             family = hurdle_negbinomial(), chains = 4, iter = 2000)

hnb.y <- brm(y ~  x + mnb + age + im, data = example_data, prior = prior.y,  
            family = bernoulli(), chains = 4, iter = 2000)
```

#### Get outputs of mediation analysis
```r
outmed <- medbayes(hnb.m, hnb.y, mediator="mnb", treat="x", outcome = "y", ind_mediator = "im",  
                   control.value = 0, treat.value = 1)
```

#### Results under Risk Difference (RD) scale
```r
outmed$Results$effects.rd
```
#### Results under Relative Risk (RR) scale
```r
outmed$Results$effects.rr
```
#### Outputs details
The output object contains several components that users can explore to inspect model parameters, extract results, or compute additional quantities of interest.
- **Parameter information:**
   Parameter summaries can be accessed through _outmed$params$..._.
   For example, _outmed$params$model.m_ returns detailed information about the mediator model.

- **Main results:**
   Computed mediation results are stored under _outmed$Results$..._, as shown in the examples above.

- **Posterior predictions:**
   Posterior samples for predicted outcomes (on both the original and transformed scales) are also available.
   For instance, outmed$Results$outcome.linpred returns posterior draws on the original linear predictive scale.
   This object is a 2 × 2 × 2 × iter_numbers × nrow(data) array.
