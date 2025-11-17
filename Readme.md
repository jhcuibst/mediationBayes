# mediaitonBayes: mediation analysis approach leveraging Bayesian framework. 

## Overview

This R package provides a flexible Bayesian mediation analysis framework particularly suitable for zero-inflated data. Both zero-inflated mediators and zero-inflated outcomes are handled by this package, compatible with a wide range of outcome/mediator distributions. This novel technique employs Bayesian models for both the mediator and the outcome, utilizing Markov Chain Monte Carlo (MCMC) algorithms for parameter estimation. This method can be used to conduct mediation analysis with not only zero-inflated mediators but also a wide range of mediator and outcome distributions that can be fitted with GLMs. 

For Zero-Inflated meditor, users can define the indicator of mediator to get decomposed mediation analysis results, or users will get overall mediaiton results if no indicator of mediator is defined. Please check the results from example below for detailed explanations.

Our package handles both continuous and binary outcomes, with zero-inflated mediator. For coninuous outcomes, our approach summarizes mediation analysis results under difference scale by default. For binary outcomes, we provide results under both risk difference (RD) scale and risk ratio (RR) scale. 

Priors in the R package _brms_ are used by default, please refer to R package _brms_ for more details. Howerver, author self-defined priors can be incorporated. For details, please check the example section below.

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
#### Explanation of each item from output

- **Indirect_treated**: This corresponds to the Natural Indirect Effect (NIE) as described in our paper. It aligns with the Total Indirect Effect (TIE) based on Robins’ (2003)<sup>1</sup> definition. Conceptually, it quantifies the average change in the outcome that would occur if the mediator were shifted from the value it would take under the treatment condition to the value it would take under the control condition, while holding the exposure fixed (to control level).
  
- **Indirect_Indicator**: If user defined the parameter **ind_mediator**, that means users are requiring estimating decomposed indirect effects. This is part of a decomposed indirect effect, which signifies the average change in outcome, if the non-zero porpotion of mediator were shifted from the proportion it would take under treatment condition to the value it would take under the control condition, while holding the exposure fixed (to control level).
  
- **Indirect_total**: If no decomposition was required, the **Indirect_treated** should be the same as **Indirect_treated** value. If decomposition was required, then the **Indrect_total** = **Indirect_treated** + **Indirect_Indicator**, under difference scale, or **Indrect_total** = **Indirect_treated** * **Indirect_Indicator** under ratio scale.
  
- **Direct**: This corresponds to the Natural Direct Effect (NDE) as described in our paper. It aligns with the Pure Direct Effect (PDE) based on Robins’ (2003)<sup>1</sup> definition. Conceptually, it represents the average change in the outcome if the exposure were shifted from treatment to control, while holding the mediator fixed at the value it would be expected to take under the control condition.

- **Total effect**: Under difference scale, total effect equals **NDE + Indirect_total**, if NIE decomposition was required. Under risk ratio scale, total effect equals **NDE * Indirect_total**.
  
- **Prop.Med**: Proportion of Mediation. This represents the proportion of total effect that is mediated through mediator. Under difference scale, **Prop.Med = Indirect_total / Total Effect**. While under ratio scale, **Prop.Med = NDE * (Indirect_total-1) / (Total Effect-1)**.<sup>2</sup>
  

#### Outputs details

The output object contains several components that users can explore to inspect model parameters, extract results, or compute additional quantities of interest.
- **Parameter information:**
   Parameter summaries can be accessed through _outmed$params$..._.
   For example, _outmed$params$model.m_ returns detailed information about the mediator model.

- **Main results:**
   Computed mediation results are stored under _outmed$Results$..._, as shown in the examples above.

- **Posterior predictions:**
   Posterior samples for predicted outcomes (on both the original and transformed scales) are also available.
   For instance, _outmed$Results$outcome.linpred_ returns posterior draws on the original linear predictive scale.
   This object is a 2 × 2 × 2 × iter_numbers × nrow(data) array.


_[1] Robins JM. Semantics of causal DAG models and the identification of direct and indirect effects. Highly structured stochastic systems. 2003 May 1:70-82._  
_[2] VanderWeele T. Explanation in causal inference: methods for mediation and interaction. Oxford University Press; 2015 Feb 13._
