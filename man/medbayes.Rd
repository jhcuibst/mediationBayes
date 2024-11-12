\name{medbayes}
\Rdversion{1.1.0}
\alias{mediationBayes:medbayes}
\docType{package}
\title{
  Main function for conducting mediation analysis approach.
}

\description{
  This the main function for conducting mediation analysis approach. Mediator model and outcome models should be fitted using package \link{brms} by users. Examples for fitting theses models are listed in the 'usage' section below. Both Zero-Inflated mediators and Zero-Inflated outcomes are accomodated by this function. Causal effects decompositions are provided. All of the GMLs are accommodated by this function by default.

  Priors can be defined by user as they need. If no priors are specified, then default priors are utilized. For details of prior distributions, please refer to R package \link{brms}.
}

\arguments{
*medbayes <- function(model.m = NULL, model.y = NULL,
                     treat = "treatment", mediator = "mediator", ind_mediator = NULL, outcome = "outcome",
                     control.value = 0, treat.value = 1,
                     logM = FALSE )*

\item{model.m}{
  mediator model: can be all the standard families defined in \code{\link{brms}}.
}

\item{model.y}{
  Outcome model: can be all the standard families defined in \code{\link{brms}}.
}

\item{treat}{
   the name of the treatment or exposure variable, need quotation to specify it.
}

\item{mediator}{
   the name of the mediator variable, need quotation to specify it.
}

\item{ind_mediator}{
   if the mediator variable is zero-inflated, this item could be specified. By specifying this item, outcomes will provide decomposed mediation analysis. Otherwise, only overall mediation analysis without decomposition is provided.
}

\item{outcome}{
   the name of the outcome variable, need quotation to specify it.
}

\item{control.value}{
  the controlled value of treatment or exposure that user interested.
}

\item{treat.value}{
  the treatment value of treatment or exposure that user interested.
}

}
\details{
  This function fits Zero-Inflated mediators via hurdle models and fits zero-infalted outcomes via either hurdle models or zero-inflated models.
  The ZI mediator can be seen as a two-sequential mediator with natural indirect effect (NIE) partitioned into two different pathways one entire indirect effect transmitted along \eqn{M}, denoted by the solid
line pathway in Figure 1(A) (we defined it as \eqn{NIE_1}) and the other indirect effect transmitted through \eqn{I_{M>0}}
only, represented by the solid line pathway in Figure 1(B) (we defined it as \eqn{NIE2}).
}

\examples{
# Loade data
data(example_data)

# set up priors
prior.m = c(
  set_prior("student_t(3, 0, 2.5)", class = "b"),
  set_prior("student_t(3, 0, 2.5)", class = "b", dpar = "hu")
)

prior.y = set_prior("student_t(3, 0, 2.5)", class="b")

# fit models
hnb.m <- brm(bf(mnb~ x + c, hu ~ x + c), data = data(example_data), prior = prior.m,
            family = hurdle_negbinomial(), chains = 4, iter = 2000)
hnb.y <- brm(y ~  x + m + c + im, data = data(example_data), prior = prior.y,
             family = bernoulli(), chains = 4, iter = 2000)

outmed <- medbayes(hnb.m, hnb.y, mediator="mnb", treat="x", outcome = "y", ind_mediator = "im",
                  control.value = 0, treat.value = 1)

# For binary outcomes, we provide results under both Relative Risk (RR) and Risk Difference (RD) scales

## to get restuls under RR  scale:
outmed$effects.rr

## to get restuls under RD scale:
outmed$effects.rd
}

\author{
  Jinhong Cui, jhcui@uab.edu
}
