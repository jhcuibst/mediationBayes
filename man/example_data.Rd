\name{datasim}
\Rdversion{1.1.0}
\alias{mediationBayes:medbayes}
\docType{package}
\title{
  The example data used for testing function.
}

\description{
  This data is simulated given the 'datasim' function in our package. Users can define the coefficients and distribution and generate data as they want. This existed example dataset only used for testing.
}

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

\author{
  Jinhong Cui, jhcui@uab.edu
}
