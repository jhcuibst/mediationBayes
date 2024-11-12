#' Example Dataset
#'
#' This dataset is provided as an example for testing the package functions.
#'
#' @format A data frame with 5 rows and 2 columns:
#' \describe{
#'   \item{x}{The treatment / exposure variable}
#'   \item{x.d}{A matrix containing exposure and all covariates variables}
#'   \item{coef.z}{Coefficients for defining zero-inflation in ZI mediator}
#'   \item{coef.d}{Coefficients for defining non-zero inflated class distribution}
#'   \item{mdist}{The distribution of mediator}
#'   \item{theta}{The dispersion parameter for mediator}
#'   \item{ydist}{The distribution of outcome}
#'   \item{coef.y}{Coefficients of outcome model}
#' }
#' @usage data(example_data)
#' @examples
#' data(example_data)
#' head(example_data)
"example_data"
