
# Sim function ----
dat.sim <- function(
    x.d,               # Predictors for log link, if missing, then "x "is generated based on "px" with level "0", "1" etc.
    x.z,               # Predictors for logit link
    coef.d = NULL,     # coefficients of log model, including intercept
    coef.z = NULL,     # coefficients of logit model, including intercept
    dist = c("ZINB", "ZIP", "NB", "Poi", "ZILN"),
    sd = NULL,
    theta = NULL     # shape parameter of NB distribution, Var = u + (u^2)/theta
)
{
  # set.seed(seed)
  if (!(ncol(x.d) == length(coef.d))){
    stop("Number of coefficients and predictors not match in log link. Coefficients should include the intercept.")
  }

  # ------------------------------------------------------------------------------------------------
  ## Simulate Covariates
  # ------------------------------------------------------------------------------------------------
  if (missing(x.z)) x.z <- x.d
  if (missing(coef.z)) coef.z = NULL

  coef.d <- as.matrix(coef.d)
  coef.z <- as.matrix(coef.z)

  # ------------------------------------------------------------------------------------------------
  ## Simulate Mediator
  # ------------------------------------------------------------------------------------------------

  etaz <- x.z %*% as.matrix(coef.z)
  m.normal <- rnorm(n, -etaz, 1.6)
  p.zero = exp(etaz)/ (1+exp(etaz))
  p.zero = round(mean(p.zero),2)
  quantiles <- quantile(m.normal, p.zero)
  m.z <-  as.numeric( factor(cut(m.normal, breaks = c(-Inf, quantiles, Inf))) ) - 1

  etad <- x.d %*% coef.d
  if (dist %in% c("ZINB", "NB") )  m.nb <- rnbinom(n, mu = exp(etad), size = theta)
  if (dist %in% c("ZIP", "Poi") )  m.nb <- rpois(n, lambda =  exp(etad))
  if (dist %in% c("ZILN") )        m.nb <- exp(rnorm(n, mean = etad, sd = sd))
  # if (dist %in% c("ZILN") )        m.nb <- rlnorm(n, meanlog = etad, sdlog = exp(sd))

  zi.y = ifelse(m.z == 0, 0, m.nb)

  df.y <- data.frame(x = x, m = m, c = c, zi.y = zi.y)
  df.y$ind.y = ifelse(df.y$zi.y == 0, 0, 1)

  return(df.y)
}


