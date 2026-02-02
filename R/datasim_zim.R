
# Sim function ----
dat.sim <- function(
    n.ind = 1,         # indicator of number of clusters (fix to 1 by default at present)
    n.measure = n,     # number of observations for each group (total sample size given cluster =1)
    px = c(0.5, 0.5),  # Vector of Proportion(s) for (each) level of exposure variable
    ## ("0.5,0.5 "by default, can be multiple numbers for ordinal variable)
    x.d,               # exposure variable, if missing, then "x.d "is generated based on "px" with level "0", "1" etc.
    x.z,               # Covariates (should be a vector or matrix)
    coef.d = 0,        # coefficients of log model, including intercept
    coef.z = NULL,     # coefficients of logit model, including intercept
    mdist = c("ZINB", "ZIP", "NB", "Poi"),
    tau.d = 0, tau.z = 0, # Var of random effects = 0 by default
    theta = NULL,     # shape parameter of NB distribution, Var = u + (u^2)/theta
    p.zero = NULL,     # P of excess-zeros for ctrl and trmt groups
    ydist,             # distribution of outcome model ("gaussian"/"binom")
    sd,                # if dist = "guassian", sd denotes the SD of the outcome variable
    py = NULL,         # P(Y=1) if outcome is binary
    coef.y = 0         # vector of coefficients in outcome model
)
{
  # set.seed(seed)

  if (missing(n.ind)) n.ind = 1
  if (missing(n.measure)) n.measure = n

  n <- n.ind * n.measure                                  # total number of observations
  ind <- sort(rep(1:n.ind, n.measure))                    # subject ID

  # ------------------------------------------------------------------------------------------------
  ## Simulate exposure variable
  # ------------------------------------------------------------------------------------------------
  if (missing(px) & missing(x.d)) {
    warning("Note: a bianry distribution and p=0.5 will be used by default")
    px = c(0.5, 0.5)                       # binary exposure with p=0.5 by default
    x.d <- c(rep(0, n/2), rep(1, n/2))     # generate a binary exposure with p=0.5 by default
  }
  if (!missing(px) & missing(x.d)){
    x.d <- list()
    for (i in 1:length(px)){
      x.d[[i]] <- rep((i-1), n*px[i])
    }
    x.d <- unlist(x.d)
  }

  x.d <- as.matrix(cbind(1, x.d) )      # combine all predictors in M model (both logit and log models)
  coef.d <- as.matrix(coef.d)

  if (!(ncol(x.d) == length(coef.d))){
    stop("Number of coefficients and predictors not match in log link")
  }

  # ------------------------------------------------------------------------------------------------
  ## Simulate Covariates
  # ------------------------------------------------------------------------------------------------
  if (missing(x.z) & is.null(p.zero)) x.z <- x.d
  if (missing(x.z) & !is.null(p.zero)) x.z <- 0
  if (missing(coef.z)) coef.z = NULL

  if (is.null(coef.z)){
    x.z = 0
  } else {
    x.z <- x.z
    coef.z <- as.matrix(coef.z)
    if (!(ncol(x.z) == length(coef.z))){
      stop("Number of coefficients and predictors not match in logit link")
    }
  }

  if (missing(tau.d)) tau.d <- 0   # Var of random effects = 0 by default at present
  if (missing(tau.z)) tau.z <- 0

  # ------------------------------------------------------------------------------------------------
  ## Simulate Mediator
  # ------------------------------------------------------------------------------------------------
  if (missing(p.zero)) p.zero = NULL
  if (missing(py)) py = NULL

  varx = x.d[,2]
  x.lvl = nlevels(as.factor(varx))

  if (!is.null(p.zero)){
    if (length(p.zero) == 1 & x.lvl > 1){
      p.zero[1:x.lvl] = p.zero[1]
    } else if(length(p.zero) != x.lvl & length(p.zero) >1 ) {
      stop("Number of 'p.zero' not match 'X levels' ")
    }
  }

  ## Mediator Simulation Function
  msim <- function(x.lvl){

    if (is.null(coef.z) & is.null(p.zero)){
      x.z = 0
    } else if (!is.null(coef.z) & is.null(p.zero)){
      x.z <- x.z[varx == unique(varx)[x.lvl], ]
    } else if (is.null(coef.z) & !is.null(p.zero)){
      x.z = 0
    }
    x.d <- x.d[varx == unique(varx)[x.lvl], ]

    n = n*px[x.lvl]
    ind <- sort(rep(1:n.ind, n))

    b <- rep(NA, n.ind)                  # random effect: b ~ N(0, tau.z^2)
    for (j in 1:n.ind) b[j] <- rnorm(1, 0, tau.z)

    ## calculate the excess-zero proportion given each X-level

    if (is.null(p.zero)){
      etaz <- b[ind] + x.z %*% coef.z
      m.normal <- rnorm(n, -etaz, 1.6)

      lvl <- (as.factor(x.d[,2]))[x.lvl]
      subdat <- x.d[x.d[,2] == lvl, ]
      p.zero = exp(subdat %*% coef.z)/(1+exp(subdat %*% coef.z))

      p.zero = round(mean(p.zero),2)
      quantiles <- quantile(m.normal, p.zero)
      m.z <-  as.numeric( factor(cut(m.normal, breaks = c(-Inf, quantiles, Inf))) ) - 1
    } else{
      p.zero = p.zero[x.lvl]
      m.z = rep(1, length(x.d[,1]))
    }

    ## simulate M given zero proportion
    b <- rep(NA, n.ind)                  # random effect: b ~ N(0, tau.c^2)
    for (j in 1:n.ind) b[j] <- rnorm(1, 0, tau.d)

    etad <- b[ind] + x.d %*% coef.d
    etad <- mean(etad)
    if (mdist %in% c("ZINB", "NB") )  m.nb <- rnbinom(n, mu = exp(etad), size = theta)
    if (mdist %in% c("ZIP", "Poi") )  m.nb <- rpois(n, lambda =  exp(etad))

    mnb <- ifelse(m.z == 0, 0, m.nb)

    return(mnb)
  }


  ls <- list()
  for(i in 1:x.lvl) ls[[i]] <- msim(x.lvl = i)
  # comb_df <- do.call(cbind, ls)
  comb_m <- c(ls[[1]], ls[[2]])
  ## Consider log-transformation of M
  dfm <- cbind(int = 1, x = x.d[1:n,2], mnb = comb_m, c = x.d[1:n,-(1:2)])
  colnames(dfm) <- c("int", "x", "mnb", colnames(c))

  ## -----------------------------------------------------------------------------------------------
  ## Simulate Outcome
  ## -----------------------------------------------------------------------------------------------
  xm <- c(dfm[,"x"]*dfm[,"mnb"])
  dfm <- data.frame(cbind(dfm, xm)) %>%
    mutate(im = ifelse(mnb == 0, 0, 1)) %>% mutate(xi = x*im) %>%
    as.matrix()
  expit <- function(x) exp(x)/(1 + exp(x))

  if (!(length(coef.y) == ncol(dfm))){
    stop("Coefficients in y model does not match number of predictors")
  }

  if (missing(ydist)) ydist = "gaussian"

  if (missing(sd) & ydist == "gaussian") sd = sd
  if (ydist == "gaussian"){
    y = rnorm(n, mean = dfm %*% coef.y, sd = sd)

  }else if (ydist == "binom"){

    if ( is.null(py) ){
      py0 <- expit(dfm[x==0,] %*% coef.y)
      py1 <- expit(dfm[x==1,] %*% coef.y)
    } else if ( !(is.null(py)) ){
      py0 <- py[1]
      py1 <- py[2]
    }
    p0 = mean(py0);     p1 = mean(py1)
    eta_0 <- dfm[x==0,] %*% coef.y
    y.normal <- rnorm(n/2, -eta_0, 1.6)
    quantiles <- quantile(y.normal, p0)
    y0 <- 2- as.numeric( factor(cut(y.normal, breaks = c(-Inf, quantiles, Inf))) )

    eta_1 <- dfm[x==1,] %*% coef.y
    y.normal <- rnorm(n/2, -eta_1, 1.6)
    quantiles <- quantile(y.normal, p1)
    y1 <- 2- as.numeric( factor(cut(y.normal, breaks = c(-Inf, quantiles, Inf))) )
    y <- c(y0, y1)
  }
  dat <- data.frame(cbind(dfm, y))

  return(dat)
}

# Example for simulating a dataset
## First define the sample size n
n = 1000
## Second define a binary exposure X, p = 0.50
x <- c(rep(0, n/2), rep(1, n/2))
## Assume we have a binary variable called covs with p=0.4
c = matrix(rbinom(n, 1, 0.4), ncol = 1, dimnames = list(NULL, "covs"))
x.d <- cbind(x, c)
mdist = "ZINB"  ## Define distribution of M, should be one of ZINB, NB, ZIP, Poi

# Example: simualte a data with ZI% of M ~ 20% ----
b0 = -2.0;  bx = 0.30;  bc = 0.2; bm = 0.06; bi = 0.25; bxi = 0; bxm = 0
g0 = -1.0;  gx = -1.5;  gc = 0.4
a0 =  3.0;  ax = 0.25;  ac = 0.35; theta = 5

if(mdist == "ZINB"){
  data1 <- dat.sim(
    n.ind = 1, n.measure = n,
    px = c(0.5, 0.5),
    x.d = x.d,
    coef.z = c(g0, gx, gc),
    coef.d = c(a0, ax, ac),
    mdist = "ZINB",
    theta = theta,
    ydist = "binom",
    sd = 1,
    coef.y = c(b0, bx, bm, bc, bxm, bi, bxi) )
} else if(mdist == "NB"){
  data1 <- dat.sim(
    n.ind = 1, n.measure = n,
    px = c(0.5, 0.5),
    x.d = x.d,
    p.zero = c(0,0),
    coef.d = c(a0, ax, ac),
    mdist = "NB",
    theta = 5,
    ydist = "binom",
    sd = 1,
    coef.y = c(b0, bx, bm, bc, bxm, bi, bxi) )
}


