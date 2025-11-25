#'
#' @name medbayes
#' @title Run Bayesian mediation model
#'
#' @param model.m A brms model for the mediator
#' @param model.y A brms model for the outcome
#' @param dat.new New data for prediction
#' @return Posterior predictions (first column)
#' @export
medbayes_zim <- function(model.m = model.m, model.y = model.y,
                     treat = "treatment", mediator = "mediator",
                     ind_mediator = NULL, outcome = "outcome",
                     control.value = 0, treat.value = 1 )
{
  library(brms)
  # Check whether package 'brms' is installed
  if (!requireNamespace("brms", quietly = TRUE))
    stop("Package 'brms' is required but not installed.")
  message("Mediation analysis function is running...")

  if(!(grepl("hurdle", family(model.m)$family) | grepl("zero", family(model.m)$family))) ind_mediator = NULL

  value = c(control.value, treat.value)
  dat.new = model.m$data

  dpar.m = c("mu", "hu")
  predict.ms <- array(NA, dim = c(2, ndraws(model.m), nrow(dat.new)))
  for(i in 1:length(value)){
    dat.new[, treat] = value[i]
    predict.ms[i,,] = 0
    predict.ms[i,,] = posterior_epred(model.m, newdata = dat.new, dpar = NULL)
    }

  if(!is.null(ind_mediator)){
    predict.ind_ms <- array(NA, dim = c(2, ndraws(model.m), nrow(dat.new)))
    for(i in 1:length(value)){
      dat.new[, treat] = value[i]
      predict.ind_ms[i,,] = 0
      predict.ind_ms[i,,] = 1- posterior_epred(model.m, newdata = dat.new, dpar = dpar.m[2])
    }
  }else{
     predict.ind_ms <- array(0, dim = c(2, ndraws(model.m), 1))
  }

  dat.y = model.y$data
  ef_y = as_draws_df(model.y)

  depar.outcome = NULL
  if(grepl("zero", family(model.y)$family)) depar.outcome = c("mu", "zi")
  if(grepl("hurdle", family(model.y)$family)) depar.outcome = c("mu", "hu")

  # pred_y without mediator(s) effect ----
  predict.y.cov.mu <- array(NA, dim = c(2, ndraws(model.y), nrow(dat.y)) )

  for(i in 1:length(value)){
    dat.y.temp <- dat.y
    dat.y.temp[, treat] <- value[i]
    # dat.y.temp[, mediator] = 0

    if(!is.numeric(dat.y.temp[, mediator])){
      dat.y.temp[, mediator] = levels(factor(dat.y.temp[, mediator]))[1]
    }else if(is.numeric(dat.y.temp[, mediator])){
      dat.y.temp[, mediator] = 0
    }

    predict.y.cov.mu[i,,] = posterior_linpred(model.y, newdata = dat.y.temp)
  }

  # pred_y for treat and mediator(s) effect ----

  if(is.numeric(dat.y.temp[, mediator])){
    coef.mediator = paste("b_", mediator, sep = "")
  }else if(!is.numeric(dat.y.temp[, mediator])){
    coef.med.starts = paste("b_", mediator, sep = "")
    coef.mediator <- colnames(ef_y)[startsWith(colnames(ef_y), coef.med.starts)]
  }
  # coef.mediator = paste("b_", mediator, sep = "")
  bm = as.matrix(ef_y)[,coef.mediator]

  if(!is.null(ind_mediator)){
    coef.ind_mediator = paste("b_", ind_mediator, sep = "")
    b.ind_m = as.matrix(ef_y)[,coef.ind_mediator]

    if( grepl("zero", family(model.y)$family) ){
      coef.ind_mediator.zi = paste("b_zi_", ind_mediator, sep = "")
      b_zi_indm = as.matrix(ef_y)[,coef.ind_mediator.zi]
    }
    if( grepl("hurdle", family(model.y)$family) ){
      coef.ind_mediator.zi = paste("b_hu_", ind_mediator, sep = "")
      b_zi_indm = as.matrix(ef_y)[,coef.ind_mediator.zi]
    }
  }else{
    b.ind_m = 0
    b_zi_indm = 0
  }

  ## check interaction between treatment and mediator/Indicator of mediator
  INT.xm  <- c(paste("b_", treat, ":", mediator, sep="") %in% colnames(as.matrix(model.y)),
              paste("b_", mediator, ":", treat, sep="") %in% colnames(as.matrix(model.y))  )
  INT.xIm <- c(paste("b_", treat, ":", ind_mediator, sep="") %in% colnames(as.matrix(model.y)) ,
               paste("b_", ind_mediator, ":", treat, sep="") %in% colnames(as.matrix(model.y))  )

  INT.xm_hu <- c(paste("b_hu_", treat, ":", mediator, sep="") %in% colnames(as.matrix(model.y)) ,
                 paste("b_hu_", mediator, ":", treat, sep="") %in% colnames(as.matrix(model.y))  )
  INT.xm_zi <- c(paste("b_zi_", treat, ":", mediator, sep="") %in% colnames(as.matrix(model.y)) ,
                 paste("b_zi_", mediator, ":", treat, sep="") %in% colnames(as.matrix(model.y))  )

  coef.intxm = c(paste("b_", treat, ":", mediator, sep=""), paste("b_", mediator, ":", treat, sep="") )
  coef.intxm = coef.intxm[INT.xm]

  coef.intxIm = c(paste("b_", treat, ":", ind_mediator, sep=""), paste("b_", ind_mediator, ":", treat, sep=""))
  coef.intxIm = coef.intxIm[INT.xIm]

  int_of_xm <- any(coef.intxm %in% colnames(as.matrix(model.y)))
  int_of_xIm <- any(coef.intxIm %in% colnames(as.matrix(model.y)))

  if(int_of_xm) {
    bxm = as.matrix(ef_y)[,coef.intxm]
  } else(
    bxm = 0
  )

  if(int_of_xIm ) {
    bxIm = as.matrix(ef_y)[,coef.intxIm]
  } else {
    bxIm = 0
  }

  outcome.linpred.mu = array(NA, dim = c(2, 2, 2, ndraws(model.y), nrow(dat.y)))

  if(is.character(value)){
    comp.val = order(levels(factor(value)))-1
  }else{
    comp.val = value
  }

  calc_linpred.mu <- function(i, j, r, bm, bxm, int_of_xm, int_of_xIm, b.ind_m ){
    if(int_of_xm){
      int_xm = 1
    } else {
      int_xm=0
      bxm=0
    }
    if(int_of_xIm){
      int_xIm = 1
    } else {
      int_xIm=0
      bxIm=0}

    if(is.null(ind_mediator)){
      as.numeric(
        as.matrix(predict.ms[j,,]*bm) + as.matrix(bxm*int_xm*comp.val[i]*predict.ms[j,,]) )  +
        predict.y.cov.mu[i,,]
    }else{
      as.numeric(
      as.matrix(predict.ms[j,,]*bm) + as.matrix(bxm*int_xm*comp.val[i]*predict.ms[j,,]) +
      as.matrix(b.ind_m*predict.ind_ms[r,,] + bxIm*int_xIm*comp.val[i]*predict.ind_ms[r,,]) ) +
      predict.y.cov.mu[i,,]
    }
  }

  for(i in 1:length(comp.val)){
    for(j in 1:2){
      for(r in 1:2){
        outcome.linpred.mu[i,j,r,,] =  calc_linpred.mu(i,j,r, bm, bxm, int_of_xm, int_of_xIm, b.ind_m )
      }
    }
  }

  # ----
  y_link = family(model.y)$link

  if(y_link == "identity")  outcome.pred = outcome.linpred.mu
  if(y_link == "logit")     outcome.pred  = exp(outcome.linpred.mu)/(1+exp(outcome.linpred.mu))

  if(!is.null(ind_mediator)) {

    if(y_link == "logit") {
      res.rd = cal.rd.effects.ind(outcome.pred, ind_mediator)
      res.rr = cal.rr.effects.ind(outcome.pred, ind_mediator)
      rst <- list(effects.rd = res.rd, effects.rr = res.rr, outcome.linpred=outcome.linpred.mu, outcome.pred = outcome.pred)

    }else{
      res = cal.rd.effects.ind(outcome.pred, ind_mediator)
      rst <- list(effects = res, outcome.linpred=outcome.linpred.mu, outcome.pred = outcome.pred)
    }
  } else {

    if(y_link == "logit") {
      res.rd = cal.rd.effects.ind(outcome.pred, ind_mediator)
      res.rr = cal.rr.effects.ind(outcome.pred, ind_mediator)
      rst <- list(effects.rd = res.rd, effects.rr = res.rr, outcome.linpred=outcome.linpred.mu, outcome.pred = outcome.pred)
    } else{
      res.rd = cal.rd.effects.ind(outcome.pred, ind_mediator)
      res.rr = cal.rr.effects.ind(outcome.pred, ind_mediator)
      rst <- list(effects.rd = res.rd, effects.rr = res.rr, outcome.linpred=outcome.linpred.mu, outcome.pred = outcome.pred)
    }
  }

  param.terms <- list(
    call = match.call(),
    model.m=model.m,
    model.y=model.y,
    treat=treat,
    mediator = mediator,
    if(!is.null(ind_mediator)){
      ind_mediator = ind_mediator
    } else{
      ind_mediator = NULL
    },
    outcome = outcome,
    control.value = control.value,
    treat.value = treat.value
  )
  class(param.terms) <- "params"
  return(list(Results = rst, params = param.terms))
}

cal.rd.effects.ind  <- function(outcome.pred, ind_mediator)
{
  direct = outcome.pred[2,1,1,,] - outcome.pred[1,1,1,,]
  indirect_treated = outcome.pred[2,2,2,,] - outcome.pred[2,1,2,,]

  # **************************************************
  indirect_Im = outcome.pred[2,1,2,,] - outcome.pred[2,1,1,,]

  indirect_total = indirect_treated + indirect_Im
  total = direct+indirect_total

  # **************************************************
  pmed = indirect_total/total
  # pmed.mean = mean(direct)*(mean(indirect)-1)/(mean(direct)*mean(indirect)-1)
  pmed.med = median(indirect_total)/median(total)

  res = rbind(

    c( median(indirect_treated), sd(indirect_treated),
       quantile(indirect_treated, probs=c(0.025,0.975), na.rm = T),
       2*min(mean(indirect_treated<0), mean(indirect_treated>0))),

    # *******************************************************
    c(median(indirect_Im), sd(indirect_Im),
      quantile(indirect_Im, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_Im<0), mean(indirect_Im>0))),
    # *******************************************************

    c(median(indirect_total), sd(indirect_total),
          quantile(indirect_total, probs=c(0.025,0.975), na.rm = T),
          2*min(mean(indirect_total<0), mean(indirect_total>0))),

    c(median(direct), sd(direct),
      quantile(direct, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct<0), mean(direct>0))),

    c(median(total), sd(total),
      quantile(total, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(total<0), mean(total>0))),

    c(median(pmed.med), sd(pmed),
      quantile(pmed, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(pmed<0), mean(pmed>0)))
  ) # Bayes p-value: tail probability (see JMbayes), 2*min{pr(b<0), pr(b>0))}
  res[,1:4] = round(res[,1:4], digits=3)
  res[,5] = round(res[,5], digits=4)
  res[6,1] = ifelse(res[6,1] < 0, 0, res[6,1] )

  rownames(res) = c("Indirect_treated", "Indirect_Indicator",
                    "Indirect_total", "Direct", "Total Effect", "Prop.Med")
  colnames(res) = c("Estimate", "SD", "l-95% CI", "u-95% CI", "Bayes_p")
  if(is.null(ind_mediator)){
    res = res[rownames(res) != "Indirect_Indicator", ]
  }else{
    res = res
  }
  res
}

cal.rr.effects.ind <- function(outcome.pred, ind_mediator)
{
  direct = outcome.pred[2,1,1,,] / outcome.pred[1,1,1,,]
  indirect_treated = outcome.pred[2,2,2,,] / outcome.pred[2,1,2,,]

  # **************************************************
  indirect_Im = outcome.pred[2,1,2,,] / outcome.pred[2,1,1,,]

  indirect_total = indirect_treated * indirect_Im
  total = direct*indirect_total

  # **************************************************
  pmed = direct*(indirect_total-1)/(total-1)
  # pmed.mean = mean(direct)*(mean(indirect)-1)/(mean(direct)*mean(indirect)-1)
  pmed.med = median(direct)*(median(indirect_total)-1)/(median(direct)*median(indirect_total)-1)

  res = rbind(

    c( median(indirect_treated), sd(indirect_treated),
      quantile(indirect_treated, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_treated<1), mean(indirect_treated>1))),

    # *******************************************************
    c(median(indirect_Im), sd(indirect_Im),
      quantile(indirect_Im, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_Im<1), mean(indirect_Im>1))),
    # *******************************************************

    c(median(indirect_total), sd(indirect_total),
      quantile(indirect_total, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_total<1), mean(indirect_total>1))),

    c(median(direct), sd(direct),
      quantile(direct, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct<1), mean(direct>1))),

    c(median(total), sd(total),
      quantile(total, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(total<1), mean(total>1))),

    c(median(pmed.med), sd(pmed),
      quantile(pmed, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(pmed<0), mean(pmed>0)))
  ) # Bayes p-value: tail probability (see JMbayes), 2*min{pr(b<0), pr(b>0))}
  res[,1:4] = round(res[,1:4], digits=3)
  res[,5] = round(res[,5], digits=4)
  res[6,1] = ifelse(res[6,1] < 0, 0, res[6,1] )

  rownames(res) = c("Indirect_treated", "Indirect_Indicator",
                    "Indirect_total", "Direct", "Total Effect", "Prop.Med")
  colnames(res) = c("Estimate", "SD", "l-95% CI", "u-95% CI", "Bayes_p")
  if(is.null(ind_mediator)){
    res = res[rownames(res) != "Indirect_Indicator", ]
  }else{
    res = res
  }
  res
}





