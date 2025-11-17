#'
#' @name medbayes
#' @title Run Bayesian mediation model
#'
#' @param model.m A brms model for the mediator
#' @param dat.new New data for prediction
#' @return Posterior predictions (first column)
#' @export
medbayes_zim <- function(model.m = model.m, model.y = model.y,
                     treat = "treatment", mediator = "mediator", ind_mediator = NULL, outcome = "outcome",
                     control.value = 0, treat.value = 1,
                     logM = FALSE )
{
  library(brms)
  # Check whether package 'brms' is installed
  if (!requireNamespace("brms", quietly = TRUE))
    stop("Package 'brms' is required but not installed.")
  message("Mediation analysis function is running...")

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
    dat.y.temp[, mediator] = 0
    predict.y.cov.mu[i,,] = posterior_linpred(model.y, newdata = dat.y.temp)
  }

  # pred_y for treat and mediator(s) effect ----
  coef.mediator = paste("b_", mediator, sep = "")
  bm = as.matrix(ef_y)[,coef.mediator]

  if( grepl("zero", family(model.y)$family) ){
    coef.mediator.zi = paste("b_zi_", mediator, sep = "")
    b_zi_m = as.matrix(ef_y)[,coef.mediator.zi]
  }
  if( grepl("hurdle", family(model.y)$family) ){
    coef.mediator.zi = paste("b_hu_", mediator, sep = "")
    b_zi_m = as.matrix(ef_y)[,coef.mediator.zi]
  }

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
        as.matrix(predict.ms[j,,]*bm) + as.matrix(bxm*int_xm*value[i]*predict.ms[j,,]) )  +
        predict.y.cov.mu[i,,]
    }else{
      as.numeric(
      as.matrix(predict.ms[j,,]*bm) + as.matrix(bxm*int_xm*value[i]*predict.ms[j,,]) +
      as.matrix(b.ind_m*predict.ind_ms[r,,] + bxIm*int_xIm*value[i]*predict.ind_ms[r,,]) ) +
      predict.y.cov.mu[i,,]
    }
  }

  for(i in 1:length(value)){
    for(j in 1:2){
      for(r in 1:2){
        outcome.linpred.mu[i,j,r,,] =  calc_linpred.mu(i,j,r, bm, bxm, int_of_xm, int_of_xIm, b.ind_m )
      }
    }
  }

  # ----
  ### for i, j, m:
  ### 1,1,1,, == x0,m0,Im0
  ### 1,1,2,, == x0,m0,Im1
  ### 1,2,2,, == x0,m1,Im1
  ### 2,1,2,, == x1,m0,Im1 etc.

  # ----
  y_link = family(model.y)$link

  zi.outcome = grepl("zero", family(model.y)$family) | grepl("hurdle", family(model.y)$family)
  if(y_link == "identity" & !(zi.outcome))  outcome.pred = outcome.linpred.mu
  if(y_link == "logit"& !(zi.outcome))     outcome.pred  = exp(outcome.linpred.mu)/(1+exp(outcome.linpred.mu))

  if(!is.null(ind_mediator)) {

    if(y_link == "logit") {
      res.rd = cal.rd.effects.ind(outcome.pred)
      res.rr = cal.rr.effects.ind(outcome.pred)
      rst <- list(effects.rd = res.rd, effects.rr = res.rr, outcome.pred = outcome.pred)

    } else{
      res = cal.rd.effects.ind(outcome.pred)
      rst <- list(effects = res, outcome.pred = outcome.pred)
    }
  } else {

    if(y_link == "logit") {
      res.rd = cal.rd.effects(outcome.pred)
      res.rr = cal.rr.effects(outcome.pred)
      rst <- list(effects.rd = res.rd, effects.rr = res.rr, outcome.linpred=outcome.linpred.mu, outcome.pred = outcome.pred)
    } else if(zi.outcome){

      res.rr  = cal.rr.effects.y(outcome.pred.mu, outcome.pred.zi)

      res.mu.zi = cal.rd.effects.y(outcome.pred.mu, outcome.pred.zi)
      res.rd <- list(effects.mu.rd = res.mu.zi)
      rst <- list(effects.rd = res.rd, effects.rr = res.rr, outcome.linpred.mu=outcome.linpred.mu,
                  outcome.linpred.zi=outcome.linpred.zi,
                  outcome.pred.mu= outcome.pred.mu,  outcome.pred.zi =  outcome.pred.zi)

    } else{
      res = cal.rd.effects(outcome.pred)
      rst <- list(effects = res, outcome.linpred=outcome.linpred.mu, outcome.pred = outcome.pred)
    }
  }
  # rst

  param.terms <- list(
    call = match.call(),
    model.m=model.m,
    model.y=model.y,
    treat=treat,
    mediator = mediator,
    ind_mediator = ind_mediator,
    outcome = outcome,
    control.value = control.value,
    treat.value = treat.value
    # params = list(...),   # if you want to capture extra arguments
    # estimates = compute_effects(model_m, model_y)
  )
  class(param.terms) <- "params"

  return(list(Results = rst, params = param.terms))
}

cal.rd.effects.ind  <- function(outcome.pred)
{
  direct_control = outcome.pred[2,1,1,,] - outcome.pred[1,1,1,,]
  direct_treated = outcome.pred[2,2,1,,] - outcome.pred[1,2,1,,]
  indirect_control = outcome.pred[1,2,2,,] - outcome.pred[1,1,2,,]
  indirect_treated = outcome.pred[2,2,2,,] - outcome.pred[2,1,2,,]

  direct = (direct_control + direct_treated)/2
  indirect = (indirect_control + indirect_treated)/2

  indirect_Im = outcome.pred[2,1,2,,] - outcome.pred[2,1,1,,]
  indirect = indirect + indirect_Im

  total = direct + indirect
  pmed = indirect/total
  pmed.mean = mean(indirect)/mean(total)
  pmed.med = median(indirect)/median(total)

  res = rbind(
    c(mean(indirect_control), median(indirect_control), sd(indirect_control),
      quantile(indirect_control, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_control<0), mean(indirect_control>0))),

    c(mean(indirect_treated), median(indirect_treated), sd(indirect_treated),
      quantile(indirect_treated, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_treated<0), mean(indirect_treated>0))),

    # *******************************************************
    c(mean(indirect_Im), median(indirect_Im), sd(indirect_Im),
      quantile(indirect_Im, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_Im<0), mean(indirect_Im>0))),
    # *******************************************************

    c(mean(direct_control), median(direct_control), sd(direct_control),
      quantile(direct_control, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct_control<0), mean(direct_control>0))),

    c(mean(direct_treated), median(direct_treated), sd(direct_treated),
      quantile(direct_treated, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct_treated<0), mean(direct_treated>0))),

    c(mean(total), median(total), sd(total),
      quantile(total, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(total<0), mean(total>0))),

    c(mean(indirect), median(indirect), sd(indirect),
      quantile(indirect, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect<0), mean(indirect>0))),

    c(mean(direct), median(direct), sd(direct),
      quantile(direct, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct<0), mean(direct>0))),

    c(mean(pmed.mean), median(pmed.med), sd(pmed),
      quantile(pmed, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(pmed<0), mean(pmed>0)))
  ) # Bayes p-value: tail probability (see JMbayes), 2*min{pr(b<0), pr(b>0))}
  res[,1:5] = round(res[,1:5], digits=3)
  res[,6] = signif(res[,6], digits=4)
  res[9,1] = ifelse(res[9,1] < 0, 0, res[9,1] )
  rownames(res) = c("Indirect_control", "Indirect_treated", "Indirect_Indicator",
                    "Direct_control", "Direct_treated",
                    "Total Effect", "Indirect", "Direct", "Prop.Med")
  colnames(res) = c("Mean", "Median", "SD", "l-95% CI", "u-95% CI", "Bayes_p")
  # if(is.null(ind_mediator)){
  #   res = res[-c("Indirect_Indicator")]
  # }else{
    res = res
  # }
  res
}
cal.rr.effects.ind <- function(outcome.pred)
{
  direct_control = outcome.pred[2,1,1,,] / outcome.pred[1,1,1,,]
  direct_treated = outcome.pred[2,2,1,,] / outcome.pred[1,2,1,,]
  indirect_control = outcome.pred[1,2,2,,] / outcome.pred[1,1,2,,]
  indirect_treated = outcome.pred[2,2,2,,] / outcome.pred[2,1,2,,]

  direct = (direct_control + direct_treated)/2
  indirect = (indirect_control + indirect_treated)/2

  # **************************************************
  indirect_Im = outcome.pred[2,1,2,,] / outcome.pred[2,1,1,,]
  indirect = indirect * indirect_Im
  total = direct*indirect

  # **************************************************
  pmed = direct*(indirect-1)/(total-1)
  pmed.mean = mean(direct)*(mean(indirect)-1)/(mean(direct)*mean(indirect)-1)
  pmed.med = median(direct)*(median(indirect)-1)/(median(direct)*median(indirect)-1)

  res = rbind(
    c(mean(indirect_control), median(indirect_control), sd(indirect_control),
      quantile(indirect_control, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_control<1), mean(indirect_control>1))),

    c(mean(indirect_treated), median(indirect_treated), sd(indirect_treated),
      quantile(indirect_treated, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_treated<1), mean(indirect_treated>1))),

    # *******************************************************
    c(mean(indirect_Im), median(indirect_Im), sd(indirect_Im),
      quantile(indirect_Im, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_Im<1), mean(indirect_Im>1))),
    # *******************************************************

    c(mean(direct_control), median(direct_control), sd(direct_control),
      quantile(direct_control, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct_control<1), mean(direct_control>1))),

    c(mean(direct_treated), median(direct_treated), sd(direct_treated),
      quantile(direct_treated, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct_treated<1), mean(direct_treated>1))),

    c(mean(total), median(total), sd(total),
      quantile(total, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(total<1), mean(total>1))),

    c(mean(indirect), median(indirect), sd(indirect),
      quantile(indirect, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect<1), mean(indirect>1))),

    c(mean(direct), median(direct), sd(direct),
      quantile(direct, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct<1), mean(direct>1))),

    c(mean(pmed.mean), median(pmed.med), sd(pmed),
      quantile(pmed, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(pmed<0), mean(pmed>0)))
  ) # Bayes p-value: tail probability (see JMbayes), 2*min{pr(b<0), pr(b>0))}
  res[,1:5] = round(res[,1:5], digits=3)
  res[,6] = signif(res[,6], digits=4)
  res[9,1] = ifelse(res[9,1] < 0, 0, res[9,1] )

  rownames(res) = c("Indirect_control", "Indirect_treated", "Indirect_Indicator",
                    "Direct_control", "Direct_treated",
                    "Total Effect", "Indirect", "Direct", "Prop.Med")
  colnames(res) = c("Mean", "Median", "SD", "l-95% CI", "u-95% CI", "Bayes_p")
  # if(is.null(ind_mediator)){
  #   res = res[-c("Indirect_Indicator")]
  # }else{
    res = res
  # }
  res
}



