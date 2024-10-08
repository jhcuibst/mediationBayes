
library(brms)
# **************************************************
# This function only for zero-inflated mediator iwth hurdle model
medbayes <- function(model.m = hnb.m, model.y = hnb.y,
                     treat = "treatment", mediator = "mediator", ind_mediator = NULL, outcome = "outcome",
                     control.value = 0, treat.value = 1,
                     logM = FALSE )
{
  # Check whether package 'brms' is installed

  if (!require("brms", character.only = TRUE)) {
    # Ask the user if they want to install the package
    answer <- readline(prompt = "Package 'brms' is not installed. Do you want to install it? (yes/no): ")

    if (tolower(answer) == "yes" || answer == "y" || answer == "") {
      # Install the package
      install.packages("brms")
      # Load the package after installation
      library(brms, character.only = TRUE)
      message("Package 'brms' has been installed and loaded.")
    } else {
      stop("Package 'brms' is required but not installed. Exiting.")
    }
  } else {
    message("Package 'brms' is already installed and loaded.")
  }

  message("Your main function is running...")


  value = c(control.value, treat.value)

  # *********************************************************************
  # ----------------
  dat.new = model.m$data

  if (family(model.m)$family %in% c("hurdle_negbinomial", "hurdle_poisson")){
    dpar.m = c("mu", "hu")
    predict.ms <- array(NA, dim = c(2, ndraws(model.m), nrow(dat.new)))
    for(i in 1:length(value)){
      dat.new[, treat] = value[i]
      predict.ms[i,,] = 0
      # predict.ms[i,,] = posterior_epred(model.m, newdata = dat.new, dpar = dpar.m[1])
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
      predict.ind_ms <- array(0, dim = c(2, ndraws(model.m), nrow(dat.new)))
    }
  }

  # *********************************************************************
  # TBD ----------------
  if (family(model.m)$family %in% c("zero_inflated_negbinomial", "zero_inflated_poisson")){
    dpar.m = c("mu", "zi")
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
        predict.ind_ms[i,,] = 1 - posterior_epred(model.m, newdata = dat.new, dpar = dpar.m[2])
      }
    }else{
      predict.ind_ms <- array(0, dim = c(2, ndraws(model.m), nrow(dat.new)))
    }
  }

  # *********************************************************************
  # If model.m not belong to ZI/HU models ----------------
  hu.models = c("hurdle_negbinomial", "hurdle_poisson","zero_inflated_negbinomial", "zero_inflated_poisson")
  if (!(family(model.m)$family %in% hu.models)) {
    predict.ms <- array(NA, dim = c(2, ndraws(model.m), nrow(dat.new)))
    predict.ind_ms <- array(0, dim = c(2, ndraws(model.m), nrow(dat.new)))
    for(i in 1:length(value)){
      dat.new[, treat] = value[i]
      predict.ms[i,,] = 0
      predict.ms[i,,] = posterior_epred(model.m, newdata = dat.new, resp = mediator)
    }
  }

  # *********************************************************************
  # ----------------
  # if(logM){
  #   predict.ind_ms = ifelse(predict.ind_ms == 0, 0.001, predict.ind_ms)
  #   predict.ind_ms = log(predict.ind_ms)
  # }
  # *********************************************************************-----

  dat.y = model.y$data
  ef_y = as_draws_df(model.y)

  depar.outcome = NULL
  if(grepl("zero", family(model.y)$family)) depar.outcome = c("mu", "zi")
  if(grepl("hurdle", family(model.y)$family)) depar.outcome = c("mu", "hu")

  zi.outcome = grepl("zero", family(model.y)$family) | grepl("hurdle", family(model.y)$family)
  # pred_y without mediator(s) effect ----
  predict.y.cov <- array(NA, dim = c(2, ndraws(model.y), nrow(dat.y)) )
  predict.y.cov.mu <- array(NA, dim = c(2, ndraws(model.y), nrow(dat.y)) )
  predict.y.cov.zi <- array(NA, dim = c(2, ndraws(model.y), nrow(dat.y)) )

  for(i in 1:length(value)){
    dat.y[, treat] = value[i]
    dat.y[, mediator] = 0
    predict.y.cov.mu[i,,] = 0
    predict.y.cov.zi[i,,] = 0
    predict.y.cov[i,,] = 0

    if( zi.outcome  ){
      predict.y.cov.mu[i,,] = posterior_linpred(model.y, newdata = dat.y, dpar = depar.outcome[1])
      predict.y.cov.zi[i,,] = posterior_linpred(model.y, newdata = dat.y, dpar = depar.outcome[2])
      predict.y.cov[i,,] = posterior_linpred(model.y, newdata = dat.y)
    }else{
      predict.y.cov[i,,] = posterior_linpred(model.y, newdata = dat.y)
    }
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
    b_zi_indm=0
  }

  # *********************************************************************
  ## check interaction between treatment and mediator(s)
  #******** wait to check interaction in the ZI-outcome model within ZI part *********
  #*
  INT.xm <- c(paste("b_", treat, ":", mediator, sep="") %in% colnames(as.matrix(model.y)),
              paste("b_", mediator, ":", treat, sep="") %in% colnames(as.matrix(model.y))
  )
  INT.xIm <- c(paste("b_", treat, ":", ind_mediator, sep="") %in% colnames(as.matrix(model.y)) ,
               paste("b_", ind_mediator, ":", treat, sep="") %in% colnames(as.matrix(model.y))
  )

  # *********************************************************************
  # if(any(INT.xm, INT.xIm)){
  coef.intxm = c(paste("b_", treat, ":", mediator, sep=""), paste("b_", mediator, ":", treat, sep="")  )
  coef.intxm = coef.intxm[INT.xm]

  coef.intxIm = c(paste("b_", treat, ":", ind_mediator, sep=""), paste("b_", ind_mediator, ":", treat, sep=""))
  coef.intxIm = coef.intxIm[INT.xIm]

  int_of_xm <- any(coef.intxm %in% colnames(as.matrix(model.y)))
  int_of_xIm <- any(coef.intxIm %in% colnames(as.matrix(model.y)))

  if( int_of_xm) {
    bxm = as.matrix(ef_y)[,coef.intxm]
  } else{
    bxm = 0
  }
  if( int_of_xIm ) {
    bxIm = as.matrix(ef_y)[,coef.intxIm]
  } else {
    bxIm = 0
  }

  #
  outcome.linpred = array(NA, dim = c(2, 2, 2, ndraws(model.y), nrow(dat.y)))
  outcome.linpred.mu = array(NA, dim = c(2, 2, 2, ndraws(model.y), nrow(dat.y)))
  outcome.linpred.zi = array(NA, dim = c(2, 2, 2, ndraws(model.y), nrow(dat.y)))

  if( int_of_xIm & int_of_xm ) {

    for(i in 1:length(value)){
      for(j in 1:2){
        for(r in 1:2){

          if(zi.outcome ){
            outcome.linpred.mu[i,j,r,,] =
              suppressWarnings(as.numeric(as.matrix(predict.ms[j,,]*bm) + as.matrix(bxm*value[i]*predict.ms[j,,]) +
                                            as.matrix(b.ind_m *predict.ind_ms[r,,] + bxIm*value[i]*predict.ind_ms[r,,]) ) +
                                 predict.y.cov.mu[i,,])
            outcome.linpred.zi[i,j,r,,] =
              suppressWarnings(as.numeric(as.matrix(predict.ms[j,,]*b_zi_m) + as.matrix(bxm*value[i]*predict.ms[j,,]) +
                                            as.matrix(b_zi_indm *predict.ind_ms[r,,] + bxIm*value[i]*predict.ind_ms[r,,]) ) +
                                 predict.y.cov.zi[i,,])
            outcome.linpred[i,j,r,,] = outcome.linpred.mu[i,j,r,,]*(1-outcome.linpred.zi[i,j,r,,])
          }else{
            outcome.linpred[i,j,r,,] =
              suppressWarnings(as.numeric(as.matrix(predict.ms[j,,]*bm) + as.matrix(bxm*value[i]*predict.ms[j,,]) +
                                            as.matrix(b.ind_m *predict.ind_ms[r,,] + bxIm*value[i]*predict.ind_ms[r,,]) ) +
                                 predict.y.cov[i,,])
          }
        }
      }
    }
  }
  #
  if(!(int_of_xIm) & int_of_xm) {

    for(i in 1:length(value)){
      for(j in 1:2){
        for(r in 1:2){

          if(zi.outcome){
            outcome.linpred.mu[i,j,r,,] =
              suppressWarnings(as.numeric(as.matrix(predict.ms[j,,]*bm) + as.matrix(bxm*value[i]*predict.ms[j,,]) +
                                            as.matrix(b.ind_m *predict.ind_ms[r,,] ) ) +  predict.y.cov.mu[i,,])
            outcome.linpred.zi[i,j,r,,] =
              suppressWarnings(as.numeric(as.matrix(predict.ms[j,,]*b_zi_m) + as.matrix(bxm*value[i]*predict.ms[j,,]) +
                                            as.matrix(b_zi_indm *predict.ind_ms[r,,])) + predict.y.cov.zi[i,,])
            outcome.linpred[i,j,r,,] = outcome.linpred.mu[i,j,r,,]*(1-outcome.linpred.zi[i,j,r,,])
          }else{
            outcome.linpred[i,j,r,,] =
              suppressWarnings(as.numeric(as.matrix(predict.ms[j,,]*bm) + as.matrix(bxm*value[i]*predict.ms[j,,]) +
                                            as.matrix(b.ind_m *predict.ind_ms[r,,] + bxIm*value[i]*predict.ind_ms[r,,]) ) +
                                 predict.y.cov[i,,])
          }
        }
      }
    }
  }
  #
  if( int_of_xIm & !(int_of_xm)) {

    for(i in 1:length(value)){
      for(j in 1:2){
        for(r in 1:2){

          if(zi.outcome){
            outcome.linpred.mu[i,j,r,,] =
              suppressWarnings(
                as.numeric(as.matrix(predict.ms[j,,]*bm) + as.matrix(b.ind_m *predict.ind_ms[r,,] + bxIm*value[i]*predict.ind_ms[r,,]) ) +
                  predict.y.cov.mu[i,,] )
            outcome.linpred.zi[i,j,r,,] =
              suppressWarnings(as.numeric(as.matrix(predict.ms[j,,]*b_zi_m) + as.matrix(b_zi_indm *predict.ind_ms[r,,]
                                                                                        + bxIm*value[i]*predict.ind_ms[r,,]) ) + predict.y.cov.zi[i,,])
            outcome.linpred[i,j,r,,] = outcome.linpred.mu[i,j,r,,]*(1-outcome.linpred.zi[i,j,r,,])
          }else{
            outcome.linpred[i,j,r,,] =
              suppressWarnings(as.numeric(as.matrix(predict.ms[j,,]*bm) + as.matrix(bxm*value[i]*predict.ms[j,,]) +
                                            as.matrix(b.ind_m *predict.ind_ms[r,,] + bxIm*value[i]*predict.ind_ms[r,,]) ) +
                                 predict.y.cov[i,,])
          }
        }
      }
    }
  }
  #
  if(!(int_of_xIm) & !(int_of_xm) ) {

    for(i in 1:length(value)){
      for(j in 1:2){
        outcome.linpred[i,j,,,] = 0
        for(r in 1:2){

          if(zi.outcome){
            outcome.linpred.mu[i,j,r,,] = suppressWarnings(
              as.numeric(as.matrix(predict.ms[j,,]*bm) + as.matrix(b.ind_m*predict.ind_ms[r,,]) ) + predict.y.cov.mu[i,,])
            outcome.linpred.zi[i,j,r,,] = suppressWarnings(
              as.numeric(as.matrix(predict.ms[j,,]*b_zi_m) + as.matrix(b_zi_indm *predict.ind_ms[r,,])) + predict.y.cov.zi[i,,])
            outcome.linpred[i,j,r,,] = outcome.linpred.mu[i,j,r,,]*(1-outcome.linpred.zi[i,j,r,,])
          }else{
            outcome.linpred[i,j,r,,] =
              suppressWarnings(as.numeric(as.matrix(predict.ms[j,,]*bm) + as.matrix(bxm*value[i]*predict.ms[j,,]) +
                                            as.matrix(b.ind_m *predict.ind_ms[r,,] + bxIm*value[i]*predict.ind_ms[r,,]) ) +
                                 predict.y.cov[i,,])
          }
        }
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

  if(y_link == "identity")  outcome.pred = outcome.linpred
  if(y_link == "logit")     outcome.pred  = exp(outcome.linpred)/(1+exp(outcome.linpred))

  zi.outcome = grepl("zero", family(model.y)$family) | grepl("hurdle", family(model.y)$family)
  if(zi.outcome & y_link == "log"){
    outcome.pred.mu  = exp(outcome.linpred.mu)
    outcome.pred.zi  = exp(outcome.linpred.zi)/(1+exp(outcome.linpred.zi))
    outcome.pred.overall  = outcome.pred.mu*(1-outcome.pred.zi)
  }else if(!(zi.outcome) & y_link == "log"){
    outcome.pred  = exp(outcome.linpred)
  }

  # res.all = cal.effects (outcome.pred)
  if(!is.null(ind_mediator)) {

    if(y_link == "logit") {
      res.rd = cal.effects.ind(outcome.pred)
      res.rr = cal.rr.effects.ind(outcome.pred)
      rst <- list(effects.rd = res.rd, effects.rr = res.rr, outcome.linpred=outcome.linpred, outcome.pred = outcome.pred)
    } else if(zi.outcome){
      res.mu.rr = cal.rr.effects.ind(outcome.pred.mu)
      res.zi.rr = cal.rr.effects.ind(outcome.pred.zi)
      res.overall.rr = cal.rr.effects.ind(outcome.pred.overall)
      res.rr <- list(effects.mu.rr = res.mu.rr, effects.zi.rr = res.zi.rr, effects.overall.rr = res.overall.rr)

      res.mu.rd = cal.effects.ind(outcome.pred.mu)
      res.zi.rd = cal.effects.ind(outcome.pred.zi)
      res.overall.rd = cal.effects.ind(outcome.pred.overall)
      res.rd <- list(effects.mu.rd = res.mu.rd, effects.zi.rd = res.zi.rd, effects.overall.rd = res.overall.rd)
      rst <- list(effects.rd = res.rd, effects.rr = res.rr, outcome.linpred.mu=outcome.linpred.mu, outcome.linpred.zi=outcome.linpred.zi)

    } else{
      res = cal.effects.ind(outcome.pred)
      rst <- list(effects = res, outcome.linpred=outcome.linpred, outcome.pred = outcome.pred)
    }

  } else {

    if(y_link == "logit") {
      res.rd = cal.diff.effects(outcome.pred)
      res.rr = cal.rr.effects(outcome.pred)
      rst <- list(effects.rd = res.rd, effects.rr = res.rr, outcome.linpred=outcome.linpred, outcome.pred = outcome.pred)
    } else if(zi.outcome){
      res.mu = cal.diff.effects(outcome.pred.mu)
      res.zi = cal.diff.effects(outcome.pred.zi)
      res.overall = cal.diff.effects(outcome.pred.overall)
      rst <- list(effects.mu = res.mu, effects.zi = res.zi, effects.overall = res.overall,
                  outcome.linpred.mu=outcome.linpred.mu, outcome.linpred.zi=outcome.linpred.zi)
    } else{
      res = cal.diff.effects(outcome.pred)
      rst <- list(effects = res, outcome.linpred=outcome.linpred, outcome.pred = outcome.pred)
    }
  }
  rst
}

cal.diff.effects <- function(outcome.pred)
{
  direct_control = outcome.pred[2,1,1,,] - outcome.pred[1,1,1,,]
  # treated: Y(1,M(1)) - Y(0,M(1))
  direct_treated = outcome.pred[2,2,1,,] - outcome.pred[1,2,1,,]
  # mediation effect: Y(t,M(1)) - Y(t,M(0))
  # control: Y(0,M(1)) - Y(0,M(0))
  indirect_control = outcome.pred[1,2,2,,] - outcome.pred[1,1,2,,]
  # treated: Y(1,M(1)) - Y(1,M(0))
  indirect_treated = outcome.pred[2,2,2,,] - outcome.pred[2,1,2,,]

  direct = (direct_control + direct_treated)/2
  indirect = (indirect_control + indirect_treated)/2

  total = direct + indirect
  # **************************************************
  pmed = indirect/total

  res = rbind(
    c(mean(indirect_control), median(indirect_control), sd(indirect_control),
      quantile(indirect_control, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_control<0), mean(indirect_control>0))),

    c(mean(indirect_treated), median(indirect_treated), sd(indirect_treated),
      quantile(indirect_treated, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_treated<0), mean(indirect_treated>0))),

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

    c(mean(pmed), median(pmed), sd(pmed),
      quantile(pmed, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(pmed<0), mean(pmed>0)))
  ) # Bayes p-value: tail probability (see JMbayes), 2*min{pr(b<0), pr(b>0))}
  res[,1:5] = round(res[,1:5], digits=3)
  res[,6] = signif(res[,6], digits=2)
  res[8,1] = ifelse(res[8,1] < 0, 0, res[8,1] )
  rownames(res) = c("Indirect_control", "Indirect_treated", "Direct_control", "Direct_treated",
                    "Total Effect", "Indirect", "Direct", "Prop.Med")
  colnames(res) = c("Mean", "Median", "Sd", "l-95% CI", "u-95% CI", "Bayes_p")

  res
}
cal.effects.ind  <- function(outcome.pred)
{
  direct_control = outcome.pred[2,1,1,,] - outcome.pred[1,1,1,,]
  # treated: Y(1,M(1)) - Y(0,M(1))
  direct_treated = outcome.pred[2,2,1,,] - outcome.pred[1,2,1,,]
  # mediation effect: Y(t,M(1)) - Y(t,M(0))
  # control: Y(0,M(1)) - Y(0,M(0))
  indirect_control = outcome.pred[1,2,2,,] - outcome.pred[1,1,2,,]
  # treated: Y(1,M(1)) - Y(1,M(0))
  indirect_treated = outcome.pred[2,2,2,,] - outcome.pred[2,1,2,,]

  direct = (direct_control + direct_treated)/2
  indirect = (indirect_control + indirect_treated)/2

  # **************************************************
  # treated: Y(1,M(0), I(m1)) - Y(1,M(0), I(m0))
  indirect_Im = outcome.pred[2,1,2,,] - outcome.pred[2,1,1,,]
  indirect = indirect + indirect_Im

  total = direct + indirect
  # **************************************************
  # t1 = outcome.pred[2,2,2,,] - outcome.pred[1,1,1,,]  # include Im
  # t2 = outcome.pred[2,2, ,,] - outcome.pred[1,1, ,,]  # does not include Im
  # **************************************************
  pmed = indirect/total

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

    c(mean(pmed), median(pmed), sd(pmed),
      quantile(pmed, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(pmed<0), mean(pmed>0)))
  ) # Bayes p-value: tail probability (see JMbayes), 2*min{pr(b<0), pr(b>0))}
  res[,1:5] = round(res[,1:5], digits=3)
  res[,6] = signif(res[,6], digits=2)
  res[9,1] = ifelse(res[9,1] < 0, 0, res[9,1] )
  rownames(res) = c("Indirect_control", "Indirect_treated", "Indirect_Indicator",
                    "Direct_control", "Direct_treated",
                    "Total Effect", "Indirect", "Direct", "Prop.Med")
  colnames(res) = c("Mean", "Median", "SD", "l-95% CI", "u-95% CI", "Bayes_p")

  res
}
cal.rr.effects <- function(outcome.pred)
{
  direct_control = outcome.pred[2,1,1,,] / outcome.pred[1,1,1,,]
  # treated: Y(1,M(1)) - Y(0,M(1))
  direct_treated = outcome.pred[2,2,1,,] / outcome.pred[1,2,1,,]
  # mediation effect: Y(t,M(1)) - Y(t,M(0))
  # control: Y(0,M(1)) - Y(0,M(0))
  indirect_control = outcome.pred[1,2,2,,] / outcome.pred[1,1,2,,]
  # treated: Y(1,M(1)) - Y(1,M(0))
  indirect_treated = outcome.pred[2,2,2,,] / outcome.pred[2,1,2,,]

  direct = (direct_control + direct_treated)/2
  indirect = (indirect_control + indirect_treated)/2

  total = direct * indirect
  # **************************************************
  pmed = direct*(indirect-1)/(total-1)

  res = rbind(
    c(mean(indirect_control), median(indirect_control), sd(indirect_control),
      quantile(indirect_control, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_control<1), mean(indirect_control>1))),

    c(mean(indirect_treated), median(indirect_treated), sd(indirect_treated),
      quantile(indirect_treated, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_treated<1), mean(indirect_treated>1))),

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

    c(mean(pmed), median(pmed), sd(pmed),
      quantile(pmed, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(pmed<0), mean(pmed>0)))
  ) # Bayes p-value: tail probability (see JMbayes), 2*min{pr(b<0), pr(b>0))}
  res[,1:5] = round(res[,1:5], digits=3)
  res[,6] = signif(res[,6], digits=2)
  res[8,1] = ifelse(res[8,1] < 0, 0, res[8,1] )
  rownames(res) = c("Indirect_control", "Indirect_treated", "Direct_control", "Direct_treated",
                    "Total Effect", "Indirect", "Direct", "Prop.Med")
  colnames(res) = c("Mean", "Median", "SD", "l-95% CI", "u-95% CI", "Bayes_p")

  res
}
cal.rr.effects.ind <- function(outcome.pred)
{
  direct_control = outcome.pred[2,1,1,,] / outcome.pred[1,1,1,,]
  # treated: Y(1,M(1)) - Y(0,M(1))
  direct_treated = outcome.pred[2,2,1,,] / outcome.pred[1,2,1,,]
  # mediation effect: Y(t,M(1)) - Y(t,M(0))
  # control: Y(0,M(1)) - Y(0,M(0))
  indirect_control = outcome.pred[1,2,2,,] / outcome.pred[1,1,2,,]
  # treated: Y(1,M(1)) - Y(1,M(0))
  indirect_treated = outcome.pred[2,2,2,,] / outcome.pred[2,1,2,,]

  direct = (direct_control + direct_treated)/2
  indirect = (indirect_control + indirect_treated)/2

  # **************************************************
  # treated: Y(1,M(0), I(m1)) - Y(1,M(0), I(m0))
  indirect_Im = outcome.pred[2,1,2,,] / outcome.pred[2,1,1,,]
  indirect = indirect * indirect_Im

  total = direct*indirect

  # **************************************************
  pmed = direct*(indirect-1)/(total-1)

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

    c(mean(pmed), median(pmed), sd(pmed),
      quantile(pmed, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(pmed<0), mean(pmed>0)))
  ) # Bayes p-value: tail probability (see JMbayes), 2*min{pr(b<0), pr(b>0))}
  res[,1:5] = round(res[,1:5], digits=3)
  res[,6] = signif(res[,6], digits=2)
  res[9,1] = ifelse(res[9,1] < 0, 0, res[9,1] )
  rownames(res) = c("Indirect_control", "Indirect_treated", "Indirect_Indicator",
                    "Direct_control", "Direct_treated",
                    "Total Effect", "Indirect", "Direct", "Prop.Med")
  colnames(res) = c("Mean", "Median", "SD", "l-95% CI", "u-95% CI", "Bayes_p")

  res
}

