#'
#' @name medbayes
#' @title Run Bayesian mediation model
#'
#' @param model.m A brms model for the mediator
#' @param model.y A brms model for the outcome
#' @param dat.new New data for prediction
#' @return Posterior predictions (first column)
#' @export
medbayes.y <- function(model.m = model.m, model.y = model.y, treat = "treatment",
                     mediator = "mediator", ind_mediator = NULL,
                     outcome = "outcome", # ind_outcome = NULL,
                     control.value = 0, treat.value = 1)
{
  library(brms)
  # Check whether package 'brms' is installed
  if (!requireNamespace("brms", quietly = TRUE))
    stop("Package 'brms' is required but not installed.")
  message("Mediation analysis function is running...")

  value = c(control.value, treat.value)

  # *********************************************************************
  dat.new = model.m$data
  predict.ms <- array(NA, dim = c(2, ndraws(model.m), nrow(dat.new)))
  for(i in 1:length(value)){
    dat.new[, treat] = value[i]
    predict.ms[i,,] = 0
    predict.ms[i,,] = posterior_epred(model.m, newdata = dat.new, dpar = NULL)
  }

  dat.y = model.y$data
  depar.outcome = NULL
  if(grepl("zero", family(model.y)$family)) depar.outcome = c("mu", "zi")
  if(grepl("hurdle", family(model.y)$family)) depar.outcome = c("mu", "hu")
  zi.outcome = grepl("zero", family(model.y)$family) | grepl("hurdle", family(model.y)$family)

  # pred_y without mediator(s) effect ----
  # predict.y.cov <- array(NA, dim = c(2, ndraws(model.y), nrow(dat.y)) )
  predict.y.cov.mu <- array(NA, dim = c(2, ndraws(model.y), nrow(dat.y)) )
  predict.y.cov.zi <- array(NA, dim = c(2, ndraws(model.y), nrow(dat.y)) )

  for(i in 1:length(value)){
    dat.y.temp <- dat.y
    dat.y.temp[, treat] <- value[i]
    dat.y.temp[, mediator] = 0

    predict.y.cov.mu[i,,] = posterior_linpred(model.y, newdata = dat.y.temp, dpar = depar.outcome[1])
    predict.y.cov.zi[i,,] = posterior_linpred(model.y, newdata = dat.y.temp, dpar = depar.outcome[2])
  }

  # pred_y for treat and mediator(s) effect ----
  ef_y = as_draws_df(model.y)
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

  ## check interaction between treatment and mediator(s)
  INT.xm <- c(paste("b_", treat, ":", mediator, sep="") %in% colnames(as.matrix(model.y)),
              paste("b_", mediator, ":", treat, sep="") %in% colnames(as.matrix(model.y))  )

  INT.xm_hu <- c(paste("b_hu_", treat, ":", mediator, sep="") %in% colnames(as.matrix(model.y)) ,
                 paste("b_hu_", mediator, ":", treat, sep="") %in% colnames(as.matrix(model.y))  )
  INT.xm_zi <- c(paste("b_zi_", treat, ":", mediator, sep="") %in% colnames(as.matrix(model.y)) ,
                 paste("b_zi_", mediator, ":", treat, sep="") %in% colnames(as.matrix(model.y))  )

  coef.intxm = c(paste("b_", treat, ":", mediator, sep=""), paste("b_", mediator, ":", treat, sep="")  )
  coef.intxm = coef.intxm[INT.xm]
  int_of_xm <- any(coef.intxm %in% colnames(as.matrix(model.y)))

  # ******** check if interaction in zi/hu y-model
  coef.inthuxm = c(paste("b_hu_", treat, ":", mediator, sep=""), paste("b_hu_", mediator, ":", treat, sep="")  )
  coef.intzixm = c(paste("b_zi_", treat, ":", mediator, sep=""), paste("b_zi_", mediator, ":", treat, sep="")  )
  coef.intxmhu = coef.intxm[INT.xm_hu]
  coef.intxmzi = coef.intxm[INT.xm_zi]

  int_of_xm_hu <- any(coef.intxmhu %in% colnames(as.matrix(model.y)))
  int_of_xm_zi <- any(coef.intxmzi %in% colnames(as.matrix(model.y)))

  # ********
  if(int_of_xm) {
    bxm = as.matrix(ef_y)[,coef.intxm]
  } else(
    bxm = 0
  )

  if(int_of_xm_hu){
    bxm_zi = as.matrix(ef_y)[,coef.intxmhu]
  } else if(int_of_xm_zi){
    bxm_zi = as.matrix(ef_y)[,coef.intxmzi]
  } else{
    bxm_zi = 0
  }

  outcome.linpred.mu = array(NA, dim = c(2, 2, 2, ndraws(model.y), nrow(dat.y)))
  outcome.linpred.zi = array(NA, dim = c(2, 2, 2, ndraws(model.y), nrow(dat.y)))

  # ****************************************************************************;
  # calculate pred mu of outcome
  calc_linpred.mu <- function(i, j, r, bm, bxm, int_xm_zi){
    if(any(INT.xm)){
      int_xm_zi = 1
    } else {
      int_xm_zi=0
      bxm=0
    }

    as.numeric(
      as.matrix(predict.ms[j,,]*bm) + as.matrix(bxm*int_xm_zi*value[i]*predict.ms[j,,]) +
        predict.y.cov.mu[i,,] )
  }

  for(i in 1:length(value)){
    for(j in 1:2){
      for(r in 1:2){
        outcome.linpred.mu[i,j,r,,] =  calc_linpred.mu(i, j, r, bm, bxm, int_xm_zi)
      }
    }
  }

  # ****************************************************************************;
  # calculate pred zi of outcome
  calc_linpred.zi <- function(i, j, r, bm, bxm, int_of_xm_hu_zi ){

    if(any(int_of_xm_hu, int_of_xm_zi)){
      int_xm_hu_zi = 1
    } else {
      int_xm_hu_zi=0
      bxm_zi=0}

    as.numeric(
      as.matrix(predict.ms[j,,]*b_zi_m) + as.matrix(bxm_zi*int_xm_hu_zi*value[i]*predict.ms[j,,]) +
        predict.y.cov.zi[i,,] )
  }

  for(i in 1:length(value)){
    for(j in 1:2){
      for(r in 1:2){
        outcome.linpred.zi[i,j,r,,] =  calc_linpred.zi(i, j, r, bm, bxm, int_of_xm_hu_zi )
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

  outcome.pred.mu  = exp(outcome.linpred.mu)
  outcome.pred.zi  = exp(outcome.linpred.zi)/(1+exp(outcome.linpred.zi))
  outcome.pred.overall  = outcome.pred.mu*(1-outcome.pred.zi)

  res.rr = cal.rr.effects.y(outcome.pred.mu, outcome.pred.zi)
  res.rd = cal.rd.effects.y(outcome.pred.mu, outcome.pred.zi)
  rst <- list(effects.rd = res.rd, effects.rr = res.rr, outcome.linpred.mu=outcome.linpred.mu,
                  outcome.linpred.zi=outcome.linpred.zi,
                  outcome.pred.mu= outcome.pred.mu,  outcome.pred.zi =  outcome.pred.zi)

  rst
}

cal.rr.effects.y <- function(outcome.pred = outcome.pred.mu, outcome.pred.zi = outcome.pred.zi)
{
  direct_control = (outcome.pred [2,1,2,,] / outcome.pred [1,1,2,,])
  direct_treated = (outcome.pred[2,2,2,,] / outcome.pred[1,2,2,,])
  direct_Im_c = (1-outcome.pred.zi[2,1,2,,]) / (1-outcome.pred.zi[1,1,2,,])
  direct_Im_t = (1-outcome.pred.zi[2,2,2,,]) / (1-outcome.pred.zi[1,2,2,,])

  indirect_control = outcome.pred[1,2,2,,] / outcome.pred[1,1,2,,]
  indirect_treated = outcome.pred[2,2,2,,] / outcome.pred[2,1,2,,]
  indirect_Im_t = (1-outcome.pred.zi[2,2,2,,]) / (1-outcome.pred.zi[2,1,2,,])
  indirect_Im_c = (1-outcome.pred.zi[1,2,2,,]) / (1-outcome.pred.zi[1,1,2,,])

  direct.c_total = median(direct_control) * median(direct_Im_c)
  direct.t_total = median(direct_treated) * median(direct_Im_t)
  indirect.c_total = mean(indirect_control)*mean(indirect_Im_c)
  indirect.t_total = mean(indirect_treated)*mean(indirect_Im_t)

  direct_avg = (direct.c_total + direct.t_total)/2
  indirect_nz_avg = (indirect_control + indirect_treated)/2
  indirect_z_avg  = (indirect_Im_t + indirect_Im_c)/2

  total = direct.c_total*indirect.t_total
  total_avg = direct_avg*indirect_nz_avg*indirect_z_avg

  d_avg = (direct_control*direct_Im_c + direct_treated*direct_Im_t)/2
  i_avg = (indirect_control*indirect_Im_c + indirect_treated*indirect_Im_t)/2

  # **************************************************
  pmed = direct.c_total*(indirect_treated*indirect_Im_t-1)/(direct.c_total*indirect_treated*indirect_Im_t-1)
  pmed.mean = mean(direct.c_total)*(mean(indirect_treated*indirect_Im_t)-1)/(mean(direct.c_total)*mean(indirect_treated*indirect_Im_t)-1)
  pmed.med = median(direct.c_total)*(median(indirect_treated*indirect_Im_t)-1)/(median(direct.c_total)*median(indirect_treated*indirect_Im_t)-1)

  res = rbind(
    c(mean(direct_control), median(direct_control), sd(direct_control),
      quantile(direct_control, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct_control<1), mean(direct_control>1))),
    c(mean(direct_treated), median(direct_treated), sd(direct_treated),
      quantile(direct_treated, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct_treated<1), mean(direct_treated>1))),
    c(mean(direct_Im_c), median(direct_Im_c), sd(direct_Im_c),
      quantile(direct_Im_c, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct_Im_c<1), mean(direct_Im_c>1))),
    c(mean(direct_Im_t), median(direct_Im_t), sd(direct_Im_t),
      quantile(direct_Im_t, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct_Im_t<1), mean(direct_Im_t>1))),
    # **************************************************

    c(mean(indirect_control), median(indirect_control), sd(indirect_control),
      quantile(indirect_control, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_control<1), mean(indirect_control>1))),
    c(mean(indirect_treated), median(indirect_treated), sd(indirect_treated),
      quantile(indirect_treated, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_treated<1), mean(indirect_treated>1))),
    c(mean(indirect_Im_c), median(indirect_Im_c), sd(indirect_Im_c),
      quantile(indirect_Im_c, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_Im_c<1), mean(indirect_Im_c>1))),
    c(mean(indirect_Im_t), median(indirect_Im_t), sd(indirect_Im_t),
      quantile(indirect_Im_t, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_Im_t<1), mean(indirect_Im_t>1))),
    # **************************************************

    c(mean(direct.c_total), median(direct.c_total), sd(direct_control*direct_Im_c),
      quantile(direct_control*direct_Im_c, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct_control*direct_Im_c<1), mean(direct_control*direct_Im_c>1))),
    c(mean(direct.t_total), median(direct.t_total), sd(direct_treated*direct_Im_t),
      quantile(direct_treated*direct_Im_t, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct_treated*direct_Im_t<1), mean(direct_treated*direct_Im_t>1))),
    c(mean(indirect.c_total), median(indirect.c_total), sd(indirect_control*indirect_Im_c),
      quantile(indirect_control*indirect_Im_c, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_control*indirect_Im_c<1), mean(indirect_control*indirect_Im_c>1))),
    c(mean(indirect.t_total), median(indirect.t_total), sd(indirect_treated*indirect_Im_t),
      quantile(indirect_treated*indirect_Im_t, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_treated*indirect_Im_t<1), mean(indirect_treated*indirect_Im_t>1))),
    # **************************************************

    c(mean(direct_avg), median(direct_avg), sd((direct_control*direct_Im_c+  direct_treated*direct_Im_t)/2),
      quantile((direct_control*direct_Im_c+  direct_treated*direct_Im_t)/2, probs=c(0.025,0.975), na.rm = T),
      2*min(mean((direct_control*direct_Im_c+  direct_treated*direct_Im_t)/2<1),
            mean((direct_control*direct_Im_c+  direct_treated*direct_Im_t)/2>1))),
    c(mean(indirect_nz_avg), median(indirect_nz_avg), sd(indirect_nz_avg),
      quantile(indirect_nz_avg, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_nz_avg<1), mean(indirect_nz_avg>1))),
    c(mean(indirect_z_avg), median(indirect_z_avg), sd(indirect_z_avg),
      quantile(indirect_z_avg, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_z_avg<1), mean(indirect_z_avg>1))),
    # **************************************************

    c(mean(total), median(total), sd(direct_control*direct_Im_c*indirect_treated*indirect_Im_t),
      quantile(direct_control*direct_Im_c*indirect_treated*indirect_Im_t, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct_control*direct_Im_c*indirect_treated*indirect_Im_t<1),
            mean(direct_control*direct_Im_c*indirect_treated*indirect_Im_t>1))),
    c(mean(total_avg), median(total_avg), sd(d_avg*i_avg),
      quantile(d_avg*i_avg, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(d_avg*i_avg<1), mean(d_avg*i_avg>1))),

    c(mean(pmed.mean), median(pmed.med), sd(pmed),
      quantile(pmed, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(pmed<1), mean(pmed>1)))

  ) # Bayes p-value: tail probability (see JMbayes), 2*min{pr(b<0), pr(b>0))}
  res[,1:5] = round(res[,1:5], digits=3)
  res[,6] = signif(res[,6], digits=4)
  res[9,1] = ifelse(res[9,1] < 0, 0, res[9,1] )

  colnames(res) = c("Mean", "Median", "SD", "l-95% CI", "u-95% CI", "Bayes_p")
  rownames(res) = c("NDE_NZ.C",	"NDE_NZ.T",	"NDE_Z.C",	"NDE_Z.T", "NIE_NZ.C",	"NIE_NZ.T",	"NIE_Z.C",	"NIE_Z.T",
                    "NDE.C",	"NDE.T",	"NIE.C",	"NIE.T",
                    "NDE.AVG",	"NIE_NZ.AVG",	"NIE_Z.AVG", "TOTAL",	"TOTAL_AVG", "PMed")
  res
}
cal.rd.effects.y <- function(outcome.pred = outcome.pred.mu, outcome.pred.zi = outcome.pred.zi)
{
  direct_control = (outcome.pred [2,1,2,,] - outcome.pred [1,1,2,,])
  direct_treated = (outcome.pred[2,2,2,,] - outcome.pred[1,2,2,,])
  direct_Im_c = (1-outcome.pred.zi[2,1,2,,]) - (1-outcome.pred.zi[1,1,2,,])
  direct_Im_t = (1-outcome.pred.zi[2,2,2,,]) - (1-outcome.pred.zi[1,2,2,,])

  indirect_control = outcome.pred[1,2,2,,] - outcome.pred[1,1,2,,]
  indirect_treated = outcome.pred[2,2,2,,] - outcome.pred[2,1,2,,]
  indirect_Im_t = (1-outcome.pred.zi[2,2,2,,]) / (1-outcome.pred.zi[2,1,2,,])
  indirect_Im_c = (1-outcome.pred.zi[1,2,2,,]) / (1-outcome.pred.zi[1,1,2,,])

  direct.c_total = median(direct_control) + median(direct_Im_c)
  direct.t_total = median(direct_treated) + median(direct_Im_t)
  indirect.c_total = mean(indirect_control)+ mean(indirect_Im_c)
  indirect.t_total = mean(indirect_treated)+ mean(indirect_Im_t)

  direct_avg = (direct.c_total + direct.t_total)/2
  indirect_nz_avg = (indirect_control + indirect_treated)/2
  indirect_z_avg  = (indirect_Im_t + indirect_Im_c)/2

  total = direct.c_total+indirect.t_total
  total_avg = direct_avg+indirect_nz_avg+indirect_z_avg

  d_avg = (direct_control+direct_Im_c + direct_treated+direct_Im_t)/2
  i_avg = (indirect_control+indirect_Im_c + indirect_treated+indirect_Im_t)/2

  # **************************************************
  pmed = indirect_treated*indirect_Im_t/(direct.c_total+indirect_treated+indirect_Im_t)
  pmed.mean = mean(indirect_treated)*mean(indirect_Im_t)/mean(direct.c_total+indirect_treated+indirect_Im_t)
  pmed.med = median(indirect_treated)*median(indirect_Im_t)/median(direct.c_total+indirect_treated+indirect_Im_t)

  res = rbind(
    c(mean(direct_control), median(direct_control), sd(direct_control),
      quantile(direct_control, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct_control<1), mean(direct_control>1))),
    c(mean(direct_treated), median(direct_treated), sd(direct_treated),
      quantile(direct_treated, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct_treated<1), mean(direct_treated>1))),
    c(mean(direct_Im_c), median(direct_Im_c), sd(direct_Im_c),
      quantile(direct_Im_c, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct_Im_c<1), mean(direct_Im_c>1))),
    c(mean(direct_Im_t), median(direct_Im_t), sd(direct_Im_t),
      quantile(direct_Im_t, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct_Im_t<1), mean(direct_Im_t>1))),
    # **************************************************

    c(mean(indirect_control), median(indirect_control), sd(indirect_control),
      quantile(indirect_control, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_control<1), mean(indirect_control>1))),
    c(mean(indirect_treated), median(indirect_treated), sd(indirect_treated),
      quantile(indirect_treated, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_treated<1), mean(indirect_treated>1))),
    c(mean(indirect_Im_c), median(indirect_Im_c), sd(indirect_Im_c),
      quantile(indirect_Im_c, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_Im_c<1), mean(indirect_Im_c>1))),
    c(mean(indirect_Im_t), median(indirect_Im_t), sd(indirect_Im_t),
      quantile(indirect_Im_t, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_Im_t<1), mean(indirect_Im_t>1))),
    # **************************************************

    c(mean(direct.c_total), median(direct.c_total), sd(direct_control*direct_Im_c),
      quantile(direct_control*direct_Im_c, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct_control*direct_Im_c<1), mean(direct_control*direct_Im_c>1))),
    c(mean(direct.t_total), median(direct.t_total), sd(direct_treated*direct_Im_t),
      quantile(direct_treated*direct_Im_t, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct_treated*direct_Im_t<1), mean(direct_treated*direct_Im_t>1))),
    c(mean(indirect.c_total), median(indirect.c_total), sd(indirect_control*indirect_Im_c),
      quantile(indirect_control*indirect_Im_c, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_control*indirect_Im_c<1), mean(indirect_control*indirect_Im_c>1))),
    c(mean(indirect.t_total), median(indirect.t_total), sd(indirect_treated*indirect_Im_t),
      quantile(indirect_treated*indirect_Im_t, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_treated*indirect_Im_t<1), mean(indirect_treated*indirect_Im_t>1))),
    # **************************************************

    c(mean(direct_avg), median(direct_avg), sd((direct_control*direct_Im_c+  direct_treated*direct_Im_t)/2),
      quantile((direct_control*direct_Im_c+  direct_treated*direct_Im_t)/2, probs=c(0.025,0.975), na.rm = T),
      2*min(mean((direct_control*direct_Im_c+  direct_treated*direct_Im_t)/2<1),
            mean((direct_control*direct_Im_c+  direct_treated*direct_Im_t)/2>1))),
    c(mean(indirect_nz_avg), median(indirect_nz_avg), sd(indirect_nz_avg),
      quantile(indirect_nz_avg, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_nz_avg<1), mean(indirect_nz_avg>1))),
    c(mean(indirect_z_avg), median(indirect_z_avg), sd(indirect_z_avg),
      quantile(indirect_z_avg, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_z_avg<1), mean(indirect_z_avg>1))),
    # **************************************************

    c(mean(total), median(total), sd(direct_control*direct_Im_c*indirect_treated*indirect_Im_t),
      quantile(direct_control*direct_Im_c*indirect_treated*indirect_Im_t, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(direct_control*direct_Im_c*indirect_treated*indirect_Im_t<1),
            mean(direct_control*direct_Im_c*indirect_treated*indirect_Im_t>1))),
    c(mean(total_avg), median(total_avg), sd(d_avg*i_avg),
      quantile(d_avg*i_avg, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(d_avg*i_avg<1), mean(d_avg*i_avg>1))),

    c(mean(pmed), median(pmed), sd(pmed),
      quantile(pmed, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(pmed<1), mean(pmed>1)))

  ) # Bayes p-value: tail probability (see JMbayes), 2*min{pr(b<0), pr(b>0))}
  res[,1:5] = round(res[,1:5], digits=3)
  res[,6] = signif(res[,6], digits=4)
  res[9,1] = ifelse(res[9,1] < 0, 0, res[9,1] )

  colnames(res) = c("Mean", "Median", "SD", "l-95% CI", "u-95% CI", "Bayes_p")
  rownames(res) = c("NDE_NZ.C",	"NDE_NZ.T",	"NDE_Z.C",	"NDE_Z.T", "NIE_NZ.C",	"NIE_NZ.T",	"NIE_Z.C",	"NIE_Z.T",
                    "NDE.C",	"NDE.T",	"NIE.C",	"NIE.T",
                    "NDE.AVG",	"NIE_NZ.AVG",	"NIE_Z.AVG", "TOTAL",	"TOTAL_AVG", "PMed")
  res
}



