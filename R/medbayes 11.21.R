
library(brms)
# **************************************************
medbayes2 <- function(model.m = hnb.m, model.y = hnb.y,
                      treat = "treatment", mediator = "mediator", ind_mediator = NULL, outcome = "outcome",
                      control.value = 0, treat.value = 1,
                      logM = FALSE )
{

  value = c(control.value, treat.value)

  dat.new = model.m$data

  predict.ms <- array(NA, dim = c(2, nrow(dat.new), 1))
  for(i in 1:length(value)){
    dat.new[, treat] = value[i]
    predict.ms[i,,] = 0
    predict.ms[i,,] = predict.glm(fitm2, newdata = dat.new, type = "response")
  }

  dat.y = model.y$data

  outcome.linpred.mu <- array(NA, dim = c(2, 2,ndraws(model.y), nrow(dat.y)) )
  outcome.linpred.zi <- array(NA, dim = c(2, 2,ndraws(model.y), nrow(dat.y)) )

  if(grepl("zero", family(model.y)$family)) depar.outcome = c("mu", "zi")
  if(grepl("hurdle", family(model.y)$family)) depar.outcome = c("mu", "hu")

  for(i in 1:length(value)){
    for(j in 1:length(value)){
      dat = dat.y %>% mutate(x = value[i], m = predict.ms[j,,])
      outcome.linpred.mu[i,j,,] = posterior_linpred(object = fity, newdata = dat, dpar = depar.outcome[1])
      outcome.linpred.zi[i,j,,] = posterior_linpred(object = fity, newdata = dat, dpar = depar.outcome[2])
    }
  }

  outcome.pred.mu  = exp(outcome.linpred.mu)
  outcome.pred.zi  = exp(outcome.linpred.zi)/(1+exp(outcome.linpred.zi))

  res.mu.zi = cal.rr.effects.y(outcome.pred.mu, outcome.pred.zi)
  res.rr <- list(effects.mu.rr = res.mu.zi)

  rst <- list(effects.rr = res.rr, outcome.linpred.mu=outcome.linpred.mu,
              outcome.linpred.zi=outcome.linpred.zi,
              outcome.pred.mu= outcome.pred.mu,  outcome.pred.zi =  outcome.pred.zi)
  rst
}

cal.rr.effects.y <- function(outcome.pred = outcome.pred.mu, outcome.pred.zi = outcome.pred.zi)
{
  direct_control = (outcome.pred [2,1,,] / outcome.pred [1,1,,])
  direct_treated = (outcome.pred[2,2,,] / outcome.pred[1,2,,])
  direct_Im_c = (1-outcome.pred.zi[2,1,,]) / (1-outcome.pred.zi[1,1,,])
  direct_Im_t = (1-outcome.pred.zi[2,2,,]) / (1-outcome.pred.zi[1,2,,])

  indirect_control = outcome.pred[1,2,,] / outcome.pred[1,1,,]
  indirect_treated = outcome.pred[2,2,,] / outcome.pred[2,1,,]
  indirect_Im_t = (1-outcome.pred.zi[2,2,,]) / (1-outcome.pred.zi[2,1,,])
  indirect_Im_c = (1-outcome.pred.zi[1,2,,]) / (1-outcome.pred.zi[1,1,,])

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
  pmed = direct.c_total*(indirect_treated*indirect_Im_t-1)/(total-1)


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
      2*min(mean(d_avg*i_avg<1), mean(d_avg*i_avg>1)))

  ) # Bayes p-value: tail probability (see JMbayes), 2*min{pr(b<0), pr(b>0))}
  res[,1:5] = round(res[,1:5], digits=3)
  res[,6] = signif(res[,6], digits=2)
  res[9,1] = ifelse(res[9,1] < 0, 0, res[9,1] )

  colnames(res) = c("Mean", "Median", "SD", "l-95% CI", "u-95% CI", "Bayes_p")
  rownames(res) = c("NDE_NZ.C",	"NDE_NZ.T",	"NDE_Z.C",	"NDE_Z.T", "NIE_NZ.C",	"NIE_NZ.T",	"NIE_Z.C",	"NIE_Z.T",
                    "NDE.C",	"NDE.T",	"NIE.C",	"NIE.T",
                    "NDE.AVG",	"NIE_NZ.AVG",	"NIE_Z.AVG", "TOTAL",	"TOTAL_AVG")
  res
}
cal.rd.effects.y <- function(outcome.pred = outcome.pred.mu, outcome.pred.zi = outcome.pred.zi)
{
  direct_control = (outcome.pred [2,1,,] - outcome.pred [1,1,,])
  direct_treated = (outcome.pred[2,2,,] - outcome.pred[1,2,,])
  direct_Im_c = (1-outcome.pred.zi[2,1,,]) - (1-outcome.pred.zi[1,1,,])
  direct_Im_t = (1-outcome.pred.zi[2,2,,]) - (1-outcome.pred.zi[1,2,,])

  indirect_control = outcome.pred[1,2,,] - outcome.pred[1,1,,]
  indirect_treated = outcome.pred[2,2,,] - outcome.pred[2,1,,]
  indirect_Im_t = (1-outcome.pred.zi[2,2,,]) / (1-outcome.pred.zi[2,1,,])
  indirect_Im_c = (1-outcome.pred.zi[1,2,,]) / (1-outcome.pred.zi[1,1,,])

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
  pmed =(indirect_treated+indirect_Im_t)/(total)


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
      2*min(mean(d_avg*i_avg<1), mean(d_avg*i_avg>1)))

  ) # Bayes p-value: tail probability (see JMbayes), 2*min{pr(b<0), pr(b>0))}
  res[,1:5] = round(res[,1:5], digits=3)
  res[,6] = signif(res[,6], digits=2)
  res[9,1] = ifelse(res[9,1] < 0, 0, res[9,1] )

  colnames(res) = c("Mean", "Median", "SD", "l-95% CI", "u-95% CI", "Bayes_p")
  rownames(res) = c("NDE_NZ.C",	"NDE_NZ.T",	"NDE_Z.C",	"NDE_Z.T", "NIE_NZ.C",	"NIE_NZ.T",	"NIE_Z.C",	"NIE_Z.T",
                    "NDE.C",	"NDE.T",	"NIE.C",	"NIE.T",
                    "NDE.AVG",	"NIE_NZ.AVG",	"NIE_Z.AVG", "TOTAL",	"TOTAL_AVG")
  res
}



