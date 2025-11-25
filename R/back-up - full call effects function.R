# For mediator ----
cal.rd.effects.ind  <- function(outcome.pred, ind_mediator)
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
  res[,6] = round(res[,6], digits=4)
  res[9,1] = ifelse(res[9,1] < 0, 0, res[9,1] )
  rownames(res) = c("Indirect_control", "Indirect_treated", "Indirect_Indicator",
                    "Direct_control", "Direct_treated",
                    "Total Effect", "Indirect", "Direct", "Prop.Med")
  colnames(res) = c("Mean", "Median", "SD", "l-95% CI", "u-95% CI", "Bayes_p")
  if(is.null(ind_mediator)){
    res = res[rownames(res) != "Indirect_Indicator", ]
  }else{
    res = res
  }
  res
}
cal.rr.effects.ind <- function(outcome.pred, ind_mediator)
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
  res[,6] = round(res[,6], digits=4)
  res[9,1] = ifelse(res[9,1] < 0, 0, res[9,1] )

  rownames(res) = c("Indirect_control", "Indirect_treated", "Indirect_Indicator",
                    "Direct_control", "Direct_treated",
                    "Total Effect", "Indirect", "Direct", "Prop.Med")
  colnames(res) = c("Mean", "Median", "SD", "l-95% CI", "u-95% CI", "Bayes_p")
  if(is.null(ind_mediator)){
    res = res[rownames(res) != "Indirect_Indicator", ]
  }else{
    res = res
  }
  res
}

# for outcome
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


