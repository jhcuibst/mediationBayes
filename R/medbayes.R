
#'
#' @name medbayes
#' @title Run Bayesian mediation model
#'
#' @param model.m A brms model for the mediator
#' @param model.y A brms model for the outcome
#' @param dat.new New data for prediction
#' @return Posterior predictions (first column)
#' @export
medbayes <- function(model.m = model.m,
                     model.y = model.y,
                     treat = "treatment",
                     mediator = "mediator",
                     ind_mediator = NULL,
                     outcome = "outcome",
                     control.value = 0,
                     treat.value = 1){

  zi.mediator = grepl("zero", family(model.m)$family) | grepl("hurdle", family(model.m)$family)
  zi.outcome = grepl("zero", family(model.y)$family) | grepl("hurdle", family(model.y)$family)

  if(zi.outcome){
    result <- medbayes_ziy(model.m=model.m, model.y=model.y, treat=treat,
                           mediator = mediator, ind_mediator = ind_mediator,
                           outcome = outcome,
                           control.value = control.value, treat.value = treat.value)
  }else{
    result <- medbayes_zim(model.m=model.m, model.y=model.y, treat=treat,
                           mediator = mediator, ind_mediator = ind_mediator,
                           outcome = outcome,
                           control.value = control.value, treat.value = treat.value)
  }
  return(result)
}
