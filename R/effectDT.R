# Adapted from Thorson and Kristensen 2024

# Change marginaleffects options to define `custom_tmb` class
options("marginaleffects_model_classes" = "custom_tmb")
quant = function(x) seq(min(x),max(x),length=100)

# Function to get coefficients for TMB model
get_coef.custom_tmb = function(model, param, ...){
  out = model$parhat[[param]]
  names(out) = rep(param, length(out))
  return(out)
}

# Function to get variance-covariance for TMB model
get_vcov.custom_tmb = function(model, param, ...){
  rows = which( names(model$opt$par) == param )
  array( model$sd$cov.fixed[rows,rows],
         dim = rep(length(rows),2),
         dimnames = list(rep(param,length(rows)),rep(param,length(rows))) )
}

# get_vcov(model = fit, param = "beta_k")

# Function to change coefficients for TMB model
set_coef.custom_tmb = function(model, newpar, param, ...){
  model$parhat[[param]] <- newpar
  return(model)
}

# set_coef(model = fit, newpar = rep(0, 11), param = "beta_k")

# Function to get predictions when changing coefficients
get_predict.custom_tmb = function(model, newdata, param, center=FALSE, ...){
  Xpred <- predict(model$pref_mod, newdata = newdata, type = "lpmatrix")
  beta_k = get_coef.custom_tmb(model, param)
  yhat_i = Xpred %*% beta_k
  if(center==TRUE) yhat_i = yhat_i - mean(yhat_i)
  out = data.frame( rowid=seq_along(yhat_i[,1]), estimate=yhat_i )
  return(out)
}

