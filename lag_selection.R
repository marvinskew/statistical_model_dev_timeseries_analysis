
lagSelection<-function(y, var.name, type = c("none", "drift", "trend"), lags=1){
  #'@Descrption: This function computes both the AIC and BIC for each AR model of 1,..., lags order
  #'@param : y - a numeric vector or time series.
  #'@param : type - Test type, either "none", "drift" or "trend"
  #'@param : lags - Number of lags for endogenous variable to be included.
  
  #'@return : dataframe with model lag lengths and corresponding AIC and BIC.
  #'
  lags <- lags + 1
  z <- diff(y)
  n <- length(z)
  x <- embed(z, lags)
  z.diff <- x[, 1]
  z.lag.1 <- y[lags:n]
  tt <- lags:n
  result.aic <- result.bic <- c()
  result <- list()
  res.names  <- c("VariableName", "Model.Type", "Lag", "AIC", "BIC")
  res        <- data.frame(matrix(NA, nrow=lags, ncol=length(res.names)))
  for(i in 2:lags){
    j= i-1; z.diff.lag = x[, 2:i]
    if (type == "none"){result[[j]]  <- lm(z.diff ~ z.lag.1 - 1 + z.diff.lag)}
    if (type == "drift"){result[[j]] <- lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)} 
    if (type == "trend"){result[[j]] <- lm(z.diff ~ z.lag.1 + 1 + tt + z.diff.lag)}
    res.aic[j] <- AIC(result[[j]])
    res.bic[j] <- BIC(result[[j]])
    res[j, ]<-list(model.name, var.name, type, j, res.aic[j], res.bic[j])
  }
  return(na.omit(res))
}


