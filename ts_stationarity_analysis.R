
#rm(list=ls())
pkages  <-c("urca","tseries","data.table","dplyr","broom","lubridate","car","lmtest","sandwich","orcutt", "ggplot2","zoo", "tidyr", "astsa")
for(pkage in pkages){
  if(!pkage %in% installed.packages())install.packages(pkage)
  
}
lapply(pkages, require, character.only= TRUE)

runStationarityTest <- function(var.name, var.series, max.lag){
  #'@author : Marvin Akuffo
  #'@description : This function conducts stationarity test on time series data of a variable using the following stationarityy test types: 
  #'               Augmented Dickey Fuller(ADF), Phillip Perron (PP), KPSS, "Zivot & Andrews Unit Root Test" and "DF-GLS" tests.
  #'               The function utilizes only 'urca' package in R i.e. only dependency package to test for stationarity.

  #'@param  :  var.name - The name of the variable
  #'@param  :  var.series - the time series of the variable 
  #'@param  :  lag.max - The highest number of lagged endogenous differenced variables to be included in the test regression

  #'@return  : stationarity output for all test methods in urca package. 
  
  num.rows= 3*max.lag
  if(!inherits(var.series, "ts")){stop("The data type of variable series must be a timeseries")}
  if(length(var.series)==0){stop("The varable timeseries cannot be empty or have zero length")}
  
  #Augmented-Dickey-Fuller Unit Root Test
  stationaritytest.result <- rbindlist(lapply(1:max.lag, function(s, var.name, var.series, structure.type="none", select.lags= "Fixed"){ 
    test.adf.none         <- ur.df(y= var.series, type = "none", lags=s, selectlags= "Fixed")
    return(list(Var.Name= var.name, Test.Type=test.adf.none@test.name,Stationarity.Type= "None", Variant.of.Test.Type="tau1",
                Lag= s, CritVal.1pct= test.adf.none@cval[1], CritVal.5pct= test.adf.none@cval[2], CritVal.10pct= test.adf.none@cval[3],
                test.stats.val= test.adf.none@teststat[1], p.value=NA))
    }, var.name= var.name, var.series= var.series, structure.type=structure.type, select.lags=select.lags), use.names=TRUE, fill=TRUE)
  stationaritytest.result <- rbindlist(list(stationaritytest.result, rbindlist(lapply(1:max.lag, function(s, var.name, var.series, type, selectlags){
    test.adf.drift      <- ur.df(y= var.series, type = "drift", lags=s, selectlags= "Fixed")
    out.adf.drift.tau2  <- list(Var.Name= var.name, Test.Type=test.adf.drift@test.name, Stationarity.Type= "Drift", Variant.of.Test.Type="tau2", Lag= s,  
                                CritVal.1pct= test.adf.drift@cval[1,1], CritVal.5pct= test.adf.drift@cval[1,2], CritVal.10pct= test.adf.drift@cval[1,3],
                                test.stats.val= test.adf.drift@teststat[1,1], p.value=NA)
    out.adf.drift.phi1  <- list(Var.Name= var.name, Test.Type=test.adf.drift@test.name, Stationarity.Type= "Drift", Variant.of.Test.Type="phi1", Lag= s,   
                                CritVal.1pct= test.adf.drift@cval[2,1], CritVal.5pct= test.adf.drift@cval[2,2], CritVal.10pct= test.adf.drift@cval[2,3],      
                                test.stats.val= test.adf.drift@teststat[1,2], p.value=NA)
    return(do.call(rbindlist, list(list(out.adf.drift.tau2, out.adf.drift.phi1), use.names=TRUE, fill=TRUE)))}, var.name= var.name, var.series= var.series, 
    type= type, selectlags=selectlags), use.names=TRUE, fill=TRUE)))
  stationaritytest.result <- rbindlist(list(stationaritytest.result, rbindlist(lapply(1:max.lag, function(s, var.name, var.series, type, selectlags){
    test.adf.trend        <- ur.df(y= var.series, type = "trend", lags=s, selectlags= "Fixed")
    test.adf.trend.tau3   <- list(Var.Name= var.name, Test.Type=test.adf.trend@test.name, Stationarity.Type= "Trend", Variant.of.Test.Type="tau3", Lag= s,  
                                  CritVal.1pct= test.adf.trend@cval[1,1], CritVal.5pct= test.adf.trend@cval[1,2], CritVal.10pct= test.adf.trend@cval[1,3],
                                  test.stats.val= test.adf.trend@teststat[1,1], p.value=NA)
    test.adf.trend.phi2   <- list(Var.Name= var.name, Test.Type= test.adf.trend@test.name, Stationarity.Type= "Trend", Variant.of.Test.Type="phi2", Lag= s,   
                                  CritVal.1pct= test.adf.trend@cval[2,1], CritVal.5pct= test.adf.trend@cval[2,2], CritVal.10pct= test.adf.trend@cval[2,3],       
                                  test.stats.val= test.adf.trend@teststat[1,2], p.value=NA)
    test.adf.trend.phi3   <- list(Var.Name= var.name, Test.Type=test.adf.trend@test.name, Stationarity.Type= "Trend", Variant.of.Test.Type="phi3", Lag= s,   
                                  CritVal.1pct= test.adf.trend@cval[3,1], CritVal.5pct= test.adf.trend@cval[3,2], CritVal.10pct= test.adf.trend@cval[3,3],       
                                  test.stats.val= test.adf.trend@teststat[1,3], p.value=NA)
    return(do.call(rbindlist, list(list(test.adf.trend.tau3, test.adf.trend.phi2, test.adf.trend.phi3), use.names=TRUE, fill=TRUE)))}, var.name= var.name, 
    var.series= var.series, type= type, selectlags=selectlags), use.names=TRUE, fill=TRUE)))
  
  #Phillips \& Perron Unit Root Test
  for(k in c("constant", "trend")){
    for(i in c("Z-alpha", "Z-tau")){
      for(j in 1:max.lag){
        test.pp<-ur.pp(var.series, type = i, model = k, lags= "short", use.lag= j)
        out.test.pp<-list(Var.Name= var.name, Test.Type= test.pp@test.name, Stationarity.Type= k, Variant.of.Test.Type= i, Lag= j, 
                          CritVal.1pct= test.pp@cval[1,1], CritVal.5pct= test.pp@cval[1,2], CritVal.10pct= test.pp@cval[1,3],
                          test.stats.val= test.pp@teststat[1], p.value=NA)
        stationaritytest.result <- dplyr::bind_rows(stationaritytest.result, out.test.pp)
      }
    }
  }
  
  #Kwiatkowski et al. Unit Root Test
  for(k in c("mu", "tau")){
    for(j in 1:max.lag){
      test.kpss      <-ur.kpss(var.series, type = k, lags = "short", use.lag = j)
      out.test.kpss  <-list(Var.Name= var.name, Test.Type=test.kpss@test.name, Stationarity.Type= "Deterministic", Variant.of.Test.Type=k, Lag= j, CritVal.1pct= test.kpss@cval[1,1], 
                            CritVal.5pct= test.kpss@cval[1,2], CritVal.10pct= test.kpss@cval[1,3], test.stats.val= test.kpss@teststat[1], p.value=NA)
      stationaritytest.result <- dplyr::bind_rows(stationaritytest.result, out.test.kpss)
    }
  }
  
  #Elliott, Rothenberg \& Stock Unit Root Test
  for(k in c("constant", "trend")){
    for(i in c("P-test", "DF-GLS")){
      for(j in 1:max.lag){
        test.df.gls      <-ur.ers(var.series, type = i , model =k , lag.max= j)
        out.test.df_gls  <-list(Var.Name= var.name, Test.Type= test.df.gls@test.name, Stationarity.Type= k, Variant.of.Test.Type= i, Lag= j, CritVal.1pct= test.df.gls@cval[1,1], 
                                CritVal.5pct= test.df.gls@cval[1,2], CritVal.10pct= test.df.gls@cval[1,3], test.stats.val= test.df.gls@teststat[1], p.value=NA)
        stationaritytest.result <- dplyr::bind_rows(stationaritytest.result, out.test.kpss)
      }
    }
  }
  
  #Zivot \& Andrews Unit Root Test
  for(i in c("intercept", "trend", "both")){
    for(j in 1:max.lag){
      test.za<-ur.za(var.series, model = i, lag=j)
      out.test.za  <-list(Var.Name= var.name, Test.Type=test.za@test.name,Stationarity.Type= i, Variant.of.Test.Type= "t-stats=min{coeffs}", Lag= j, CritVal.1pct= test.za@cval[1], 
                          CritVal.5pct= test.za@cval[2], CritVal.10pct= test.za@cval[3], test.stats.val= test.kpss@teststat[1], p.value=NA)
      stationaritytest.result <- dplyr::bind_rows(stationaritytest.result, out.test.za)
    }
  }
  
  #write to a csv file
  if(!dir.exists("./output/"))dir.create("output")
  write.csv(stationaritytest.result, file= paste0("./output/",var.name,"_StationarityTest.csv"))
  
  return(stationaritytest.result)
}

#Example
ts.data <- ts(1000*rnorm(100), frequency = 4, start = c(2000, 1)) 
var.series <- ts.data
model.name = "DCM"
var.name= "SP500"
runStationarityTest("SP500", var.series, 4)
