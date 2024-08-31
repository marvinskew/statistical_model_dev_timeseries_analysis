
#rm(list=ls())
#rm(list=ls())
pkages  <-c("urca","tseries","data.table","dplyr","broom","lubridate","car","lmtest","sandwich","orcutt", "ggplot2","zoo", "tidyr", "astsa")
for(pkage in pkages){
  if(!pkage %in% installed.packages())install.packages(pkage)
  
}
lapply(pkages, require, character.only= TRUE)

cointegrationTest<-function(lm.obj, data, max.lag, method=c("Engle-Granger","Johansen"), johansen.type= c("eigen","trace"), plot.it=TRUE){
  #'@author : Marvin Akuffo
  #'@description : This function test for cointegration - the long run equilibrium of variables used for running the model.
  #'@param : lm.obj - linear regression model
  #'@param : data - modeling data with dependent and independent variables 
  #'@param : max.lag - The highest number of lagged endogenous differenced variables to be included in the test regression
  #'@param : method - the method to use to test for co-integration
  #'@param : johansen.type - the test statistic type to use for Johansen type.
  #'@param : plot.it - logical to plot Long-Run- Equilibrium Relationship Among Model Variables

  #'@return : co-integration test output.
  
  if(!inherits(data, "data.frame"))stop("Data type must be a dataframe")
  if(!inherits(lm.obj, "lm"))stop("The class of 'lm.obj' must be 'lm' ")
  
  for(j in colnames(data)){
    test<-adf.test(data[,j])
    if(test$p.value < 0.05){cat(j, " series is stationary", "\n Series must be non-stationary or I(1) to use for the purposes of co-integration test\n")}
    else{cat("\n",j,"is non-stationary or I(1)\n")}
  }
  
  num.eigen= length(var.names) + 1     
  eg.column.names                     <-  c("Coint.Test.Type","Test.Type.On.Residuals", "Variant.of.Test.Type","Lag", "CritVal.5pct", "test.stats.val")
  jo.column.names                     <-  c("Coint.Test.Type", "Model.Type", "Model.Type.Desc", "Test.Stats.Type","Test.Stats.Type.Desc","Lag","Rank","Test.stats.val",  "CritVal.5pct") 
  eg.co.integration.out               <-  data.frame(matrix(NA, nrow=3, ncol= length(eg.column.names)))
  jo.co.integration.out               <-  data.frame(matrix(NA, nrow=3, ncol= length(jo.column.names)))
  jo.eigenval.out<-jo.eigenval.df     <-  jo.eigenval.out<-data.frame(matrix(NA, nrow=3, ncol= num.eigen))
  colnames(jo.eigenval.out)           <-  c(paste0("Eigenvalue", seq(num.eigen)))
  colnames(jo.eigenval.out)           <-  colnames(jo.eigenval.df)<- c(paste0("Eigenvalue", seq(num.eigen)))
  colnames(eg.co.integration.out)     <-  eg.column.names
  colnames(jo.co.integration.out)     <-  jo.column.names
  
  if(method == "Engle-Granger"){
    if(inherits(lm.obj, "lm")){
      resid.adf.none    <-ur.df(lm.residuals, type = "none", lags = max.lag, selectlags = "Fixed")
      resid.adf.drift   <-ur.df(lm.residuals, type = "drift", lags = max.lag,selectlags = "Fixed")
      resid.adf.trend   <-ur.df(lm.residuals, type = "trend", lags = max.lag,selectlags = "Fixed")
      eg.co.integration.out[1,]  <- list(Coint.Test.Type="Engle-Granger", Test.Type.On.Residuals= paste0(resid.adf.none@test.name, "-","None"), 
                                         Variant.of.Test.Type="tau1", Lag= max.lag, CritVal.5pct= resid.adf.none@cval[2], test.stats.val= resid.adf.none@teststat[1])
      eg.co.integration.out[2,]  <- list(Coint.Test.Type="Engle-Granger", Test.Type.On.Residuals= paste0(resid.adf.drift@test.name, "-", "Drift"), 
                                         Variant.of.Test.Type="tau2", Lag= max.lag, CritVal.5pct= resid.adf.drift@cval[1,2], test.stats.val= resid.adf.drift@teststat[1,1])
      eg.co.integration.out[3,]  <- list(Coint.Test.Type="Engle-Granger", Test.Type.On.Residuals= paste0(resid.adf.trend@test.name, "-", "Trend"), 
                                         Variant.of.Test.Type="tau3", Lag= max.lag, CritVal.5pct= resid.adf.trend@cval[1,2], test.stats.val= resid.adf.trend@teststat[1,1])
      return(eg.co.integration.out)
    }
  }else{
    if(johansen.type=="eigen"){
      johansen.test.none<-ca.jo(data, type ="eigen", ecdet = "none", K = max.lag, spec= "longrun")
      for(ii in 1:ncol(data)){
        jo.co.integration.out[ii,]<-list(Coint.Test.Type= "Johansen", Model.Type="none", Model.Type.Desc= johansen.test.none@model, Test.Stats.Type= "eigen", 
                                         Test.Stats.Type.Desc= johansen.test.none@type, Lag= max.lag, Rank= rownames(johansen.test.none@cval)[ii], 
                                         Test.stats.val= johansen.test.none@teststat[ii], CritVal.5pct= johansen.test.none@cval[ii,2])}
      for(ee in 1:num.eigen){jo.eigenval.df[1,ee] <- johansen.test.none@lambda[ee]}
      for(jj in c("const", "trend")){
        johansen.test<-ca.jo(data, type ="eigen", ecdet = jj, K = max.lag, spec= "longrun")
        for(ii in 1:ncol(data)){
          johansen.test.out <- list(Coint.Test.Type= "Johansen", Model.Type=jj , Model.Type.Desc= johansen.test@model, Test.Stats.Type= "eigen",
                                    Test.Stats.Type.Desc= johansen.test@type, Lag= max.lag, Rank= rownames(johansen.test@cval)[ii], 
                                    Test.stats.val= johansen.test@teststat[ii], CritVal.5pct= johansen.test@cval[ii,2])
          jo.co.integration.out  <- dplyr::bind_rows(jo.co.integration.out, johansen.test.out)} 
        for(ee in 1:num.eigen){jo.eigenval.out[1,ee] <- johansen.test@lambda[ee]}
        jo.eigenval.df   <- dplyr::bind_rows(jo.eigenval.df, jo.eigenval.out)
        jo.out           <- dplyr::bind_cols(jo.co.integration.out, jo.eigenval.df)
      }
    }else{
      johansen.test.none  <-ca.jo(data, type ="trace", ecdet = "none", K = max.lag, spec= "longrun")
      for(ii in 1:ncol(data)){
        jo.co.integration.out[ii,] <-list(Coint.Test.Type= "Johansen", Model.Type='none', Model.Type.Desc= johansen.test.none@model, 
                                          Test.Stats.Type= "trace", Test.Stats.Type.Desc= johansen.test.none@type, Lag= max.lag, 
                                          Rank= rownames(johansen.test.none@cval)[ii], Test.stats.val= johansen.test.none@teststat[ii], 
                                          CritVal.5pct= johansen.test.none@cval[ii,2])}
      for(ee in 1:num.eigen){jo.eigenval.df[1,ee] <- johansen.test.none@lambda[ee]}
      for(jj in c("const", "trend")){
        johansen.test <- ca.jo(data, type ="trace", ecdet = jj, K = max.lag, spec= "longrun")
        for(ii in 1:ncol(data)){
          johansen.test.out<-list(Coint.Test.Type= "Johansen", Model.Type=jj, Model.Type.Desc= johansen.test@model, Test.Stats.Type= "trace",
                                  Test.Stats.Type.Desc= johansen.test@type, Lag= max.lag, Rank= rownames(johansen.test@cval)[ii], 
                                  Test.stats.val= johansen.test@teststat[ii], CritVal.5pct= johansen.test@cval[ii,2])
          jo.co.integration.out <- dplyr::bind_rows(jo.co.integration.out, johansen.test.out)} 
        for(ee in 1:num.eigen){jo.eigenval.out[1,ee] <- johansen.test@lambda[ee]}
        jo.eigenval.df <- dplyr::bind_rows(jo.eigenval.df, jo.eigenval.out)
        jo.out <- dplyr::bind_cols(jo.co.integration.out, jo.eigenval.df)
      }
    } 
    
    if(plot.it){
      png(file= paste0("./output/","cointegration_plot.png"))
      p <- data %>% 
              gather(key="Series", value="value") %>%
              ggplot(aes(x=Date, y=value)) + geom_line(aes(color=Series), size=0.5) +
              ggtitle("Co-integration : Long-Run- Equilibrium Relationship Among Model Variables") + 
              ylab(label="Values") + 
              theme(plot.title=  element_text(size = 10, face = "bold"),
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    plot.background = element_rect(fill = "lightblue"),
                    panel.background = element_rect(fill = "light grey")
               ) + theme(legend.position = "bottom", legend.title=element_blank())
      dev.off()
    }
    print(p)
    return(jo.out)
  }
}

