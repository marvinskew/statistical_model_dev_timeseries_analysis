
pkages  <-c("urca","tseries","data.table","dplyr","broom","lubridate","car","lmtest","sandwich","orcutt", "ggplot2","zoo", "tidyr", "astsa")
for(pkage in pkages){
  if(!pkage %in% installed.packages())install.packages(pkage)
  
}
lapply(pkages, require, character.only= TRUE)

serial.correlationTest<-function(lm.obj, max.lag, plot.it=TRUE){
  stopifnot(class(lm.obj)=="lm")
  
  serial.corr.test.results <- rbindlist(lapply(1:max.lag, function(s, lm.obj, order, type){ 
    bg.tests             <- bgtest(lm.obj, order= s,  type = "Chisq")
    lj.box.tests         <- Box.test(resid(lm.obj), lag= s, type ="Ljung-Box")
    dw.tests             <- dwtest(model.obj, alternative = c("two.sided"))
    bg.test.result       <- list(OLS.Assumption= "No Serial Correlation", AutoReg.Functional.Form= paste0("AR(",s,")"),Test.Name="Breusch-Godfrey Test",
                                 Lag=s, Test.Statistic.Value= bg.tests$statistic ,Degres.of.Freedom= bg.tests$parameter, P.Value= bg.tests$p.value)
    lj.box.test.result   <- list(OLS.Assumption= "No Serial Correlation", AutoReg.Functional.Form= paste0("AR(",s,")"),Test.Name="Box-Ljung Test",
                                 Lag=s, Test.Statistic.Value= lj.box.tests$statistic ,Degres.of.Freedom= lj.box.tests$parameter, P.Value= lj.box.tests$p.value)
    dw.test.result       <- list(OLS.Assumption= "No Serial Correlation", AutoReg.Functional.Form= "AR(1)",Test.Name= "Dublin-Watson Test",
                                 Lag= 1, Test.Statistic.Value= dw.tests$statistic[[1]] ,Degres.of.Freedom= NA, P.Value= dw.tests$p.value)
    return(do.call(rbindlist, list(dw.test.result, bg.test.result, lj.box.test.result), use.names=TRUE, fill=TRUE))
  }, lm.obj = lm.obj, order = order, type =type), use.names=TRUE, fill=TRUE)
  
  lm.residuals <- 
  
  if(plot.it){
    p<-residual.df %>% 
      rename(residuals= .resid) %>%
      ggplot(aes(x=Date, y= residuals)) + geom_line(size=0.5) + geom_point(col="blue") +
      ggtitle("Serial Correlation - Residual Plot") + 
      ylab(label="Residuals") + 
      theme(  plot.title=  element_text(hjust= 0.5, size = 8, face = "bold"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.background = element_rect(fill = "lightblue"),
              panel.background = element_rect(fill = "light grey")
      ) + theme(legend.position = "bottom")
    
  }
  
  return(serial.corr.test.results)
}