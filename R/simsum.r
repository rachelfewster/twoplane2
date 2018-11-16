#---------------------------- simulation summary functions ---------------------
simsum=function(ests,true=NULL,simtimes=NULL) {
  if(is.null(true)) stop("Need to pass true values")
  D=true$D
  sigma=true$sigma
  E1=true$E1
  
  na.se=which(is.na(ests$se) | is.nan(ests$se))
  if(length(na.se)>0) {
    warning(paste(length(na.se)," se estimates were NA or NaN; these have not been included.",sep=""))
    ests=ests[-na.se,]
  }
  huge.se=which(ests$se > true$D*100)
  if(length(huge.se)>0) {
    warning(paste(length(huge.se)," se estimates were > 100 times bigger than true D; these have not been included.",sep=""))
    ests=ests[-huge.se,]
  }
  
  corrmat=cor(ests[c(1,4,2)])
  covmat=cov(ests[c(1,4,2)])
  se=sqrt(covmat[1,1])
  cv.n1=round(100*sd(ests$n1)/mean(ests$n1),1)
  cv.n2=round(100*sd(ests$n2)/mean(ests$n2),1)
  pc.cvhat.Dhat=100*ests$se/ests$Dhat
  pc.cv.pc.cvhat.Dhat=round(100*sd(pc.cvhat.Dhat)/mean(pc.cvhat.Dhat),1)
  pc.cv.Dhat=round(100*se/D,1)
  
  cat("------------------\n")
  cat(" Nsim = ",Nsim,"\n")
  cat("------------------\n")
  if(!is.null(simtimes)) cat("Mean time for fit : ",signif(mean(simtimes),5),"\n\n")

  cat(paste("Obs 1 sample size: ",round(mean(ests$n1),1),"  (%CV=",cv.n1,")",sep=""),"\n")
  cat(paste("Obs 2 sample size: ",round(mean(ests$n2),1),"  (%CV=",cv.n2,")",sep=""),"\n\n")

  simbias.Dhat=round(mean((ests$Dhat-D)/D)*100,2)
  simbias.sehat=round(mean((ests$se-se)/se)*100,2)
  simse=apply(ests[c(1,4,2)],2,sd)[1]
  simsebias.Dhat=round(mean((ests$se-simse)/simse)*100,2)
  cat(paste("Bias in Dhat       : ",simbias.Dhat,"%",sep=""),"\n")
  cat(paste("Bias in sigma      : ",round(mean((ests$sigma-sigma)/sigma)*100,2),"%",sep=""),"\n")
  cat(paste("Coverage prob      : ",round(mean(ests$inci)*100,1),"%",sep=""),"\n")
  cat(paste("Bias in E(up)      : ",round(mean((ests$E1-E1)/E1)*100,2),"%",sep=""),"\n")
  cat(paste("Bias in sehat(Dhat): ",simbias.sehat,"%",sep=""),"\n\n")
  cat(paste("%CV(Dhat): ",pc.cv.Dhat,"%",sep=""),"\n")
  cat(paste("%CVhat(Dhat): ",round(mean(pc.cvhat.Dhat),1),"  (%CV=",pc.cv.pc.cvhat.Dhat,")",sep=""),"\n")
  
  cat("\n          Correlation matrix:\n","         -------------------\n")
  print(corrmat)
  cat("\n          Covariance matrix:\n","        -------------------\n")
  print(covmat)
}



simsum.nspp=function(ests,true=NULL,simtimes=NULL) {
  if(is.null(true)) stop("Need to pass true values")
  D=true$D;sigma=true$sigma;E1=true$E1
  cv.n1=round(100*sd(ests$n1)/mean(ests$n1),1)
  cv.n2=round(100*sd(ests$n2)/mean(ests$n2),1)

  corrmat=cor(ests[c("Dhat","sigma")])
  covmat=cov(ests[,c("Dhat","sigma")])
  se=sqrt(covmat[1,1])
  pc.cv.Dhat=round(100*se/D,1)

  cat("------------------\n")
  cat(" Nsim = ",Nsim,"\n")
  cat("------------------\n")
  if(!is.null(simtimes)) cat("Mean time for fit : ",signif(mean(simtimes),5),"\n\n")
  
  cat(paste("Obs 1 sample size: ",round(mean(ests$n1),1),"  (%CV=",cv.n1,")",sep=""),"\n")
  cat(paste("Obs 2 sample size: ",round(mean(ests$n2),1),"  (%CV=",cv.n2,")",sep=""),"\n\n")
  
  simbias.Dhat=round(mean((ests$Dhat-Dstrip)/Dstrip)*100,2)

  cat(paste("Bias in Dhat: ",simbias.Dhat,"%",sep=""),"\n")
  cat(paste("Bias in sigma: ",round(mean((ests$sigma-sigma)/sigma)*100,2),"%",sep=""),"\n")
  cat(paste("%CV(Dhat): ",pc.cv.Dhat,"%",sep=""),"\n\n")

  cat("\n          Correlation matrix:\n","         -------------------\n")
  print(cor(ests[c("Dhat","sigma")]))
  cat("\n          Covariance matrix:\n","        -------------------\n")
  print(cov(ests[c("Dhat","sigma")]))
}
#---------------------------- simulation summary functions ---------------------
