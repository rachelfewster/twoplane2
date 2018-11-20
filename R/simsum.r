
dosim = function(alpha,Ec,k,w,sigmarate,planespd,D,En=100,Nsim=100,writeout=TRUE,seed=1,iomvt=FALSE,sigma.mult=5,L=NULL) {
  
  p.up = alpha # proportion of time up
  E1 = alpha*Ec
  p=c(1, 1) # definitely see if available in searched strip
  dmax.t=sigma.mult*(sigmarate*sqrt(k))/planespd # max time apart (in seconds) observations could be considered duplicates
  dmax.d = dmax.t*planespd # max distance apart observations could be considered duplicates
  ps = p.t(E1,Ec,p,sigmarate,k,dmax.t,planespd,w/2,io=iomvt) # capture history probabilities
  p. = sum(ps) # prob detect
  
  b = w
  if(iomvt) b = w + 2*dmax.d
  #  L=1100  # length in km from porpoise data
  #  L=En/(D*w*alpha) # set L to get desired sample size
  if(is.null(L)) L=En/(D*b*p.) # set L to get desired sample size
  
  p.see.up=c(1,1) # Prob see if up and in
  
  N=D*(L*b)
  
  Dstrip=N/(L*b) # density in number per sq km
  Dstrip.t=D*(planespd^2) # density in planespd units
  D.line.t=Dstrip.t*b/planespd # density in planespd along LINE units (1-dimensional)
  control.opt=list(trace=0,maxit=1000)
  
  estsio=ests=data.frame(Dhat=rep(0,Nsim),E1=rep(0,Nsim),E2=rep(0,Nsim),sigma=rep(0,Nsim),
                         n1=rep(0,Nsim),n2=rep(0,Nsim),m=rep(0,Nsim),mu_c=rep(0,Nsim),
                         se=rep(0,Nsim),inci=rep(0,Nsim))
  estsnspp=estsna=data.frame(Dhat=rep(0,Nsim),n1=rep(0,Nsim),n2=rep(0,Nsim))
  
  ests.kd=ests
  plot.sample=FALSE
  plot.displacement=FALSE
  plot.cuts=FALSE
  checkdists=FALSE
  fromfile=FALSE
  segiotime=segtime=nspptime=rep(NA,Nsim)
  
  estimate=c("D","sigma","E1") # parameters to estimate
  true=list(D=Dstrip,sigma=sigmarate,E1=E1) # parameters to use in simulation
  
  
  set.seed(seed) # initialise random number sequence (for repeatability)
  skip=c()
  startime=date()
  for(sim in 1:Nsim) {
    sdat=sim.2plane(N,L,w,sigmarate,k,planespd,p.up,Ec,p=p.see.up,
                    sigma.mult=sigma.mult,movement=list(forward=TRUE,sideways=iomvt))
    
    # fit accounting for leakage of animals in and out of strip
    segiotime[sim]=system.time(fitio<-segfit(sdat,D.line.t,E1=E1,Ec=Ec,sigmarate=sigmarate,planespd=planespd,p=c(1,1),sigma.mult=sigma.mult,
                                             control.opt=control.opt,method="BFGS",estimate=estimate,set.parscale=TRUE,
                                             io=iomvt,Dbound=NULL,hessian=TRUE))[3]
    
    # Palm
    pdat = format4palm(sdat,planespd,Ec,sigmarate,sigma.mult=sigma.mult)
    nspptime[sim]=system.time(palmfit<-fit.twocamera(pdat$points,pdat$cameras,pdat$d,pdat$w,pdat$b,pdat$l,pdat$tau,pdat$R,trace=FALSE))[3]
    est.palm=coef(palmfit)
    #  Convert 'activity centre' sigma of Palm into our sigma by *sqrt(2)  and convert to sigmarate.
    est.palm[3]=est.palm[3]*sqrt(2)/sqrt(sdat$k)
    #  
    estsio$n1[sim]=length(sdat$s1)
    estsio$n2[sim]=length(sdat$s2)
    estsnspp$n1[sim]=length(sdat$s1)
    estsnspp$n2[sim]=length(sdat$s2)
    ests$m[sim]=sdat$n11 # record number of actual duplicates
    estsio$m[sim]=sdat$n11 # record number of actual duplicates
    estsnspp$m[sim]=sdat$n11 # record number of actual duplicates
    
    infmat=try(solve(fitio$hessian),silent=TRUE)
    if(!inherits(infmat, "try-error")) {
      intest=logn.seci(log(fitio$D),sqrt(infmat[1,1]))
      estsio$se[sim]=intest$se/(b*planespd)
      estsio$inci[sim]=(intest$lower/(b*planespd)<=Dstrip & Dstrip<=intest$upper/(b*planespd))
    } else skip=c(skip,sim)
    estsio$Dhat[sim]=fitio$D/(b*planespd)
    estsio$E1[sim]=fitio$E[1]
    estsio$E2[sim]=fitio$E[2]
    estsio$sigma[sim]=fitio$sigmarate
    estsio$mu_c[sim]=fitio$mu_c
    
    estsnspp$Dhat[sim]=est.palm[1]
    estsnspp$sigma[sim]=est.palm[3]
    
    #    # Petersen estimator for known recaptures, within strip:
    #    ests.kd$Dhat[sim] = ((sdat$n1+1)*(sdat$n2+1)/(sdat$n11+1) - 1)/(L/planespd)
    
    #    if(sim==1) cat("\nCounter: \n")
    #    if(sim %% 50 == 0) cat(sim,"\n")
    #    else if(sim %% 10 == 0) cat(sim)
    #    else if(sim %% 5 == 0) cat("+")
    #    else cat("-")
  }  # End sim loop
  
  results = list(mle=estsio,palm=estsnspp)
  dir = "./inst/results/"
  fn = paste("sim-alpha_",alpha,"-k_",k,"-sigmarate_",signif(sigmarate,3),"-D_",D,"-En_",En,"-Nsim_",Nsim,".Rds",sep="")
  dirfn = paste(dir,fn,sep="")
  if(writeout) saveRDS(results,file=dirfn)
  else return(list(fn=fn,sim=results))
  
}  


harvestsim = function(alpha,k,sigmarate,D,En,Nsim,badcut=10,simresults=NULL) {
  
  if(is.null(simresults)) {
    fn = paste("./inst/results/sim-alpha_",alpha,"-k_",k,"-sigmarate_",signif(sigmarate,3),"-D_",D,"-En_",En,"-Nsim_",Nsim,".Rds",sep="")
    sim = readRDS(fn)
  } else {
    fn=simresults$fn
    sim=simresults$sim
  }
  
  alpha = as.numeric(strsplit(strsplit(fn,"alpha_")[[1]][2],"-")[[1]][1])
  k = as.numeric(strsplit(strsplit(fn,"k_")[[1]][2],"-")[[1]][1])
  sigmarate = as.numeric(strsplit(strsplit(fn,"sigmarate_")[[1]][2],"-")[[1]][1])
  D = as.numeric(strsplit(strsplit(fn,"D_")[[1]][2],"-")[[1]][1])
  En = as.numeric(strsplit(strsplit(fn,"En_")[[1]][2],"-")[[1]][1])
  Nsim = as.numeric(strsplit(strsplit(fn,"Nsim_")[[1]][2],".Rds")[[1]][1])
  
  badD = abs(sim[[1]]$Dhat-D)>=badcut | abs(sim[[2]]$Dhat)>=badcut 
  nbadD = sum(badD)
  bad = abs(sim[[1]]$Dhat-D)>=badcut | abs(sim[[2]]$Dhat)>=badcut | is.na(sim[[1]]$se)
  nbad = sum(bad)
  nbadse = nbad - nbadD
  
  pc.bias.mle = 100*(mean(sim[[1]]$Dhat[!bad])-D)/D
  pc.cv.mle = 100*sqrt(var(sim[[1]]$Dhat[!bad]))/mean(sim[[1]]$Dhat[!bad])
  cover.mle = sum(sim[[1]]$inci[!bad])/Nsim
  
  pc.bias.palm = 100*(mean(sim[[2]]$Dhat[!bad])-D)/D
  pc.cv.palm = 100*sqrt(var(sim[[2]]$Dhat[!bad]))/mean(sim[[1]]$Dhat[!bad])
  
  n1 = mean(sim[[1]]$n1[!bad])
  n2 = mean(sim[[1]]$n2[!bad])
  m = mean(sim[[1]]$m[!bad])
  
  E1 = mean(sim[[1]]$E1[!bad])
  
  sigmrate = mean(sim[[1]]$sigma[!bad])
  
  sehat.Dhat = mean(sim[[1]]$se[!bad])
  
  out = data.frame(Nsim=Nsim,alpha=alpha,k=k,speed=getspeed(sigmarate),D=D,
                   pc.bias.mle=pc.bias.mle,pc.cv.mle=pc.cv.mle,cover.mle=cover.mle,
                   pc.bias.palm=pc.bias.palm,pc.cv.palm=pc.cv.palm,
                   n1=n1,n2=n2,m=m,
                   E1=E1,
                   sigmarate=sigmarate,
                   sehat.Dhat=sehat.Dhat,
                   nbadD=nbadD,nbadse=nbadse)
  
  return(out)
}  

getspeed = function(sigmarate) return(1000*sigmarate*(sqrt(2)*gamma(1)/gamma(0.5)))

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
