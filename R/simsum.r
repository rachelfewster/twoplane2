
dosim = function(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=c(1,1),movement=list(forward=TRUE,sideways=TRUE),
                 fix.N=TRUE,En=NULL,Nsim=100,writeout=TRUE,seed=1,simethod="MLE",control.opt=control.opt,hessian=TRUE) {

  dmax.t = (b-w)/planespd
  ps = p.t(kappa,tau,p,sigmarate,k,dmax.t,planespd,w,io=movement$sideways) # capture history probabilities
#  p. = sum(ps) # prob at least one observer detects
  p1 = mean((ps$ch01+ps$ch11),(ps$ch10+ps$ch11)) # mean prob single observer detects
  
  N=D.2D*(L*2*b)
  
  Ltype = "FixedL"
  Ntype = "FixedN"
  if(!fix.N) Ntype = "RandomN"
  if(!is.null(En)) {
    if(!is.null(L)) warning("Ignoring L because En was specified; L has been calculated to give En with given D.2D.")
    L=En/(D.2D*b*p.) # set L to get desired sample size
    Ltype = "RandomL"
  }
  if(is.null(En)) {
    if(is.null(L)) stop("Need to specify one of L and En.")
    En = round(N*p1)
  } 
  
  mlests=data.frame(D.est=rep(0,Nsim),D.se=rep(0,Nsim),D.cv=rep(0,Nsim),D.inci=rep(0,Nsim),
                    gamma.est=rep(0,Nsim),gamma.se=rep(0,Nsim),gamma.cv=rep(0,Nsim),gamma.inci=rep(0,Nsim),
                    sigmarate.est=rep(0,Nsim),sigmarate.se=rep(0,Nsim),sigmarate.cv=rep(0,Nsim),sigmarate.inci=rep(0,Nsim),
                    n1=rep(0,Nsim),n2=rep(0,Nsim),m=rep(0,Nsim),tau=rep(0,Nsim))
  palmests=data.frame(Dhat=rep(0,Nsim),kappa=rep(0,Nsim),sigmarate=rep(0,Nsim),n1=rep(0,Nsim),n2=rep(0,Nsim))
  
#  ests.kd=ests
#  plot.sample=FALSE
#  plot.displacement=FALSE
#  plot.cuts=FALSE
#  checkdists=FALSE
#  fromfile=FALSE
  segiotime=segtime=nspptime=rep(NA,Nsim)
  
  estimate=c("D","sigma","E1") # parameters to estimate
#  true=list(D=D.2D,sigma=sigmarate,E1=kappa) # parameters to use in simulation
  
# Convert from sigmarate to Ben's sigma
  sigma = sigmarate/(sqrt(2)/sqrt(k))
  
  set.seed(seed) # initialise random number sequence (for repeatability)
#  skip=c()
  startime=date()
  for(sim in 1:Nsim) {
    if(simethod=="MLE") {
      sdat = sim.2plane(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=c(1,1),movement=movement,fix.N=fix.N)
    } else if(simethod=="Palm") {
      Ntype = "RandomN"
      palmdat <- sim.twocamera(c(D.2D=D.2D, kappa=kappa, sigma=sigma),d=L, w=w, b=b, l=k, tau=tau)
      sdat = Palm2mleData(palmdat$points,palmdat$sibling.list$cameras,d=L,l=k,w=w,b=b)
    } else stop("simethod must be 'MLE' or 'Palm'.")
    
    # MLE
    mlefit<-segfit(sdat,D.2D,E1=kappa,Ec=tau,sigmarate=sigmarate,planespd=planespd,p=c(1,1),
                   sigma.mult=sigma.mult,control.opt=control.opt,method="BFGS",estimate=estimate,
                   set.parscale=TRUE,io=TRUE,Dbound=NULL,hessian=hessian)
    if(hessian) {
      # Density
      mlests$D.se[sim]=mlefit$D["est"]
      mlests$D.cv[sim]=mlefit$D["cv"]
      mlests$D.inci[sim]=(mlefit$D["lcl"]<=D.2D & D.2D<=mlefit$D["ucl"])
      # gamma
      gamma = kappa/tau
      mlests$gamma.se[sim]=mlefit$gamma["est"]
      mlests$gamma.cv[sim]=mlefit$gamma["cv"]
      mlests$gamma.inci[sim]=(mlefit$gamma["lcl"]<=gamma & gamma<=mlefit$gamma["ucl"])
      # sigmarate
      mlests$sigmarate.se[sim]=mlefit$sigmarate["est"]
      mlests$sigmarate.cv[sim]=mlefit$sigmarate["cv"]
      mlests$sigmarate.inci[sim]=(mlefit$sigmarate["lcl"]<=sigmarate & sigmarate<=mlefit$sigmarate["ucl"])
    } #else skip=c(skip,sim)
    mlests$D.est[sim]=mlefit$D["est"]
    mlests$gamma.est[sim]=mlefit$gamma["est"]
    mlests$sigmarate.est[sim]=mlefit$sigmarate["est"]
    mlests$tau[sim]=mlefit["tau"]

    # Palm
    palmfit<-twoplane.fit(sdat,tau=tau,R=550,all=TRUE)
    palmests$Dhat[sim]=coef(palmfit)["D.2D"]
    palmests$kappa[sim]=coef(palmfit)["kappa"]
    palmests$sigmarate[sim]=coef(palmfit)["sigma"]*sqrt(2)/sqrt(sdat$k) # Convert to Brownian sigmarate
    
    # Sample size statistics
    mlests$n1[sim]=length(sdat$y1)
    mlests$n2[sim]=length(sdat$y2)
    palmests$n1[sim]=length(sdat$y1)
    palmests$n2[sim]=length(sdat$y2)
    if(simethod=="MLE") mlests$m[sim]=sdat$n11 # record number of actual duplicates

  }  # End sim loop
  
  results = list(mle=mlests,palm=palmests)
  dir = "./inst/results/"
  fn = paste("sim-gamma_",signif(gamma,3),"-tau_",signif(tau,3),"-k_",k,"-sigmarate_",signif(sigmarate,3),"-D_",signif(D.2D,3),
             "-En_",signif(En,3),"-",Ltype,"-",Ntype,"-simethod_",simethod,"-Nsim_",Nsim,".Rds",sep="")
  dirfn = paste(dir,fn,sep="")
  if(writeout) {
    saveRDS(results,file=dirfn)
    invisible(list(file=dirfn))
  } else return(list(fn=fn,sim=results))
  
}  


harvestsim = function(fn,badcut=100) {
  
  sim = readRDS(fn)

  gamma = as.numeric(strsplit(strsplit(fn,"gamma_")[[1]][2],"-")[[1]][1])
  tau = as.numeric(strsplit(strsplit(fn,"tau_")[[1]][2],"-")[[1]][1])
  k = as.numeric(strsplit(strsplit(fn,"k_")[[1]][2],"-")[[1]][1])
  sigmarate = as.numeric(strsplit(strsplit(fn,"sigmarate_")[[1]][2],"-")[[1]][1])
  D.2D = as.numeric(strsplit(strsplit(fn,"D_")[[1]][2],"-")[[1]][1])
  En = as.numeric(strsplit(strsplit(fn,"En_")[[1]][2],"-")[[1]][1])
  Nsim = as.numeric(strsplit(strsplit(fn,"Nsim_")[[1]][2],".Rds")[[1]][1])
  
  badD = abs(sim[[1]]$D.est)>=badcut*D.2D | abs(sim[[2]]$Dhat)>=badcut*D.2D
  nbadD = sum(badD)
  bad = abs(sim[[1]]$D.est)>=badcut*D.2D | abs(sim[[2]]$Dhat)>=badcut*D.2D | is.na(sim[[1]]$D.se)
  nbad = sum(bad)
  nbadse = nbad - nbadD
  
  pc.bias.mle = 100*(mean(sim[[1]]$D.est[!bad])-D.2D)/D.2D
  pc.cv.mle = 100*sqrt(var(sim[[1]]$D.est[!bad]))/mean(sim[[1]]$D.est[!bad])
  cover.mle = sum(sim[[1]]$D.inci[!bad])/(Nsim-nbad)
  
  pc.bias.palm = 100*(mean(sim[[2]]$Dhat[!bad])-D.2D)/D.2D
  pc.cv.palm = 100*sqrt(var(sim[[2]]$Dhat[!bad]))/mean(sim[[2]]$Dhat[!bad])
  
  Dhat.cor = cor(sim[[1]]$D.est[!bad],sim[[2]]$Dhat[!bad])
  
  n1 = mean(sim[[1]]$n1[!bad])
  n2 = mean(sim[[1]]$n2[!bad])
  m = mean(sim[[1]]$m[!bad])
  
  gamma = mean(sim[[1]]$gamma.est[!bad])
  
  sigmrate = mean(sim[[1]]$sigmarate.est[!bad])
  
  sehat.Dhat = mean(sim[[1]]$D.se[!bad])
  
  out = data.frame(Nsim=Nsim,gamma=gamma,k=k,speed=getspeed(sigmarate),D.2D=D.2D,
                   pc.bias.mle=pc.bias.mle,pc.cv.mle=pc.cv.mle,cover.mle=cover.mle,
                   pc.bias.palm=pc.bias.palm,pc.cv.palm=pc.cv.palm,
                   Dhat.cor=Dhat.cor,
                   n1=n1,n2=n2,m=m,
                   gamma=gamma,
                   sigmarate=sigmarate,
                   sehat.Dhat=sehat.Dhat,
                   nbadD=nbadD,nbadse=nbadse)
  
  return(out)
}  

getspeed = function(sigmarate) return(1000*sigmarate*(sqrt(2)*gamma(1)/gamma(0.5)))



boxplotsim = function(fns,method="mle",stat="D.est",diff.method=NULL,diff.from=NULL,sortby=" ",...) {
  
  if(method=="MLE") method = "mle"
  if(method=="Palm") method = "palm"
  if(!is.element(method,c("mle","palm"))) stop("method must be `mle` or `palm`")
  if(!is.null(diff.method)) if(!is.element(diff.method,c("mle","palm"))) stop("diff.method must be `mle` or `palm`")
  diff.statname = statname = stat
  if(method=="palm" & stat=="D.est") statname = "Dhat" # palm name is different
  if(!is.null(diff.method)) if(diff.method=="palm" & stat=="D.est") diff.statname = "Dhat" # palm name is different
  if(sortby!=" ") if(sortby!="n" & sortby!="m") stop("sortby must b n or m")
  
  nscenarios = length(fns)
  sim1 = readRDS(fns[[1]])
  nsims = dim(sim1$mle)[1]
  plotstat = matrix(rep(NA,nsims*nscenarios),nrow=nsims)
  n = m = rep(NA,nscenarios)
  
  for(ns in 1:nscenarios) {
    sim = readRDS(fns[[ns]])
    if(method=="mle" & !is.element(stat,names(sim$mle))) stop("Invalid stat for method mle")
    if(method=="palm" & !is.element(statname,names(sim$palm))) stop("Invalid stat for method palm")
    if(is.null(diff.method)) {
      plotstat[,ns] = sim[[method]][,statname]
    } else if(is.numeric(diff.from)) {
      plotstat[,ns] = sim[[method]][,statname] - diff.from
    } else {
      plotstat[,ns] = sim[[method]][,statname] - sim[[diff.method]][,diff.statname]
    }
    n[ns] = (mean(sim$mle$n1) + mean(sim$mle$n2))/2
    m[ns] = mean(sim$mle$m)
  }
  ord = simnum = 1:nscenarios
  if(sortby=="n") ord = order(n)
  if(sortby=="m") ord = order(m)
  boxplot(plotstat[,ord],names=as.character(ord),xlab="Scenario",...)
}
  
#quartz(h=8,w=10)
#boxplotsim (fns,method="palm",stat="D.est",ylab=expression(hat(D)))
#boxplotsim (fns,method="mle",stat="D.est",diff.method="palm")
#lines(c(0,37),rep(0,2))


calcsimbias = function(fns,method="mle",stat="D.est",truth) {
  
  if(method=="MLE") method = "mle"
  if(method=="Palm") method = "palm"
  statname = stat
  if(!is.element(method,c("mle","palm"))) stop("method must be `mle` or `palm`")
  if(method=="palm" & stat=="D.est") statname = "Dhat" # palm name is different

  nscenarios = length(fns)
  sim1 = readRDS(fns[[1]])
  nsims = dim(sim1$mle)[1]
  nas = rep(NA,nscenarios)
  eststats = data.frame(mean=nas,se.mean=nas,lcl=nas,ucl=nas,biased=nas)

  for(ns in 1:nscenarios) {
    sim = readRDS(fns[[ns]])
    if(method=="mle" & !is.element(stat,names(sim$mle))) stop("Invalid stat for method mle")
    if(method=="palm" & !is.element(statname,names(sim$palm))) stop("Invalid stat for method palm")
    eststats$mean[ns] = mean(sim[[method]][,statname])
    eststats$se.mean[ns] = sqrt(var(sim[[method]][,statname])/nsims)
    eststats$lcl[ns] = eststats$mean[ns]-1.96*eststats$se.mean[ns]
    eststats$ucl[ns] = eststats$mean[ns]+1.96*eststats$se.mean[ns]
    eststats$biased[ns] = truth<eststats$lcl[ns] | eststats$ucl[ns]<truth
  }
  return(eststats)
}

#biassum = calcsimbias (fns,method="mle",stat="D.est",truth=1.24)
#sum(biassum$biased)/36

#biassum = calcsimbias (fns,method="palm",stat="D.est",truth=1.24)
#sum(biassum$biased)/36



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
