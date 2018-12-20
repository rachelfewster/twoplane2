sim.2plane=function(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=c(1,1),
                    movement=list(forward=TRUE,sideways=TRUE),fix.N=FALSE,sim.ft.normal=TRUE){
  
  if(!movement$sideways) b = w
  E.N = D.2D * 2*b*L
  if(fix.N) {
    NumSimAnimals = round(E.N) # Set N to expected number
  } else {
    NumSimAnimals = rpois(1,E.N) # Draw N from Poisson
  }
  p1=p[1]
  p2=p[2]
  
#  b=sigma.mult*sigmarate*sqrt(k) # max dist apart (in km) observations could be considered duplicates; 
  # Also width of buffer around strip.
  # Brownian mvt has sd of sigmarate*sqrt(time passed)
  tb=b/planespd   # half-width of buffered stip, in plane seconds

  tL=L/planespd # transect length, in plane seconds
  tw=w/planespd # half-width of searched stip, in plane seconds

  
  l=sort(runif(NumSimAnimals,0,tL)) # generate animal locations along line (in plane seconds)
  lhoriz=runif(NumSimAnimals,-tb,tb)  #  Place animals uniformly in buffered strip.
  if(movement$forward) {
    if(sim.ft.normal) simt = rnorm(NumSimAnimals,mean=0,sd=sqrt(k)*sigmarate/planespd)
    else simt = rft(NumSimAnimals,k,planespd,log(sigmarate)) - k # deviation from lag of k
  } else { # no forward movement
    simt=rep(0,NumSimAnimals)
  }
  l2=l+simt # locations of animals (in plane time since start of transect, time measured in plane seconds) when 2nd observer passes
  l2horiz=lhoriz  #  Horizontal position when second plane passes. 
  if(movement$sideways) {
    simthoriz=rnorm(NumSimAnimals, 0, sqrt(k+simt)*sigmarate/planespd) # horizontal movement (in plane seconds)
    l2horiz=l2horiz+simthoriz
  }
  # simulate detection by leading observer:
  up.1=1*(runif(NumSimAnimals)<=(kappa/tau)) # random variable with 1=available, 0=unavailable; kappa/tau is stationary dbn prob of being up
  in.1=(lhoriz<=tw)*(lhoriz>=-tw)   #  Is it in the strip?
  see.1=(runif(NumSimAnimals)<=p[1])*up.1*in.1 # binary variable indicating whether or not plane 1 sees
  
  in.2=(l2horiz<=tw)*(l2horiz>=-tw)   #  Is it in the strip.
  p.01.2=p.11.2=rep(0,NumSimAnimals)
  t12=k+simt    # seconds elapsing between 1st and 2nd observer passing animals
  q11=-1/kappa
  q22=-1/(tau-kappa)
  Qmat=matrix(c(q11,-q22,-q11,q22),nrow=2)
  for(i in 1:NumSimAnimals) {
    p.01.2[i]=p.omega.t(t12[i], idbn=c(up.1[i],(1-up.1[i])), p1, p2, Qmat, omega=01)  # Prob(01 | in and up/down state when 1st passes)
    p.11.2[i]=p.omega.t(t12[i], idbn=c(up.1[i],(1-up.1[i])), p1, p2, Qmat, omega=11)  # Prob(11 | in and up/down state when 1st passes)
  }
  p.see.2.given.in = (p.01.2+p.11.2) # prob 2nd sees, given in for first
  see.2=(runif(NumSimAnimals)<=p.see.2.given.in*in.2) # binary variable indicating whether or not plane 2 sees
  n1=sum(see.1)
  n2=sum(see.2)
  dups=(see.1*see.2==1)
  m=sum(dups)
  
  sno=1:length(see.1)
  dno=which(dups)
  sno1=sno[see.1==1]
  sno2=sno[see.2==1]
  
  s1=l[sno1]
  s2=l2[sno2]
  sh1=lhoriz[sno1]
  sh2=l2horiz[sno2]
  # convert from times to distances
  x1 = sh1*planespd
  x2 = sh2*planespd
  y1 = s1*planespd
  y2 = s2*planespd
  # convert from times to distances
  ally1 = l*planespd
  ally2 = l2*planespd
  allx1 = lhoriz*planespd
  allx2 = l2horiz*planespd
  # mark duplicates:
  d1=which(is.element(sno1,dno))
  d2=which(is.element(sno2,dno))
  n1=length(s1)
  n2=length(s2)
  
  sim.t.area=tL*2*tb
  sim.t.D=NumSimAnimals/sim.t.area
  sim.area=sim.t.area*planespd
  sim.D=NumSimAnimals/sim.area
  
  simdat=list(s1=s1,s2=s2,sh1=sh1,sh2=sh2,x1=x1,x2=x2,y1=y1,y2=y2,tL=tL,tw=tw,tb=tb,k=k,L=L,w=w,b=b,
              sno1=sno1,sno2=sno2,d1=d1,d2=d2,n1=n1,n2=n2,n11=m,n=n1+n2-m,
              alls1=l,allsh1=lhoriz,alls2=l2,allsh2=l2horiz,ally1=ally1,ally2=ally2,allx1=allx1,allx2=allx2,
              dups=dups,sim.t.area=sim.t.area,sim.t.D=sim.t.D,sim.area=sim.area,sim.D=sim.D,
              p.01.2=p.01.2,p.11.2=p.11.2
              ,see.1=see.1, see.2=see.2,
              dy=simt*planespd,dx=simthoriz*planespd)
  class(simdat)=c("sim2plane",class(simdat))
  return(simdat)
}



plot.sim2plane=function(simdat,all=TRUE,recaps=TRUE,joindups=recaps,w=NULL,L=NULL,xlim=NULL){
  if(!is.element("sim2plane",class(simdat))) stop("simdat must be of class 'sim2plane'.")
  if(is.null(xlim)) xlim=range(simdat$alls1,simdat$alls2)
    if(!all) {
    plot(c(simdat$s1,simdat$s2),c(simdat$x1,simdat$x2),pch="+",col=c(rep("blue",simdat$n1),rep("red",simdat$n2)),
         xlab="Observer seconds",ylab="Observer seconds",xlim=xlim)
  } else {
    plot(simdat$alls1,simdat$allx1,pch="+",col="gray",xlab="Observer seconds",ylab="Observer seconds"
         ,xlim=xlim)
    points(c(simdat$s1,simdat$s2),c(simdat$x1,simdat$x2),pch="+",col=c(rep("blue",simdat$n1),rep("red",simdat$n2)))
  }
  if(recaps) {
    points(simdat$s1[simdat$d1],simdat$x1[simdat$d1],col="blue")
    points(simdat$s2[simdat$d2],simdat$x2[simdat$d2],col="red")
    segments(simdat$alls1[simdat$dups],simdat$allx1[simdat$dups],simdat$alls2[simdat$dups],simdat$allx2[simdat$dups],col="gray")
  }
  if(is.null(L)) smax=max(simdat$s1,simdat$s2)
  else smax=L
  if(!is.null(w)) segments(rep(0,2),c(-w,w),rep(L,2),c(-w,w),col="gray")
}



plotsimloc=function(l,l2,dups,tm=range(l,l2),newindow=FALSE){  
  #----------------------------------------------------------------------------------------
  # Plots locations of simulated detections in l and l2, linking duplicates via dups; 
  # tm is time range.
  #----------------------------------------------------------------------------------------
  if(newindow) quartz(h=3)
  t1=l[tm[1]<=l & l<=tm[length(tm)] & see.1==1]
  t2=l2[tm[1]<=l2 & l2<=tm[length(tm)]& see.2==1]
  d1=l[dups]
  d2=l2[dups]
  plot(t1,rep(1,length(t1)),xlim=range(t1,t2),ylim=c(1,2),xlab="Time along transect",ylab="Observer",yaxp=c(1,2,1))
  points(t2,rep(2,length(t2)),col="blue")
  nds=length(d1)
  segments(d1,rep(1,nds),d2,rep(2,nds),lty=2,col="red")
  points(d1,rep(1,nds),pch=19,col="black")
  points(d2,rep(2,nds),pch=19,col="blue")
}



plot.ests=function(ests,D,kappa,sigma, mu_c, sigma.est=TRUE,nclass=20,keep=NULL)
#-------------------------------------------------------------------------------
# NB: Not sure this is an up-to-date function.
# Plots histograms of simulated estimates and parameters and prints bias, CV
# and parameter correlations.
#
# Inputs: 
# ------
#  ests: dataframe with columns 
#        Dhat : Density estimate 
#        kappa   : Estimated expected time available
#        E2   : Estimated expected time UNavailable
#        sigma: std. err. of (normal) animal movement (dispersion) model
#        n1   : number of observer 1 detections
#        n2   : number of observer 2 detections
#        m    : number of duplicate detections
#  sigma.est: TRUE if sigma was estimated (if not, does not do sigma plot)
#  nclass: number of histogram bars in plots
#-------------------------------------------------------------------------------
{
  # windows()
  par(mfrow=c(2,3))
  
  nest=length(ests$Dhat)
  
  if(is.null(keep)) keep=1:length(ests$Dhat)
  counts=hist(ests$Dhat[keep],plot=FALSE)$counts
  hist(ests$Dhat[keep],main="(a)",xlab=expression(hat(D)),xlim=range(D,ests$Dhat[keep]),nclass=nclass)
  lines(rep(mean(ests$Dhat[keep]),2),range(c(0,counts)),col="blue",lwd=2)
  lines(rep(D,2),range(0,counts),lty=2,col="red",lwd=2)
  E.Dhat=mean(ests$Dhat[keep])
  se.Dhat=sqrt(var(ests$Dhat[keep]))
  RMSE.Dhat=sqrt(mean((ests$Dhat[keep]-D)^2))
  cv.Dhat=se.Dhat/E.Dhat
  pcbias.Dhat=100*(E.Dhat-D)/D
  pcbias.Dhat.CI=pcbias.Dhat+100*c(-2,2)*se.Dhat/sqrt(nest)
  
  counts=hist(ests$kappa[keep],plot=FALSE)$counts
  hist(ests$kappa[keep],main="(b)",xlab=expression(paste(hat(E),"[up]",sep="")),xlim=range(kappa,ests$kappa[keep]),nclass=nclass)
  lines(rep(mean(ests$kappa[keep]),2),range(c(0,counts)),col="blue",lwd=2)
  lines(rep(kappa,2),range(c(0,counts)),lty=2,col="red")
  E.kappa=mean(ests$kappa[keep])
  se.kappa=sqrt(var(ests$kappa[keep]))
  RMSE.kappa=sqrt(mean((ests$kappa[keep]-kappa)^2))
  cv.kappa=se.kappa/E.kappa
  pcbias.kappa=100*(E.kappa-kappa)/kappa
  
  pcbias.sigma=NA
  cv.sigma=NA
  if(sigma.est) {
    counts=hist(ests$sigma[keep],plot=FALSE)$counts
    hist(ests$sigma[keep],main="(c)",xlab=expression(hat(sigma)),xlim=range(sigma,ests$sigma[keep]),nclass=nclass)
    lines(rep(mean(ests$sigma[keep]),2),range(c(0,counts)),col="blue",lwd=2)
    lines(rep(sigma,2),range(c(0,counts)),lty=2,col="red")
    E.sigma=mean(ests$sigma[keep])
    se.sigma=sqrt(var(ests$sigma[keep]))
    RMSE.sigma=sqrt(mean((ests$sigma[keep]-sigma)^2))
    cv.sigma=se.sigma/E.sigma
    pcbias.sigma=100*(E.sigma-sigma)/sigma
  }
  
  counts=hist(ests$mu_c[keep],plot=FALSE)$counts
  hist(ests$mu_c[keep],main="(c)",xlab=expression(hat(mu_c)),xlim=range(mu_c,ests$mu_c[keep]),nclass=nclass)
  lines(rep(mean(ests$mu_c[keep]),2),range(c(0,counts)),col="blue",lwd=2)
  lines(rep(mu_c,2),range(c(0,counts)),lty=2,col="red")
  E.mu_c=mean(ests$mu_c[keep])
  se.mu_c=sqrt(var(ests$mu_c[keep]))
  RMSE.mu_c=sqrt(mean((ests$mu_c[keep]-mu_c)^2))
  cv.mu_c=se.mu_c/E.mu_c
  pcbias.mu_c=100*(E.mu_c-mu_c)/mu_c
  
  #counts=hist(ests$E2,plot=FALSE)$counts
  #hist(ests$E2,main="E[down]",xlim=range(E2,ests$E2))
  #lines(rep(mean(ests$E2),2),range(c(0,counts)),col="blue",lwd=2)
  #lines(rep(E2,2),range(c(0,counts)),lty=2,col="red")
  
  # Sample sizes
  counts=hist(ests$n1[keep],plot=FALSE)$counts
  hist(ests$n1[keep],main="(d)",xlab=expression(n[1]),xlim=range(ests$n1[keep],ests$n2[keep]),nclass=nclass)
  lines(rep(mean(ests$n1[keep]),2),range(c(0,counts)),col="blue",lwd=2)
  mean(ests$n1[keep])
  
  counts=hist(ests$n2[keep],plot=FALSE)$counts
  hist(ests$n2[keep],main="(e)",xlab=expression(n[2]),xlim=range(ests$n1[keep],ests$n2[keep]),nclass=nclass)
  lines(rep(mean(ests$n2[keep]),2),range(c(0,counts)),col="blue",lwd=2)
  mean(ests$n2[keep])
  
  counts=hist(ests$m[keep],plot=FALSE)$counts
  hist(ests$m[keep],main="(f)",xlab=expression(m),xlim=range(ests$m[keep]),nclass=nclass)
  lines(rep(mean(ests$m[keep]),2),range(c(0,counts)),col="blue",lwd=2)
  mean(ests$m[keep])
  
  if(var(ests$sigma[keep])==0 & var(ests$kappa[keep])==0) cormat=cov2cor(cov(ests[keep,c("Dhat","n2","n1","m")]))
  else if(var(ests$sigma[keep])==0 & var(ests$kappa[keep])>0) cormat=cov2cor(cov(ests[keep,c("Dhat","kappa","n2","n1","m")]))
  else if(var(ests$sigma[keep])>0 & var(ests$kappa[keep])==0) cormat=cov2cor(cov(ests[keep,c("Dhat","sigma","n2","n1","m")]))
  else cormat=cov2cor(cov(ests[keep,c("Dhat","kappa","sigma","n2","n1","m")]))
  #windows()
  #pairs(ests[,c("n1","n2","m")])
  return(list(
    pcbias.Dhat=pcbias.Dhat,
    cv.Dhat=cv.Dhat,
    pcbias.Dhat.CI=pcbias.Dhat.CI,
    pcbias.kappa=pcbias.kappa,
    cv.kappa=cv.kappa,
    pcbias.sigma=pcbias.sigma,
    cv.sigma=cv.sigma,
    cor=cormat
  ))
}

#' @title Simulates 2-plane data & fits LCE and CCR models to this
#'
#' @description
#'  Simulates 2-plane data & fits LCE and CCR models to this
#'  
#' @param D.2D animal density (number per sq km)
#' @param L transect length (km)
#' @param w searched strip width (km)
#' @param b buffered strip width (km)
#' @param sigmarate the Brownian movement rate (sigma) parameter
#' @param k the lag between observers (seconds)
#' @param planespd Observer speed in km/sec
#' @param kappa mean time animals on surface (seconds)
#' @param tau mean dive cycle leength (seconds)
#' @param p vector of each observers' detection prob for available animals
#' @param movement list with components \code{forward} and \code{sideways}, which should be
#' TRUE or FALSE depending on what sort of movement you want to simulate with. 
#' Note: Estimation does not currently work when simulating with \code{forward=FALSE}.
#' @param fix.N Boolean variable, which if TRUE simulates with the same abundance, equal
#' to the expected, abundance, each simulation rep.
#' @param En Ignored if NULL, else the expected number of detections by each observer. 
#' Function sets transect length, \code{L}, to achieve this
#' @param Nsim Number of simulations to do
#' @param writeout Boolean, which if TRUE, causes results to be written to file and
#' the name of the file (which encodes the simulation parameters) to be returned
#' invisibly.
#' @param seed Integer seed which can be used for repoducible results
#' @param simthod Character; either `MLE` (for max likelihood simulator) or `Palm`
#' for CCR simulator from package \code{palm}.
#' @param control.opt passed to \code{optim}
#' @param hessian Boolean controlling whether or not to return Hessian 
#' (and interval estimates for MLE method).
#' @param adj.mvt If TRUE, MLE method takes account of randomness in lag between
#' encounters; if FALSE, it sets these all to the lag, \code{k}. 
#' @param ft.normal if TRUE, approximates Brownian motion hitting time with a 
#' normal distribution when estimation with MLE method, else uses exact expression
#' for Brownian hitting time
#' @param sim.ft.normal if TRUE, approximates Brownian motion hitting time with a 
#' normal distribution when simulating, else uses exact expression
#' for Brownian hitting time in simulation. 
#' NOTE: Exact simulator seems to have too low variance; with fast observers 
#' relative to animals, normal approx is very good (see function \code{dft} and
#' code in file \code{ft_plot.R}).
#' @param progbar If TRUE, puts a progress bar on screen so you can monitor simulation
#' progress.
dosim = function(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=c(1,1),movement=list(forward=TRUE,sideways=TRUE),
                 fix.N=TRUE,En=NULL,Nsim=100,writeout=TRUE,seed=1,simethod="MLE",control.opt=control.opt,
                 hessian=TRUE,adj.mvt=FALSE,ft.normal=FALSE,sim.ft.normal=FALSE,progbar=TRUE) {
  
  # Create progress bar
  if(progbar) pb <- tkProgressBar(title=paste("Function dosim Progress (Nsim=",Nsim,")",sep=""), min=0, max=Nsim, width=400)
  
  dmax.t = (b-w)/planespd
  ps = p.t(kappa,tau,p,sigmarate,k,dmax.t,planespd,w,adj.mvt=adj.mvt,io=movement$sideways) # capture history probabilities
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
  sigmapalm = sigmarate2sigmapalm(sigmarate,k)
  
  set.seed(seed) # initialise random number sequence (for repeatability)
  #  skip=c()
  startime=date()
  for(sim in 1:Nsim) {
    if(simethod=="MLE") {
      sdat = sim.2plane(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=c(1,1),movement=movement,fix.N=fix.N,
                        sim.ft.normal=sim.ft.normal)
    } else if(simethod=="Palm") {
      Ntype = "RandomN"
      palmdat <- sim.twocamera(c(D.2D=D.2D, kappa=kappa, sigma=sigmapalm),d=L, w=w, b=b, l=k, tau=tau)
      sdat = Palm2mleData(palmdat$points,palmdat$sibling.list$cameras,d=L,l=k,w=w,b=b)
    } else stop("simethod must be 'MLE' or 'Palm'.")
    
    # MLE
    mlefit<-segfit(sdat,D.2D,E1=kappa,Ec=tau,sigmarate=sigmarate,planespd=planespd,p=c(1,1),
                   sigma.mult=sigma.mult,control.opt=control.opt,method="BFGS",estimate=estimate,
                   set.parscale=TRUE,io=TRUE,Dbound=NULL,hessian=hessian,adj.mvt=adj.mvt,ft.normal=ft.normal)
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
    palmests$sigmarate[sim]=sigmarate2sigmapalm(coef(palmfit)["sigma"],sdat$k) # Convert to Brownian sigmarate
    
    # Sample size statistics
    mlests$n1[sim]=length(sdat$y1)
    mlests$n2[sim]=length(sdat$y2)
    palmests$n1[sim]=length(sdat$y1)
    palmests$n2[sim]=length(sdat$y2)
    if(simethod=="MLE") mlests$m[sim]=sdat$n11 # record number of actual duplicates
    
    # Progress bar stuff
    if(progbar) setTkProgressBar(pb, sim, label=paste( round(sim/Nsim*100, 0),"% done"))
    
  }  # End sim loop
  
  # Close progress bar
  if(progbar) close(pb)
  
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


#' @title Gets, and summarises simulation results
#'
#' @description
#'  Reads simulation results from \code{dosim}, using appropriate file name 
#'  (see details in code) and returns summary of bias, cv, coverage etc. 
#'  
#' @param fn Nam of the file with the output from \code(dosim).
#' You should use the filename returned by \code(dosim);  see \code(dosim)
#' code for details. 
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
  pc.bias.gamma = 100*(mean(sim[[1]]$gamma.est[!bad])-gamma)/gamma
  pc.cv.gamma = 100*sqrt(var(sim[[1]]$gamma.est[!bad]))/mean(sim[[1]]$gamma.est[!bad])
  cover.gamma = sum(sim[[1]]$gamma.inci[!bad])/(Nsim-nbad)
  pc.bias.sigmarate = 100*(mean(sim[[1]]$sigmarate.est[!bad])-sigmarate)/sigmarate
  pc.cv.sigmarate = 100*sqrt(var(sim[[1]]$sigmarate.est[!bad]))/mean(sim[[1]]$sigmarate.est[!bad])
  cover.sigmarate = sum(sim[[1]]$sigmarate.inci[!bad])/(Nsim-nbad)
  
  pc.bias.palm = 100*(mean(sim[[2]]$Dhat[!bad])-D.2D)/D.2D
  pc.cv.palm = 100*sqrt(var(sim[[2]]$Dhat[!bad]))/mean(sim[[2]]$Dhat[!bad])
  
  Dhat.cor = cor(sim[[1]]$D.est[!bad],sim[[2]]$Dhat[!bad])
  
  n1 = mean(sim[[1]]$n1[!bad])
  n2 = mean(sim[[1]]$n2[!bad])
  m = mean(sim[[1]]$m[!bad])
  
  mean.gamma.est = mean(sim[[1]]$gamma.est[!bad])
  mean.sigmarate.est = mean(sim[[1]]$sigmarate.est[!bad])
  mean.sehat.Dhat = mean(sim[[1]]$D.se[!bad])
  
  out = data.frame(Nsim=Nsim,gamma=gamma,k=k,sigmarate=sigmarate,avg.mps=getspeed(sigmarate,k)*1000,D.2D=D.2D,
                   pc.bias.mle=pc.bias.mle,pc.cv.mle=pc.cv.mle,cover.mle=cover.mle,
                   pc.bias.palm=pc.bias.palm,pc.cv.palm=pc.cv.palm,
                   pc.bias.gamma=pc.bias.gamma,pc.cv.gamma=pc.cv.gamma,cover.gamma=cover.gamma,
                   pc.bias.sigmarate=pc.bias.sigmarate,pc.cv.sigmarate=pc.cv.sigmarate,cover.sigmarate=cover.sigmarate,
                   Dhat.cor=Dhat.cor,
                   n1=n1,n2=n2,m=m,
                   mean.gamma.est=mean.gamma.est,
                   mean.sigmarate.est=mean.sigmarate.est,
                   mean.sehat.Dhat=mean.sehat.Dhat,
                   nbadD=nbadD,nbadse=nbadse)
  
  return(out)
}  



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

