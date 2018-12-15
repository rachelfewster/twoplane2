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

