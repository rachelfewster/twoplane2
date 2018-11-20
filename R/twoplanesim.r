sim.2plane=function(N,L,w,sigmarate,k,planespd,p.up,E.c,p=c(1,1),
                    sigma.mult=6,movement=list(forward=TRUE,sideways=TRUE)){
  
  E1=p.up*E.c # expected time up
  E2=E.c-E1   # expected time down
  p1=p[1]
  p2=p[2]
  
  dmax.km=sigma.mult*sigmarate*sqrt(k) # max dist apart (in km) observations could be considered duplicates; 
  # Also width of buffer around strip.
  # Brownian mvt has sd of sigmarate*sqrt(time passed)
  dmax.time=(dmax.km/planespd)   # max time apart (in seconds) observations could be considered duplicates

  tL=L/planespd # transect length, in plane seconds
  tw=w/planespd # width of searched stip, in plane seconds
  btw = tw+2*dmax.time # width of stip with buffer, in plane seconds
  
  
  if(!movement$sideways) { # Here for no sideways movement.
    NumSimAnimals=round(N)
    l=sort(runif(NumSimAnimals,0,tL)) # generate animal locations (in plane seconds)
    #    lhoriz=rep(0,NumSimAnimals)  #  Place animals in centre of strip if there is no horizontal movement.
    lhoriz=runif(NumSimAnimals,-tw/2,tw/2)  #  Place animals uniformly in strip.
    if(movement$forward) {
      simt = rft(NumSimAnimals,k,planespd,log(sigmarate)) -k # deviation from lag of k
    } else { # no forward movement
      simt=rep(0,NumSimAnimals)
    }
    l2=l+simt # locations of animals (in plane time since start of transect, time measured in plane seconds) when 2nd observer passes
    l2horiz=lhoriz  #  Horizontal position when second plane passes. 
  }  # here for sideways movement
  else {
    #  Horizontal movement. Need to simulate more animals, within a wider strip with width w+2*dmax.km
#    Nhoriz=round(N*btw/tw)
    NumSimAnimals=round(N)
    l=sort(runif(NumSimAnimals,0,tL)) # generate animal locations (in plane seconds)
    lhoriz=runif(NumSimAnimals,-btw/2, btw/2)  #  Horizontal animal locations in plane seconds
    if(movement$forward) {
      simt = rft(NumSimAnimals,k,planespd,log(sigmarate)) -k # deviation from lag of k
    } else { # no forward movement
      simt=rep(0,NumSimAnimals)
    }
#    simthoriz=rtnorm(NumSimAnimals,mean=0,sd=sqrt(k)*sigmarate/planespd,xtrunc=btw/2)
    simthoriz=rnorm(NumSimAnimals, 0, sqrt(k)*sigmarate/planespd)
    l2=l+simt
    l2horiz=lhoriz+simthoriz
  }
  # simulate detection by leading observer:
  up.1=1*(runif(NumSimAnimals)<=(E1/E.c)) # random variable with 1=available, 0=unavailable; E1/E.c is stationary dbn prob of being up
  in.1=(lhoriz<=tw/2)*(lhoriz>=-tw/2)   #  Is it in the strip?
  see.1=(runif(NumSimAnimals)<=p[1])*up.1*in.1 # binary variable indicating whether or not plane 1 sees
  
  in.2=(l2horiz<=tw/2)*(l2horiz>=-tw/2)   #  Is it in the strip.
  p.01.2=p.11.2=rep(0,NumSimAnimals)
  t12=k+simt    # seconds elapsing between 1st and 2nd observer passing animals
  q11=-1/E1
  q22=-1/E2
  Qmat=matrix(c(q11,-q22,-q11,q22),nrow=2)
  for(i in 1:NumSimAnimals) {
    p.01.2[i]=p.omega.t(t12[i], idbn=c(up.1[i],(1-up.1[i])), p1, p2, Qmat, omega=01)  # Prob(2nd sees | state when 1st passes)
    p.11.2[i]=p.omega.t(t12[i], idbn=c(up.1[i],(1-up.1[i])), p1, p2, Qmat, omega=11)  # Prob(2nd sees | state when 1st passes)
  }
##  if(!movement$sideways) { # no movement in/out of stip
##    for(i in 1:NumSimAnimals) {
##      p.01.2[i]=p.omega.t(t12[i], idbn=c(up.1[i],(1-up.1[i])), p1, p2, Qmat, omega=01)  # Prob(2nd sees | state when 1st passes)
##      p.11.2[i]=p.omega.t(t12[i], idbn=c(up.1[i],(1-up.1[i])), p1, p2, Qmat, omega=11)  # Prob(2nd sees | state when 1st passes)
##    }
##  } else { # movement in/out of stip
#### For debugging, just use k as time. Later remove line below and uncomment 2 lines down TPM = ...
###        TPM = make.inout.tpm(sigma=sigmarate*sqrt(k),dmax=dmax.km,w=w) # in-out transition probability matrix
##    for(i in 1:NumSimAnimals) {
##      TPM = make.inout.tpm(sigma=sigmarate*sqrt(k+simt[i]),dmax=dmax.km,w=w) # in-out transition probability matrix
##      idbn = c(up.1[i]*in.1[i],up.1[i]*(1-in.1[i]),(1-up.1[i])*in.1[i],(1-up.1[i])*(1-in.1[i]))
##      p.01.2[i]=p.omega.t(t12[i], idbn=idbn, p1, p2, Qmat, omega=01,IO=TPM)  # Prob(2nd sees | state when 1st passes)
##      p.11.2[i]=p.omega.t(t12[i], idbn=idbn, p1, p2, Qmat, omega=11,IO=TPM)  # Prob(2nd sees | state when 1st passes)
##    }
##  }
  
  p.see.2.given.in = (p.01.2+p.11.2)
  see.2=(runif(NumSimAnimals)<=p.see.2.given.in)*in.2 # binary variable indicating whether or not plane 2 sees
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
  x1=lhoriz[sno1]
  x2=l2horiz[sno2]
  # mark duplicates:
  d1=which(is.element(sno1,dno))
  d2=which(is.element(sno2,dno))
  n1=length(s1)
  n2=length(s2)
  
  sim.t.area=tL*btw
  sim.t.D=N/sim.t.area
  sim.area=sim.t.area*planespd
  sim.D=N/sim.area
  
  simdat=list(s1=s1,s2=s2,x1=x1,x2=x2,tL=tL,tw=tw,btw=btw,k=k,
              sno1=sno1,sno2=sno2,d1=d1,d2=d2,n1=n1,n2=n2,n11=m,n=n1+n2-m,alls1=l,allx1=lhoriz,alls2=l2,allx2=l2horiz,dups=dups,
              sim.t.area=sim.t.area,sim.t.D=sim.t.D,sim.area=sim.area,sim.D=sim.D,
              p.01.2=p.01.2,p.11.2=p.11.2
              ,see.1=see.1, see.2=see.2)
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



plot.ests=function(ests,D,E1,sigma, mu_c, sigma.est=TRUE,nclass=20,keep=NULL)
#-------------------------------------------------------------------------------
# NB: Not sure this is an up-to-date function.
# Plots histograms of simulated estimates and parameters and prints bias, CV
# and parameter correlations.
#
# Inputs: 
# ------
#  ests: dataframe with columns 
#        Dhat : Density estimate 
#        E1   : Estimated expected time available
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
  
  counts=hist(ests$E1[keep],plot=FALSE)$counts
  hist(ests$E1[keep],main="(b)",xlab=expression(paste(hat(E),"[up]",sep="")),xlim=range(E1,ests$E1[keep]),nclass=nclass)
  lines(rep(mean(ests$E1[keep]),2),range(c(0,counts)),col="blue",lwd=2)
  lines(rep(E1,2),range(c(0,counts)),lty=2,col="red")
  E.E1=mean(ests$E1[keep])
  se.E1=sqrt(var(ests$E1[keep]))
  RMSE.E1=sqrt(mean((ests$E1[keep]-E1)^2))
  cv.E1=se.E1/E.E1
  pcbias.E1=100*(E.E1-E1)/E1
  
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
  
  if(var(ests$sigma[keep])==0 & var(ests$E1[keep])==0) cormat=cov2cor(cov(ests[keep,c("Dhat","n2","n1","m")]))
  else if(var(ests$sigma[keep])==0 & var(ests$E1[keep])>0) cormat=cov2cor(cov(ests[keep,c("Dhat","E1","n2","n1","m")]))
  else if(var(ests$sigma[keep])>0 & var(ests$E1[keep])==0) cormat=cov2cor(cov(ests[keep,c("Dhat","sigma","n2","n1","m")]))
  else cormat=cov2cor(cov(ests[keep,c("Dhat","E1","sigma","n2","n1","m")]))
  #windows()
  #pairs(ests[,c("n1","n2","m")])
  return(list(
    pcbias.Dhat=pcbias.Dhat,
    cv.Dhat=cv.Dhat,
    pcbias.Dhat.CI=pcbias.Dhat.CI,
    pcbias.E1=pcbias.E1,
    cv.E1=cv.E1,
    pcbias.sigma=pcbias.sigma,
    cv.sigma=cv.sigma,
    cor=cormat
  ))
}

