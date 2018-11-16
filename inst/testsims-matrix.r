library(twoplane)
# requires HiddenMarkov, functional, Rcpp, boot, expm, Matrix
library(mvtnorm)
# library(devtools); install_github("b-steve/nspp")
library(nspp)

out=data.frame()
#   27 points
for(N in c(350,500,650)) {
  for(E1 in c(22,28,34)) {
    for(sigmaknots in c(1.8,2.5, 3.2)) {
      

nm2km=1.852 # multiplier to convert nautical miles to kilometres

E2=110  # expected length unavailable period (from Hiby&Lovell, 1998) in seconds
#E1=28   # expected length available period (from Hiby&Lovell, 1998) in seconds
E=c(E1,E2)
Ec=sum(E) # mean dive cycle length in seconds
p.up=E1/Ec # proportion of time up

#N=500 # number animals in strip
L=350;w=0.3 # length and width of strip in km
planeknots=100 # observer speed in knots
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec

k=10 # time separation of observers in seconds

# Set diffusion coefficient in 1D (f(x,t)=N(x; mu=0,var=sigma^2*k)), so that about 95%
# of animals have moved a distance of no more than 2*sigma*sqrt(k) nm in time k
#sigmaknots=2.5
sigmarate=sigmaknots*nm2km/60/60 # std dev of movement distance per second (in km)
sigma=sigmarate*sqrt(k) # std dev of movement distance after k seconds (in km)
sigma*1000*2 # 95% within this many metres of start point after k seconds
sigma.mult=5 # multiplier to set bound for maximum km moved in k seconds (which is sigma.mult*sigma)

p.see.up=c(1,1) # prob(see|up) for each observer

# prob up in k seconds, given up now:
p.k=calc.p.avail(k,E1,Ec,p=c(1,1),sigma,planespd);p.k
round(p.k/p.up*100,1)-100 # % inflation in detection prob due to recent availability

#halfw=(w/2)/planespd # strip half width in planespd units
Dstrip=N/(L*w) # density in number per sq km
Dstrip.t=Dstrip*(planespd^2) # density in planespd units
D.line.t=Dstrip.t*w/planespd # density in planespd along LINE units (1-dimensional)
Dbound=list(lower=-5*abs(log(D.line.t)),upper=5*abs(log(D.line.t))) # need bounds only if doing 1-dim estimation
method="Nelder-Mead" # this is ignored if doing 1-dim estimation
control.opt=list(trace=0,maxit=1000)

Nsim=100

fromfile = FALSE # FALSE if want to generate sim data, as opposed to retrive from file

if(fromfile) {
  load("./inst/IndepSims2.RData")
  Nsim=length(simlists)
}
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

# Comment out one of the two lines below
iomvt = TRUE # allow in-out movement in simulation and estimation
#iomvt = FALSE # do not allow in-out movement in simulation and estimation

make.sim.data=FALSE
#----------------------------------------------
if(make.sim.data) {
  seed=654321 # arbitrary fixed seed - for generation of exactly same set of data on separate occasions
  set.seed(seed) # initialise random number sequence (for repeatability)
  simlists=vector("list",length=Nsim)
  for(sim in 1:Nsim) {
    simlists[[sim]] = sim.2plane(N,L,w,sigmarate,k,planespd,p.up,Ec,p=p.see.up,
                                 sigma.mult=sigma.mult,movement=list(forward=TRUE,sideways=iomvt,
                                 correct=FALSE),cdfnormvt=NULL)
  }
}
## Uncomment line below to save simulated datasets
#save(simlists,file="./inst/IndepSims2.RData")
#----------------------------------------------

seed=654321 # arbitrary fixed seed - for generation of exactly same set of data on separate occasions
set.seed(seed) # initialise random number sequence (for repeatability)
skip=c()
startime=date()
for(sim in 1:Nsim) {
  if(fromfile) {
    sdat=simlists[[sim]]
  } else {
    sdat=sim.2plane(N,L,w,sigmarate,k,planespd,p.up,Ec,p=p.see.up,
                    sigma.mult=sigma.mult,movement=list(forward=TRUE,sideways=iomvt,correct=FALSE),
                    cdfnormvt=NULL)
  }
  
  if(checkdists) checkmvt(sigma,sdat,planespd,maxd=100)
  if(plot.sample) plot(sdat,w=w/2/planespd,L=L/planespd)
  if(plot.cuts) plotcuts(sdat,sigma.mult,planespd,w,L)
  
  
  # fit ignoring leakage of animals in and out of strip
##  segtime[sim]=system.time(fit<-segfit(sdat,D.line.t,E1=E1,Ec=Ec,sigma=sigma,planespd=planespd,p=c(1,1),sigma.mult=6,
##                                       control.opt=control.opt,method="BFGS",estimate=estimate,set.parscale=TRUE,
##                                       io=FALSE,Dbound=Dbound,hessian=TRUE))[3]

  # Debugging:
  #segfit(sdat,D.line.t,E1=E1,Ec=Ec,sigma=sigma,planespd=planespd,p=c(1,1),sigma.mult=sigma.mult,
  #       control.opt=control.opt,method="BFGS",estimate=estimate,set.parscale=TRUE,
  #       io=TRUE,Dbound=Dbound,hessian=TRUE)$D
  #sigma.mult=sigma.mult*1.2
  
  # fit accounting for leakage of animals in and out of strip
  segiotime[sim]=system.time(fitio<-segfit(sdat,D.line.t,E1=E1,Ec=Ec,sigmarate=sigmarate,planespd=planespd,p=c(1,1),sigma.mult=sigma.mult,
                                           control.opt=control.opt,method="BFGS",estimate=estimate,set.parscale=TRUE,
                                           io=iomvt,Dbound=NULL,hessian=TRUE))[3]
  
#  # fit with incorrect Ec
#  Ecerr=Ec*0.75
#  segtime[sim]=system.time(fit<-segfit(sdat,D.line.t,E1=E1,Ec=Ecerr,sigma=sigma,planespd=planespd,p=c(1,1),sigma.mult=sigma.mult,
#                                       control.opt=control.opt,method="BFGS",estimate=estimate,set.parscale=TRUE,
#                                       io=TRUE,Dbound=Dbound,hessian=TRUE))[3]
  
  #    est.na=segfit.noavail(sdat,D.line.t,E1=E1,Ec=Ec,sigma=sigma,planespd=planespd,p=c(1,1),
  #                          sigma.mult=6,control.opt=control.opt,method="BFGS",
  #                          set.parscale=TRUE,Dbound=Dbound,hessian=TRUE) 
  
#  nspptime[sim]=system.time(est.nspp<-twoplane.fit(sdat,planespd,Ec,sigma,sigma.mult=5,trace=FALSE))[3]
  
#  ests$n1[sim]=length(sdat$s1)
#  ests$n2[sim]=length(sdat$s2)
  estsio$n1[sim]=length(sdat$s1)
  estsio$n2[sim]=length(sdat$s2)
#  estsnspp$n1[sim]=length(sdat$s1)
#  estsnspp$n2[sim]=length(sdat$s2)
#  ests$m[sim]=sdat$n11 # record number of actual duplicates
  estsio$m[sim]=sdat$n11 # record number of actual duplicates
#  estsnspp$m[sim]=sdat$n11 # record number of actual duplicates
  
#  ests$Dhat[sim]=fit$D
#  ests$E2[sim]=fit$E[2]
#  ests$E1[sim]=fit$E[1]
#  ests$sigma[sim]=fit$sigma
#  ests$mu_c[sim]=fit$mu_c
#  intest=logn.seci(log(fit$D),sqrt(solve(fit$hessian)[1,1]))
#  ests$se[sim]=intest$se/(w*planespd)
#  ests$inci[sim]=(intest$lower/(w*planespd)<=Dstrip & Dstrip<=intest$upper/(w*planespd))
  
  infmat=try(solve(fitio$hessian),silent=TRUE)
  if(!inherits(infmat, "try-error")) {
    intest=logn.seci(log(fitio$D),sqrt(infmat[1,1]))
    estsio$se[sim]=intest$se/(w*planespd)
    estsio$inci[sim]=(intest$lower/(w*planespd)<=Dstrip & Dstrip<=intest$upper/(w*planespd))
  } else skip=c(skip,sim)
  estsio$Dhat[sim]=fitio$D/(w*planespd)
  estsio$E1[sim]=fitio$E[1]
  estsio$E2[sim]=fitio$E[2]
  estsio$sigma[sim]=fitio$sigmarate
  estsio$mu_c[sim]=fitio$mu_c
  #estsna$Dhat[sim]=est.na$D/(w*planespd)
  #estsna$p[sim]=est.na$p
  
#  estsnspp$Dhat[sim]=est.nspp$pars["D"]
#  estsnspp$sigma[sim]=est.nspp$pars["sigma"]
  
  # Petersen estimator for known recaptures, within strip:
  ests.kd$Dhat[sim] = ((sdat$n1+1)*(sdat$n2+1)/(sdat$n11+1) - 1)/(L/planespd)

  if(sim==1) cat("\nCounter: \n")
  if(sim %% 50 == 0) cat(sim,"\n")
  else if(sim %% 10 == 0) cat(sim)
  else if(sim %% 5 == 0) cat("+")
  else cat("-")
  
}  # End sim loop
endtime=date()

# look at distribution of actual duplicates
hist(estsio$m);mean(c(estsio$n1,estsio$n2))

# exclude stupid estimates (>= 100 times bigger than true density)
skip = c(skip,which(estsio$Dhat>10*true$D)); Nsim-length(skip)
if(length(skip)>0) simsum(estsio[-skip,],true,segiotime) else simsum(estsio,true,segiotime)
if(length(skip)>0) hist(estsio$Dhat[-skip]) else hist(estsio$Dhat)
# Not sure sigma is appropriate here:
#simsum.nspp(estsnspp,true.k,nspptime)
#true.k=true;true.k$sigma=sigma # deal with sigma for time k, not a rate

Petersen.bias=round(100*(mean(ests.kd$Dhat)-D.line.t),1)
Petersen.pc.cv=round(100*sd(ests.kd$Dhat)/mean(ests.kd$Dhat),1)
Petersen.bias;Petersen.pc.cv

#  For comparison to ML.  RMSE
(mean(( estsio$Dhat-Dstrip)^2)^0.5)
#  [1] 0.6466091   with no in/out movement, k=10
#  [1] 0.7399128   with in/out, k=10
#  Mean absolute error
mean(abs(estsio$Dhat-Dstrip))
#  [1] 0.5071363   with no in/out movement, k=10
#  [1] 0.6052138   with in/out, k=10
#  Percentage bias
(mean(estsio$Dhat)/mean(true$D))*100-100
hist(estsio$Dhat)

RMSE=(mean(( estsio$Dhat-true$D)^2)^0.5)
MAE=mean(abs(estsio$Dhat-true$D))
bias=(mean(estsio$Dhat)/mean(true$D))*100-100
out<-rbind(out, c(N, E1, sigmaknots, RMSE, MAE, bias))
}
}
}

#estsio=estsio[-skip,]
#save(estsio,ests.kd,file="./inst/IndepSims2ests.RData")
#load("./inst/IndepSims2ests.RData")
