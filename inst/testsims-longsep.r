library(twoplane)
# requires HiddenMarkov, functional, Rcpp, boot, expm, Matrix
library(mvtnorm)
# library(devtools); install_github("b-steve/palm")
library(palm)



nm2km=1.852 # multiplier to convert nautical miles to kilometres

E2=110-94
E1=94
E=c(E1,E2)
Ec=sum(E) # mean dive cycle length in seconds
p.up=E1/Ec # proportion of time up


L=1100;w=0.25 # length and width of strip in km
N=1.05*L*w
k=248 # time separation of observers in seconds

#=================================
# Just temp to make things faster:
N=1.05*L*w/2 
k=20 # time separation of observers in seconds
#=================================


planeknots=100 # observer speed in knots
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec


# Set diffusion coefficient in 1D (f(x,t)=N(x; mu=0,var=sigma^2*k)), so that about 95%
# of animals have moved a distance of no more than 2*sigma*sqrt(k) km in time k
# Hiby & Lovell use animal speed of 1.5m/s. Here's corresponding sigma such that about 95% are no faster than this

sigma=0.34
sigmarate=sigma/sqrt(k)
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

Nsim=30

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
  
  # fit accounting for leakage of animals in and out of strip
  segiotime[sim]=system.time(fitio<-segfit(sdat,D.line.t,E1=E1,Ec=Ec,sigmarate=sigmarate,planespd=planespd,p=c(1,1),sigma.mult=sigma.mult,
                                           control.opt=control.opt,method="BFGS",estimate=estimate,set.parscale=TRUE,
                                           io=iomvt,Dbound=NULL,hessian=TRUE))[3]
  
  # Palm
  nspptime[sim]=system.time(est.nspp<-twoplane.fit(sdat,planespd,Ec,sigma,sigma.mult=5,trace=FALSE))[3]
  
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
    estsio$se[sim]=intest$se/(w*planespd)
    estsio$inci[sim]=(intest$lower/(w*planespd)<=Dstrip & Dstrip<=intest$upper/(w*planespd))
  } else skip=c(skip,sim)
  estsio$Dhat[sim]=fitio$D/(w*planespd)
  estsio$E1[sim]=fitio$E[1]
  estsio$E2[sim]=fitio$E[2]
  estsio$sigma[sim]=fitio$sigmarate
  estsio$mu_c[sim]=fitio$mu_c
  
  estsnspp$Dhat[sim]=est.nspp[1]
  estsnspp$sigma[sim]=est.nspp[3]
  
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

simsum.nspp(estsnspp,true,nspptime)

Petersen.bias=round(100*(mean(ests.kd$Dhat)-D.line.t),1)
Petersen.pc.cv=round(100*sd(ests.kd$Dhat)/mean(ests.kd$Dhat),1)
Petersen.bias;Petersen.pc.cv
