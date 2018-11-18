# Simulation of "real" data scenario in Stevension e al. (2018)
library(twoplane)
library(palm)

planeknots=100 # observer speed in knots
nm2km=1.852 # multiplier to convert nautical miles to kilometres
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec


w=0.125*2 # width of strip in km from porpoise data
D = 1.05
Ec = 110 # mean dive cycle length
E1 = 94
k = 248
alpha = E1/Ec
L = 1100
sigmarate = 0.95/1000/(sqrt(2)*gamma(1)/gamma(0.5)) # consistent with avg speed of 0.95 m/s
sigma.mult=5
dmax.km=sigma.mult*sigmarate*sqrt(k) # max dist apart (in km) observations could be considered duplicates; 
dmax.time=(dmax.km/planespd)   # max time apart (in seconds) observations could be considered duplicates
tw=w/planespd # width of searched stip, in plane seconds
btw = tw+2*dmax.time # width of stip with buffer, in plane seconds
tw/btw

p=c(1, 1) # definitely see if available in searched strip
ps = p.t(E1,Ec,p,sigmarate,k,dmax.time,planespd,w/2) # capture history probabilities
p. = sum(ps) # prob detect

b = w + 2*dmax.km
En = L*(D*b*p.)
N=D*(L*b)

Dstrip=N/(L*b) # density in number per sq km
Dstrip.t=D*(planespd^2) # density in planespd units
D.line.t=Dstrip.t*b/planespd # density in planespd along LINE units (1-dimensional)
control.opt=list(trace=0,maxit=1000)

Nsim=30

seed = 1
# Do one scenaio with a few simulations to check it works:
tm   = system.time(testsim   <- dosim(alpha,Ec,k,w,sigmarate,planespd,D,En=En,sigma.mult=sigma.mult,
                                      seed=seed,Nsim=Nsim,writeout=FALSE,iomvt=FALSE))
tmvt = system.time(testsimvt <- dosim(alpha,Ec,k,w,sigmarate,planespd,D,En=En,sigma.mult=sigma.mult,
                                      seed=seed,Nsim=Nsim,writeout=FALSE,iomvt=TRUE))
harvestsim(alpha,k,sigmarate,D,En=En,Nsim=10,simresults=testsim)
harvestsim(alpha,k,sigmarate,D,En=En,Nsim=10,simresults=testsimvt)
