library(twoplane)
# requires HiddenMarkov, functional, Rcpp, boot, expm, Matrix
library(mvtnorm)
# library(devtools); install_github("b-steve/palm")
library(palm)

#  Compare MLE to palm library using simulated porpoise data

data("porpoise")

nm2km=1.852 # multiplier to convert nautical miles to kilometres

tau = 110
gamma = 86/110
kappa = gamma*tau

planeknots=100 # observer speed in knots   CHECK THIS
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec

sigma=0.15*sqrt(2)   #  sigma from Ben's paper * sqrt(2).
b = 2
sigma.mult = (b-sdat$w)/sigma

sdat = Palm2mleData(porpoise.data$points,porpoise.data$cameras,porpoise.data$d,porpoise.data$l,porpoise.data$w,b)

control.opt=list(trace=0,maxit=1000)
estimate=c("D","sigma","E1") # parameters to estimate

mlefit<-segfit(sdat,D.2D,E1=kappa,Ec=tau,sigmarate=sigmarate,planespd=planespd,p=c(1,1),
               sigma.mult=sigma.mult,control.opt=control.opt,method="BFGS",estimate=estimate,
               set.parscale=TRUE,io=TRUE,Dbound=NULL,hessian=TRUE)

mlefit


# Try with a range of taus:
taus = c(80,90,100,110,120,130,140)
mles = taus*0
for(i in 1:length(mles)) {
  mlefit<-segfit(sdat,D.2D,E1=gamma*taus[i],Ec=taus[i],sigmarate=sigmarate,planespd=planespd,p=c(1,1),
                 sigma.mult=sigma.mult,control.opt=control.opt,method="BFGS",estimate=estimate,
                 set.parscale=TRUE,io=TRUE,Dbound=NULL,hessian=TRUE)
  
  mles[i] = mlefit["Dhat"]
}
plot(taus,100*(mles-mles[4])/mles[4])

#=============== Not sure of stuff below here; might be redundant ==============







# prob up in k seconds, given up now:
#p.k=calc.p.avail(k,E1,Ec,p=c(1,1),sigma,planespd);p.k
idbn = c(1,0,0,0)
p.k = p.t(kappa,tau,p.see.up,sigmarate,k,dmax.time,planespd,halfw.dist=halfw.dist,idbn=idbn)$ch11
round(p.k/p.up*100,1)-100 # % inflation in detection prob due to recent availability
p.k.nomvt = p.t(kappa,tau,p.see.up,sigmarate,k,dmax.time,planespd,halfw.dist=halfw.dist,idbn=c(1,0),io=FALSE)$ch11
round(p.k.nomvt/p.up*100,1)-100 # % inflation in detection prob due to recent availability

N=1.05*(L*w)

Dstrip=N/(L*w) # density in number per sq km
Dstrip.t=Dstrip*(planespd^2) # density in planespd units
D.line.t=Dstrip.t*w/planespd # density in planespd along LINE units (1-dimensional)
method="Nelder-Mead" # this is ignored if doing 1-dim estimation
control.opt=list(trace=0,maxit=1000)
estimate=c("D","sigma","E1") # parameters to estimate


#  Split into camera 1 and camera 2 and convert from km to plane seconds. 
y1 = porpoise.data$points[porpoise.data$cameras==1]
s1=y1/planespd
y2 = porpoise.data$points[porpoise.data$cameras==2]
s2=y2/planespd

n1 = length(y1)
n2 = length(y2)
dists = as.matrix(dist(c(y1,y2)))[1:n1,(n1+1):(n1+n2)]
mins = apply(dists,1,min)
hist(mins[mins<70/1000])

sdat = list(s1=s1,s2=s2,k=k,dmax.t=dmax.time,tL=L/planespd,tw=w/planespd)
#k=dat$k
#dmax.t=sigma.mult*(sigmarate*sqrt(k))/planespd # max time apart (in seconds) observations could be considered duplicates
#s1=dat$s1;s2=dat$s2
#tL=dat$tL
#halfw=dat$tw/2

fitio<-segfit(sdat,D.line.t,E1=kappa,Ec=tau,sigmarate=sigmarate,planespd=planespd,p=c(1,1),sigma.mult=sigma.mult,
              control.opt=control.opt,method="BFGS",estimate=estimate,set.parscale=TRUE,
              io=TRUE,Dbound=NULL,hessian=TRUE)
estsio=data.frame(Dhat=0,E1=0,E2=0,sigma=0,n1=0,n2=0,mu_c=0,se=0,lcl=0,ucl=0)

estsio$n1=length(sdat$s1)
estsio$n2=length(sdat$s2)

infmat=try(solve(fitio$hessian),silent=TRUE)
if(!inherits(infmat, "try-error")) {
  intest=logn.seci(log(fitio$D),sqrt(infmat[1,1]))
  estsio$se=intest$se/(2*b*planespd)
  estsio$lcl=intest$lower/(2*b*planespd)
  estsio$ucl=intest$upper/(2*b*planespd)
} else skip=c(skip,sim)
estsio$Dhat=fitio$D/(2*b*planespd)
estsio$E1=fitio$E[1]
estsio$E2=fitio$E[2]
estsio$sigma=fitio$sigmarate
estsio$mu_c=fitio$mu_c

estsio




# Try some simulation:
# ====================

k=20 # time separation of observers in seconds

tau = Ec = 110
gamma = p.up = 86/110
kappa = E1 = p.up*tau

L=porpoise.data$d
w=porpoise.data$w*2 # length and width of strip in km
halfw.dist = porpoise.data$w

planeknots=100 # observer speed in knots   CHECK THIS
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec

k=248 # time separation of observers in seconds

# Choose Ben's sigma or a sigma
#bens.sigma = FALSE
bens.sigma = TRUE
if(bens.sigma) {
  sigma.ben = sqrt(((0.95*k/1000)^2)/2)/2 * sqrt(2) # Ben's sigma * sqrt(2)
  sigmarate = sigma.ben/sqrt(k)
} else {
  animalspeed = 0.95/1000 # mean speed in km/sec
  sigmarate = animalspeed/(sqrt(2)*gamma(1)/gamma(0.5))
}
sigma = sigmarate*sqrt(k)
sigma.mult=5 # multiplier to set bound for maximum distance moved
dmax.km = sigma*sigma.mult
dmax.time = dmax.km/planespd
b = w+2*dmax.km
b/w

planeknots=100 # observer speed in knots
nm2km=1.852 # multiplier to convert nautical miles to kilometres
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec

D = 1.05
#D = 0.5

ps = p.t(kappa,tau,p.see.up,sigmarate,k,dmax.time,planespd,halfw.dist=halfw.dist)
p1 = ps$ch10 + ps$ch11

En = round(mean(n1,n2))
EN = round(En/p1)
EN/(L*b) # Expected density according to parameters passed to p.t() above
N = round(D*L*b)
EN;N

Nsim=100


seed = 1
#tm   = system.time(testsim   <- dosim(gamma,Ec,k,w,sigmarate,planespd,D,En=En,sigma.mult=sigma.mult,
#                                      seed=seed,Nsim=Nsim,writeout=FALSE,iomvt=FALSE,L=L))
tmvt = system.time(testsimvt <- dosim(gamma,Ec,k,w,sigmarate,planespd,D,En=En,sigma.mult=sigma.mult,
                                      seed=seed,Nsim=Nsim,writeout=FALSE,iomvt=TRUE,L=L))
#harvestsim(gamma,k,sigmarate,D,En=En,Nsim=10,simresults=testsim)
harvestsim(gamma,k,sigmarate,D,En=En,Nsim=10,simresults=testsimvt)


sims = list()
simno = 3
sims[[simno]]=testsimvt

harvestsim(gamma,k,sigmarate,D,En=En,Nsim=10,simresults=sims[[1]])
harvestsim(gamma,k,sigmarate,D,En=En,Nsim=10,simresults=sims[[2]])

hist(sims[[1]]$sim$mle$Dhat/sims[[simno]]$sim$palm$Dhat,xlab="MLE/Palm",main="Ratio of MLE to Palm Dhat")
hist(sims[[1]]$sim$mle$Dhat,xlab="Dhat",main="MLE")
points(mean(sims[[1]]$sim$mle$Dhat),0,pch="+",col="black",cex=2)
D=as.numeric(strsplit(strsplit(sims[[1]]$fn,"-")[[1]][5],"_")[[1]][2])
points(D,0,pch="+",col="red",cex=2)

hist(sims[[2]]$sim$mle$Dhat/sims[[simno]]$sim$palm$Dhat,xlab="MLE/Palm",main="Ratio of MLE to Palm Dhat")

D=as.numeric(strsplit(strsplit(sims[[simno]]$fn,"-")[[1]][5],"_")[[1]][2])
plot(sims[[simno]]$sim$mle$Dhat,sims[[simno]]$sim$palm$Dhat,xlab="MLE Dhat",ylab="Palm Dhat")
xylim = range(sims[[simno]]$sim$mle$Dhat,sims[[simno]]$sim$palm$Dhat)
lines(xylim,xylim,lty=2)
points(D,D,pch=19,col="red")
lines(rep(D,2),xylim,lty=2,col="red")
lines(xylim,rep(D,2),lty=2,col="red")

f#######################
#
#   Compare ..
library(binhf)
library(palm)

data("porpoise")

load("./data/porpoise-2D.RData")

pairDist=porpoise.data$points-shift(porpoise.data$points, 1, dir="left")   #  Distance between adjacent observations.

sel=which(abs(pairDist)<0.1 & porpoise.data$cameras!=shift(porpoise.data$cameras, 1, dir="left"))

pairDistHoriz=porpoise.2D[,2]-shift(porpoise.2D[,2], 1, dir="left")

#  -1 if the order is reversed, i.e. camera 2 before camera 1. Otherwise 1. 
ifReversed=(porpoise.data$cameras[sel]==1) - (porpoise.data$cameras[sel]==2)

#  Negate the distance when the first observation is by camera 2. 
pairDistClose=pairDist[sel]*ifReversed

pairDistHorizClose=pairDistHoriz[sel]*ifReversed

hist(pairDistClose)

hist(pairDistHorizClose)
