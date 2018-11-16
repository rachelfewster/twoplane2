library(twoplane)
# requires HiddenMarkov, functional, Rcpp, boot, expm, Matrix
library(mvtnorm)
# library(devtools); install_github("b-steve/palm")
library(palm)

#  Compare MLE to palm library using simulated porpoise data

data("porpoise")

nm2km=1.852 # multiplier to convert nautical miles to kilometres

E2=14
E1=86
E=c(E1,E2)
Ec=sum(E) # mean dive cycle length in seconds
p.up=E1/Ec # proportion of time up

L=porpoise.data$d
w=porpoise.data$w*2 # length and width of strip in km

planeknots=100 # observer speed in knots   CHECK THIS
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec

k=248 # time separation of observers in seconds

sigma=0.15   #  Value for 248s separation.
sigmarate=sigma/sqrt(k)
sigma.mult=8 # multiplier to set bound for maximum distance moved

p.see.up=c(1,1)
#p.see.up=c(0.8,0.8) # prob(see|up) for each observer

# prob up in k seconds, given up now:
p.k=calc.p.avail(k,E1,Ec,p=c(1,1),sigma,planespd);p.k
round(p.k/p.up*100,1)-100 # % inflation in detection prob due to recent availability

N=1.05*(L*w)

Dstrip=N/(L*w) # density in number per sq km
Dstrip.t=Dstrip*(planespd^2) # density in planespd units
D.line.t=Dstrip.t*w/planespd # density in planespd along LINE units (1-dimensional)
Dbound=list(lower=-5*abs(log(D.line.t)),upper=5*abs(log(D.line.t))) # need bounds only if doing 1-dim estimation
method="Nelder-Mead" # this is ignored if doing 1-dim estimation
control.opt=list(trace=0,maxit=1000)


#  Split into camera 1 and camera 2 and convert from km to plane seconds. 
s1=porpoise.data$points[porpoise.data$cameras==1]/planespd
s2=porpoise.data$points[porpoise.data$cameras==2]/planespd


Nsim=1


estsio=ests=data.frame(Dhat=rep(0,Nsim),E1=rep(0,Nsim),E2=rep(0,Nsim),sigma=rep(0,Nsim),
                       n1=rep(0,Nsim),n2=rep(0,Nsim),m=rep(0,Nsim),mu_c=rep(0,Nsim),
                       se=rep(0,Nsim),inci=rep(0,Nsim))
estsnspp=estsna=data.frame(Dhat=rep(0,Nsim),n1=rep(0,Nsim),n2=rep(0,Nsim),sigma=rep(0,Nsim))

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

seed=654321 # arbitrary fixed seed - for generation of exactly same set of data on separate occasions
set.seed(seed) # initialise random number sequence (for repeatability)

sdat=list()
sdat$s1=s1
sdat$s2=s2
sdat$k=k
sdat$tL=L/planespd
sdat$tw=w/planespd

sim=1
# fit accounting for leakage of animals in and out of strip
segiotime[1]=system.time(fitio<-segfit(sdat,D.line.t,E1=E1,Ec=Ec,sigmarate=sigmarate,planespd=planespd,p=p.see.up,sigma.mult=sigma.mult,
                                         control.opt=control.opt,method="BFGS",estimate=estimate,set.parscale=TRUE,
                                         io=iomvt,Dbound=NULL,hessian=TRUE,krtest=FALSE))[3]
vcv.io = solve(fitio$hessian)
corr.io = cov2cor(vcv.io)
ses = sqrt(diag(vcv.io))
cvio = data.frame(D=ses[1]/(fitio$D/(w*planespd)),E1=fitio$E[2]/ses[2],sigmarate=fitio$sigmarate/ses[3])

#nspptime[sim]=system.time(est.nspp<-twoplane.fit(sdat,planespd,Ec,sigma,sigma.mult=sigma.mult,trace=FALSE,all=TRUE))[3]
pdat = format4palm(sdat,planespd,Ec,sigma,sigma.mult=sigma.mult)
nspptime[sim]=system.time(palmfit<-fit.twocamera(pdat$points,pdat$cameras,pdat$d,pdat$w,pdat$b,pdat$l,pdat$tau,pdat$R,trace=FALSE))[3]
est.palm=coef(palmfit)
#  Convert 'activity centre' sigma of Palm into our sigma by *sqrt(2)  and convert to sigmarate.
est.palm[3]=est.palm[3]*sqrt(2)/sqrt(sdat$k)
# Get palm bootstrap estimates
bootest.palm = boot.palm(palmfit,100)
bsum.palm = summary(bootest.palm)
cv.palm =  bsum.palm[[2]]/ bsum.palm[[1]]

Dhat=fitio$D/(w*planespd)
E1=fitio$E[1]
sigma=fitio$sigmarate
survey.mle = data.frame(est=c(Dhat,E1,sigma),
                 lo=exp(c(log(Dhat),log(E1),log(sigma))-1.96*ses),
                 hi=exp(c(log(Dhat),log(E1),log(sigma))+1.96*ses),
                 cv = rep(NA,3))
survey.mle$cv = (survey.mle$hi-survey.mle$lo)/(2*1.96)/survey.mle$est
row.names(survey.mle) = c("D","mu1","sigma")
# convert from sigmarate to sigma:
survey.mle["sigma",c("est","lo","hi")] = survey.mle["sigma",c("est","lo","hi")] *sqrt(k)
survey.mle
saveRDS(survey.mle,file="./inst/results/survey.mle.Rds")

estsnspp$Dhat[sim]=est.palm[1]
estsnspp$sigma[sim]=est.palm[3]

estsio
estsnspp
true



#######################
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
