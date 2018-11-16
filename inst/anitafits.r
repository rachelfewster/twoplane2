library(twoplane)
data('anita')

A=abs(diff(range(anita$x1,anita$x2)))*abs(diff(range(anita$s1,anita$s2)));A

sdat=anita125
# convert km to plane seconds:
sdat$s1=sdat$s1/planespd
sdat$s2=sdat$s2/planespd
sdat$x1=sdat$x1/planespd
sdat$x2=sdat$x2/planespd
n1=length(sdat$s1);n1
n2=length(sdat$s2);n2
sdat$tL=abs(diff(range(c(0,sdat$s1,sdat$s2))));sdat$tL
ylim=range(sdat$x1,sdat$x2);ylim
xlim=range(sdat$s1,sdat$s2);xlim
D.line.t=n1/sdat$tL;D.line.t
sdat$tw=xlim[2]

# individual and plane speeds
sdat$k=500
Ec=120
E1=0.25*Ec
sigmaknots=20
sigmarate=sigmaknots*nm2km/60/60 # std dev of movement distance per second (in m)
sigma=sigmarate*sqrt(k) # std dev of movement distance after k seconds (in m)
sigma*1000*2 # 95% within this many m of start point after k seconds
planeknots=100 # observer speed in knots
planespd=planeknots*nm2km/(60^2);planespd # observer speed in km/sec

cutstretch=1000
dmax.t=sigma.mult*sigma/planespd;dmax.t*cutstretch


cuts<-sort(unique(c(0,sdat$tL,segmentize(sdat$s1,sdat$s2,dmax=dmax.t))))
plot(sdat$s1,sdat$x1,xlim=xlim,ylim=ylim,col="blue",xlab="Along-transect distance",ylab="Perpendicular distance")
points(sdat$s2,sdat$x2,col="red")
segments(cuts,rep(-dmax.t,length(cuts)),cuts,rep(dmax.t,length(cuts)))
length(cuts); length(sdat$s1); length(sdat$s2)

control.opt=list(trace=5,maxit=1000)
fitio<-segfit(sdat,D.line.t,E1=E1,Ec=Ec,sigma=sigma,planespd=planespd,p=c(1,1),sigma.mult=sigma.mult,
              control.opt=control.opt,method="BFGS",estimate=estimate,set.parscale=TRUE,
              io=TRUE,Dbound=NULL,hessian=TRUE,cutstretch=cutstretch)
fitio
intest=logn.seci(log(fitio$D),sqrt(solve(fitio$hessian)[1,1]))

fitio$D/planespd
intest$se/planespd
c(intest$lower,intest$upper)/planespd

Nhat=fitio$D*sdat$tL;Nhat
intest$se*sdat$tL
c(intest$lower,intest$upper)*sdat$tL

Dhat=Nhat/A;Dhat

save(fitio,file="./inst/fitio.anita.RData")
load("./inst/fitio.anita.RData")
