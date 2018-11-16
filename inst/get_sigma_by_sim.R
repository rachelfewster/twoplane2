
nm2km=1.852 # multiplier to convert nautical miles to kilometres

E2=14
E1=86
E=c(E1,E2)
Ec=sum(E) # mean dive cycle length in seconds
p.up=E1/Ec # proportion of time up


L=1100;w=0.25 # length and width of strip in km
N=1.05*L*w*10 # big so get lots of observations - so mean is not noisy
planeknots=100 # observer speed in knots
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec

k=248 # time separation of observers in seconds

# Set diffusion coefficient in 1D (f(x,t)=N(x; mu=0,var=sigma^2*k)), so that about 95%
# of animals have moved a distance of no more than 2*sigma*sqrt(k) km in time k
# Hiby & Lovell use animal speed of 1.5m/s. Here's corresponding sigma such that about 95% are no faster than this

sigma1 = 0.333 # generates mean speed of 1.5m/s
sigma2 = 0.208 # generates mean speed of 0.95m/s
sigma3 = 0.1795 # what we get from MLE fitting to Ben's "real" data
sigma4 = 0.145 # what we get from MLE fitting to Ben's "real" data
sigmas = c(sigma1,sigma2,sigma3,sigma4)
names(sigmas)=c("Hiby","Westgate","mle","palm")
meanspd = sigmas*0

for(i in 1:4) {
  sigmarate=sigmas[i]/sqrt(k)
  sigma.mult=5 # multiplier to set bound for maximum km moved in k seconds (which is sigma.mult*sigma)
  
  p.see.up=c(1,1) # prob(see|up) for each observer
  
  # prob up in k seconds, given up now:
  p.k=calc.p.avail(k,E1,Ec,p=c(1,1),sigmarate,planespd);p.k
  
  #halfw=(w/2)/planespd # strip half width in planespd units
  Dstrip=N/(L*w) # density in number per sq km
  Dstrip.t=Dstrip*(planespd^2) # density in planespd units
  D.line.t=Dstrip.t*w/planespd # density in planespd along LINE units (1-dimensional)
  Dbound=list(lower=-5*abs(log(D.line.t)),upper=5*abs(log(D.line.t))) # need bounds only if doing 1-dim estimation
  method="Nelder-Mead" # this is ignored if doing 1-dim estimation
  control.opt=list(trace=0,maxit=1000)
  
  sdat=sim.2plane(N,L,w,sigmarate,k,planespd,p.up,Ec,p=p.see.up,
                  sigma.mult=sigma.mult,movement=list(forward=TRUE,sideways=iomvt,correct=FALSE),
                  cdfnormvt=NULL)
  
  mvs = sdat$alls1-sdat$alls2
#  hist(diff(mvs)*planespd*1000/k)
  n =  length(sdat$alls1);n
  meanspd[i] = mean(abs(diff(mvs)*planespd*1000/k))
}

speeds = data.frame(sigma=sigmas,meanspd=meanspd)
row.names(speeds) = names(sigmas)
saveRDS(speeds,file="./inst/results/speeds.Rds")
