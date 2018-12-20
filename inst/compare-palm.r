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

sigma=0.15*2/sqrt(248)   #  sigma from Ben's paper * sqrt(2).
b = 2
w = 0.125
sigma.mult = (b-w)/sigma

# Convert data to format required by segfit:
sdat = Palm2mleData(porpoise.data$points,porpoise.data$cameras,porpoise.data$d,porpoise.data$l,porpoise.data$w,b)

# Look at numbers within segments:
pdf(h=4,w=8,file="./inst/results/porpoise_pairingdbn.pdf")
segmentplot(sdat,planespd)
dev.off()

control.opt=list(trace=0,maxit=1000)
estimate=c("D","sigma","E1") # parameters to estimate

mlefit<-segfit(sdat,D.2D,E1=kappa,Ec=tau,sigmarate=sigmarate,planespd=planespd,p=c(1,1),
               sigma.mult=sigma.mult,control.opt=control.opt,method="BFGS",estimate=estimate,
               set.parscale=TRUE,io=TRUE,Dbound=NULL,hessian=TRUE,adj.mvt=TRUE,ft.normal=FALSE)

mlefit
cov2cor(mlefit$vcv)
sigma.palm = sigmarate2sigmapalm(mlefit$sigmarate["est"],248)
getspeed(mlefit$sigmarate["est"],248)*1000
  
# Save for use in paper:
saveRDS(mlefit,file="./inst/results/survey.mle.Rds")
