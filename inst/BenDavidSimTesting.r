library(palm)
library(twoplane)

nm2km=1.852 # multiplier to convert nautical miles to kilometres
planeknots=100 # observer speed in knots   CHECK THIS
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec

D = D.2D <- 1.05
## Time between cameras (seconds).
k = l <- 20
## Mean dive-cycle duration (seconds).
tau <- 110
## Mean duration of surface phase (seconds).
kappa <- 24
## Animal movement (roughly based on 0.95 m per s).
sigma <- gamma(0.5)*sqrt(((0.95*l/1000)^2)/2)/2
sigmarate = sigma*sqrt(2)/sqrt(k)
sigma.mult=5 # multiplier to set bound for maximum distance moved

p.see.up=c(1,1)

## Transect half-width.
halfw.dist = w <- 0.125
## Buffer distance
b <- w + sigma.mult*sigma
## Transect length.
L = d <- 1100

control.opt=list(trace=0,maxit=1000)
estimate=c("D","sigma","E1") # parameters to estimate

Nsim = 10
simests = data.frame(Dhat.MLE=rep(NA,Nsim),Dhat.Palm=rep(NA,Nsim),
                     kappa.MLE=rep(NA,Nsim),kappa.Palm=rep(NA,Nsim),
                     sigmarate.MLE=rep(NA,Nsim),sigma.Palm=rep(NA,Nsim))
set.seed(1)
for(sim in 1:Nsim) {
  #  sdat = sim.2plane(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=c(1,1),
  #                    movement=list(forward=TRUE,sideways=TRUE),fix.N=TRUE)
  #fitio<-segfit(sdat,D.2D,E1=kappa,Ec=tau,sigmarate=sigmarate,planespd=planespd,p=c(1,1),sigma.mult=sigma.mult,
  #              control.opt=control.opt,method="BFGS",estimate=estimate,set.parscale=TRUE,
  #              io=TRUE,Dbound=NULL,hessian=TRUE)
  bendata <- sim.twocamera(c(D.2D=D.2D, kappa=kappa, sigma=sigma),d=d, w=w, b=b, l=l, tau=tau)
  
  fit <- fit.twocamera(points=bendata$points, cameras=bendata$sibling.list$cameras, d=d, w=w, b=b, l=l, tau=tau, R=550)
  
  # Format for segfit, and fit
  ys = bendata$points
  obs = bendata$sibling.list$cameras
  y1 = ys[obs==1]
  y2 = ys[obs==2]
  bensdat = list(y1=y1,y2=y2,k=k,L=L,w=w,b=b)
  fitio<-segfit(bensdat,D.2D,E1=kappa,Ec=tau,sigmarate=sigmarate,planespd=planespd,p=c(1,1),sigma.mult=sigma.mult,
                control.opt=control.opt,method="BFGS",estimate=estimate,set.parscale=TRUE,
                io=TRUE,Dbound=NULL,hessian=TRUE)
  
  simests$Dhat.Palm[sim]=coef(fit)["D.2D"]
  simests$kappa.Palm[sim]=coef(fit)["kappa"]
  simests$sigma.Palm[sim]=coef(fit)["sigma"]*sqrt(2)/sqrt(l)
  
  simests$Dhat.MLE[sim]=fitio["Dhat"]
  simests$kappa.MLE[sim]=fitio["kappa"]
  simests$sigmarate.MLE[sim]=fitio["sigmarate"]
}

means = signif(apply(simests,2,mean),3)
means
100*(means[1:2]-D.2D)/D.2D
sqrt(var(simests$Dhat.MLE))/mean(simests$Dhat.MLE)
sqrt(var(simests$Dhat.Palm))/mean(simests$Dhat.Palm)



# Same, but using David's simulator instead of Ben's:
set.seed(1)
for(sim in 1:Nsim) {
  sdat = sim.2plane(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=c(1,1),
                      movement=list(forward=TRUE,sideways=TRUE),fix.N=TRUE)
  fitio<-segfit(sdat,D.2D,E1=kappa,Ec=tau,sigmarate=sigmarate,planespd=planespd,p=c(1,1),sigma.mult=sigma.mult,
                control.opt=control.opt,method="BFGS",estimate=estimate,set.parscale=TRUE,
                io=TRUE,Dbound=NULL,hessian=TRUE)
#  bendata <- sim.twocamera(c(D.2D=D.2D, kappa=kappa, sigma=sigma),d=d, w=w, b=b, l=l, tau=tau)
  
  fit = twoplane.fit(sdat,tau=tau,R=550,all=TRUE)
#  fit <- fit.twocamera(points=bendata$points, cameras=bendata$sibling.list$cameras, d=d, w=w, b=b, l=l, tau=tau, R=550)

  simests$Dhat.Palm[sim]=coef(fit)["D.2D"]
  simests$kappa.Palm[sim]=coef(fit)["kappa"]
  simests$sigma.Palm[sim]=coef(fit)["sigma"]*sqrt(2)/sqrt(sdat$k) # Convert to Brownian sigmarate
  
  simests$Dhat.MLE[sim]=fitio["Dhat"]
  simests$kappa.MLE[sim]=fitio["kappa"]
  simests$sigmarate.MLE[sim]=fitio["sigmarate"]
}

means = signif(apply(simests,2,mean),3)
means
100*(means[1:2]-D.2D)/D.2D
sqrt(var(simests$Dhat.MLE))/mean(simests$Dhat.MLE)
sqrt(var(simests$Dhat.Palm))/mean(simests$Dhat.Palm)

