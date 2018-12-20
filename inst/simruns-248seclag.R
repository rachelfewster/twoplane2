library(twoplane)
library(palm)
library(tcltk2) # for progress bar

nm2km=1.852 # multiplier to convert nautical miles to kilometres
planeknots=100 # observer speed in knots   
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec

D = D.2D <- 1.05
## Time between cameras (seconds).
k = l <- 248
## Mean dive-cycle duration (seconds).
tau <- 110
## Mean duration of surface phase (seconds).
kappa <- 0.86*tau
sigma_palm = 0.15 # estimated Palm-type sigma (in km) from porpoise data, with lag 248 seconds
sigmarate = sigmapalm2sigmarate(sigma_palm,lag=248) # Brownian sigmarate consistent with sigma_palm
sigma = sigmarate*sqrt(k) # Standard error of Brownian movement after k seconds
animalspd = getspeed(sigmarate,248)*1000 # average speed of animals in m/s, after 248 seconds
planespd/(animalspd/1000) # How much faster plane is going than average animal speed
#speed2sigmarate(.95/1000,248) # Bownian motion parameter that gives the observed speed over 248 seconds of .95 m/s
sigma.mult=5 # multiplier of sigma, to set bound for maximum distance moved in k seconds

p=c(1,1) # Probability see, given available in searched strips

## Transect half-width.
halfw.dist = w <- 0.125
## Buffer distance
b <- w + sigma.mult*sigma
## Transect length.
L = d <- 1100

control.opt=list(trace=0,maxit=1000)
estimate=c("D","sigma","E1") # parameters to estimate
movement = list(forward=TRUE,sideways=TRUE)

seed = 1 # for reproducibility
Nsim = 100 # Set number of simulations here

# Here to simulate wiht Palm simulator:
# adj.mvt=TRUE makes LCE estimator allow time between encouners to depend on animal movement (not be fixed at k)
# ft.normal=FALSE makes LCE estimator use exact expression for Brownian hitting time when adj.mvt==TRUE
# palmvt$file contains names of files that contain simulation results after dosim has finished
progbar = TRUE # set to FALSE if you don't want to generate a bar to track progress of simulations
palmtmvt = system.time(palmvt <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
                                       fix.N=TRUE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="Palm",
                                       control.opt=control.opt,adj.mvt=TRUE,ft.normal=FALSE,progbar=progbar))

# Here to simulate wiht LCE simulator:
# sim.ft.normal=TRUE makes LCE simulator use normal approxx for Brownian hitting time (exact simulator not quite working at present)
# mlemvt$file contains names of files that contain simulation results after dosim has finished
mletmvt = system.time(mlemvt <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
                                      fix.N=FALSE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="MLE",
                                      control.opt=control.opt,adj.mvt=TRUE,ft.normal=FALSE,sim.ft.normal=TRUE,
                                      progbar=progbar))

harvestsim(palmvt$file) # summarise Palm simulation results
harvestsim(mlemvt$file) # summarise LCE simulation results
msims = readRDS(mlemvt$file) # read all Palm simulation results
psims = readRDS(palmvt$file) # readall LCE simulation results

# Here's some old code, used to 
# (a) look at disribution of relative differences between Palm and LCE,
# (b) investigate bias when double and treble sample size
# Simulate from LCE simulator, EST with exact forward mvt dbn, SIM with normal approx forward mvt dbn
#Nsim = 100
#seed = 1
#tmvtadj.ft.sftnorm = system.time(mvtadj.ft.sftnorm <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
#                                                            fix.N=FALSE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="MLE",
#                                                            control.opt=control.opt,adj.mvt=TRUE,ft.normal=FALSE,sim.ft.normal=TRUE))
#harvestsim(mvtadj.ft.sftnorm$file)
#mvtadj.ft.sftnorm.sims=readRDS(mvtadj.ft.sftnorm$file)
#sim248.Nsim100 = mvtadj.ft.sftnorm.sims
## Save for use in paper:
#saveRDS(sim248.Nsim100,file="./inst/results/sim248.Nsim100.Rds")
#sim248.Nsim100 = readRDS("./inst/results/sim248.Nsim100.Rds")
#sim248.Nsim100.fn = mvtadj.ft.sftnorm$file
#saveRDS(sim248.Nsim100.fn,file="./inst/results/sim248.Nsim100.fn.Rds")
#pcDdiff = 100*(sim248.Nsim100$mle$D.est/sim248.Nsim100$palm$Dhat-1)
#hist(pcDdiff,nclass=20,probability=TRUE,xlab=expression((hat(D)[LCE]-hat(D)[Palm])/hat(D)[Palm]),main="")
#CI = quantile(pcDdiff,probs=c(0.025,0.975))
#segments(CI,rep(0,2),CI,rep(0.02,2),lwd=2)
#pcDdiff.obs = 18
#points(pcDdiff.obs,0,pch=19)
#2*(101-which(sort(c(pcDdiff,pcDdiff.obs))==pcDdiff.obs))
#simres.sim248.Nsim100=harvestsim(sim248.Nsim100.fn)
#saveRDS(simres.sim248.Nsim100,file="./inst/results/simres.sim248.Nsim100.Rds")

# Try with 2 times the sample size
#L = 2200
#Nsim = 50
#seed = 1
#tmvtadj.ft.sftnorm.2L = system.time(mvtadj.ft.sftnorm.2L <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
#                                                            fix.N=FALSE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="MLE",
#                                                            control.opt=control.opt,adj.mvt=TRUE,ft.normal=FALSE,sim.ft.normal=TRUE))
#harvestsim(mvtadj.ft.sftnorm.2L$file)
# Try with 3 times the sample size
#L = 3300
#Nsim = 300
#seed = 1
#tmvtadj.ft.sftnorm.3L = system.time(mvtadj.ft.sftnorm.3L <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
#                                                                  fix.N=FALSE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="MLE",
#                                                                  control.opt=control.opt,adj.mvt=TRUE,ft.normal=FALSE,sim.ft.normal=TRUE))
#harvestsim(mvtadj.ft.sftnorm.3L$file)

