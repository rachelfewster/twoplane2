library(twoplane)
library(palm)
library(tcltk2) # for progress bad

nm2km=1.852 # multiplier to convert nautical miles to kilometres
planeknots=100 # observer speed in knots   CHECK THIS
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec

D = D.2D <- 1.05
## Time between cameras (seconds).
k = l <- 248
## Mean dive-cycle duration (seconds).
tau <- 110
## Mean duration of surface phase (seconds).
kappa <- 86
sigma_palm = 0.15 # estimated Palm-type sigma (in km) from porpoise data, with lag 248 seconds
sigmarate = sigmapalm2sigmarate(sigma_palm,lag=248) # Brownian sigmarate consistent with sigma_palm
sigma = sigmarate*sqrt(k)
animalspd = getspeed(sigmarate,248)*1000 # in m/s
planespd/(animalspd/1000) # How much faster plane is going than average animal speed
#speed2sigmarate(.95/1000,248) # Bownian motion parameter that gives the observed speed over 248 seconds
sigma.mult=5 # multiplier to set bound for maximum distance moved

p=c(1,1)

## Transect half-width.
halfw.dist = w <- 0.125
## Buffer distance
b <- w + sigma.mult*sigma
## Transect length.
L = d <- 1100

control.opt=list(trace=0,maxit=1000)
estimate=c("D","sigma","E1") # parameters to estimate
movement = list(forward=TRUE,sideways=TRUE)

seed = 1
Nsim = 30
palmtmvt = system.time(testpalmvt <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
                                          fix.N=TRUE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="Palm",
                                          control.opt=control.opt,adj.mvt=TRUE,ft.normal=FALSE))
mletmvt = system.time(testmlemvt <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
                                          fix.N=FALSE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="MLE",
                                          control.opt=control.opt,adj.mvt=TRUE,ft.normal=FALSE,sim.ft.normal=TRUE))
harvestsim(testpalmvt$file)
harvestsim(testmlemvt$file)
msims = readRDS(testmlemvt$file)
psims = readRDS(testpalmvt$file)

pcdiff.palm=100*(psims$mle$D.est-psims$palm$Dhat)/psims$palm$Dhat
quartz(h=4,w=6)
hist(pcdiff.palm,nclass=10,xlab="%LCE bigger than Palm",main="Difference in Density estimates")
obs = 100*(1.24-1.05)/1.05
points(obs,0,pch=19)
cilims = quantile(pcdiff.palm,probs=c(0.025,0.975))
segments(cilims,rep(0,2),cilims,rep(5,2),col="red",lwd=2)


# Debugging:
Nsim = 10

# Simulate from Palm simulator
movement = list(forward=TRUE,sideways=TRUE)
seed = 1
palmtmvt = system.time(palmvt <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
                                       fix.N=TRUE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="Palm",
                                       control.opt=control.opt,adj.mvt=TRUE,ft.normal=TRUE))
harvestsim(palmvt$file)
readRDS(palmvt$file)

# Simulate from LCE simulator, EST with exact forward mvt dbn, SIM with exact forward mvt dbn
seed = 1
tmvt.adj.ft.sft = system.time(mvtadj.ft.sft <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
                                                          fix.N=FALSE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="MLE",
                                                          control.opt=control.opt,adj.mvt=TRUE,ft.normal=FALSE,sim.ft.normal=FALSE))
harvestsim(mvtadj.ft.sft$file)
readRDS(mvtadj.ft.sft$file)

# Simulate from LCE simulator, EST with exact forward mvt dbn, SIM with normal approx forward mvt dbn
seed = 1
tmvtadj.ft.sftnorm = system.time(mvtadj.ft.sftnorm <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
                                                          fix.N=FALSE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="MLE",
                                                          control.opt=control.opt,adj.mvt=TRUE,ft.normal=FALSE,sim.ft.normal=TRUE))
harvestsim(mvtadj.ft.sftnorm$file)
readRDS(mvtadj.ft.sftnorm$file)


# Simulate from LCE simulator, EST with normal approx forward mvt dbn, SIM with normal approx forward mvt dbn
seed = 1
tmvtadj.ftnorm.sftnorm = system.time(mvtadj.ftnorm.sftnorm <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
                                                            fix.N=FALSE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="MLE",
                                                            control.opt=control.opt,adj.mvt=TRUE,ft.normal=TRUE,sim.ft.normal=TRUE))
harvestsim(mvtadj.ftnorm.sftnorm$file)
readRDS(mvtadj.ftnorm.sftnorm$file)



# Try 100 simulations now, with most reliable of above options:
# Simulate from LCE simulator, EST with exact forward mvt dbn, SIM with normal approx forward mvt dbn
Nsim = 100
seed = 1
tmvtadj.ft.sftnorm = system.time(mvtadj.ft.sftnorm <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
                                                            fix.N=FALSE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="MLE",
                                                            control.opt=control.opt,adj.mvt=TRUE,ft.normal=FALSE,sim.ft.normal=TRUE))
harvestsim(mvtadj.ft.sftnorm$file)
mvtadj.ft.sftnorm.sims=readRDS(mvtadj.ft.sftnorm$file)
sim248.Nsim100 = mvtadj.ft.sftnorm.sims
# Save for use in paper:
saveRDS(sim248.Nsim100,file="./inst/results/sim248.Nsim100.Rds")
sim248.Nsim100 = readRDS("./inst/results/sim248.Nsim100.Rds")
sim248.Nsim100.fn = mvtadj.ft.sftnorm$file
saveRDS(sim248.Nsim100.fn,file="./inst/results/sim248.Nsim100.fn.Rds")
pcDdiff = 100*(sim248.Nsim100$mle$D.est/sim248.Nsim100$palm$Dhat-1)
hist(pcDdiff,nclass=20,probability=TRUE,xlab=expression((hat(D)[LCE]-hat(D)[Palm])/hat(D)[Palm]),main="")
CI = quantile(pcDdiff,probs=c(0.025,0.975))
segments(CI,rep(0,2),CI,rep(0.02,2),lwd=2)
pcDdiff.obs = 18
points(pcDdiff.obs,0,pch=19)
2*(101-which(sort(c(pcDdiff,pcDdiff.obs))==pcDdiff.obs))
simres.sim248.Nsim100=harvestsim(sim248.Nsim100.fn)
saveRDS(simres.sim248.Nsim100,file="./inst/results/simres.sim248.Nsim100.Rds")

# Try with 2 times the sample size
L = 2200
Nsim = 50
seed = 1
tmvtadj.ft.sftnorm.2L = system.time(mvtadj.ft.sftnorm.2L <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
                                                            fix.N=FALSE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="MLE",
                                                            control.opt=control.opt,adj.mvt=TRUE,ft.normal=FALSE,sim.ft.normal=TRUE))
harvestsim(mvtadj.ft.sftnorm.2L$file)
# Try with 3 times the sample size
L = 3300
Nsim = 300
seed = 1
tmvtadj.ft.sftnorm.3L = system.time(mvtadj.ft.sftnorm.3L <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
                                                                  fix.N=FALSE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="MLE",
                                                                  control.opt=control.opt,adj.mvt=TRUE,ft.normal=FALSE,sim.ft.normal=TRUE))
harvestsim(mvtadj.ft.sftnorm.3L$file)


# Simulate from LCE simulator, EST with normal approx forward mvt dbn, SIM with normal approx forward mvt dbn
seed = 1
tmvtadj.ftnorm.sftnorm = system.time(mvtadj.ftnorm.sftnorm <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
                                                                    fix.N=FALSE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="MLE",
                                                                    control.opt=control.opt,adj.mvt=TRUE,ft.normal=TRUE,sim.ft.normal=TRUE))
harvestsim(mvtadj.ftnorm.sftnorm$file)
mvtadj.ftnorm.sftnorm.sims = readRDS(mvtadj.ftnorm.sftnorm$file)







# Simulate from LCE simulator, fit with normal approx for forward mvt
seed = 1
mletmvtadj.norm = system.time(mlemvtadj.norm <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
                                          fix.N=TRUE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="MLE",
                                          control.opt=control.opt,adj.mvt=TRUE,ft.normal=TRUE,sim.ft.normal=TRUE))
harvestsim(mlemvtadj.norm$file)
readRDS(mlemvtadj.norm$file)

# Simulate from LCE simulator, fit assuming NO forward mvt
seed = 1
mletnomvt = system.time(mlenomvt <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
                                            fix.N=TRUE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="MLE",
                                            control.opt=control.opt,adj.mvt=FALSE,ft.normal=FALSE))
harvestsim(mlenomvt$file)
readRDS(mlenomvt$file)
