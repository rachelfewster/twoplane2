library(twoplane)
library(palm)

nm2km=1.852 # multiplier to convert nautical miles to kilometres
planeknots=100 # observer speed in knots   CHECK THIS
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec

# Density
D = D.2D <- 1.24

# --------------- A few simulations to chck stuff seems to be working -------------------
## Time between cameras (seconds).
k = l <- 20
## Mean dive-cycle duration (seconds).
tau <- 110
## Mean duration of surface phase (seconds).
kappa <- 80
# sigmarate
sigma_palm = 0.15 # estimated Palm-type sigma (in km) from porpoise data, with lag 248 seconds
sigmarate = sigmapalm2sigmarate(sigma_palm,lag=248) # Brownian sigmarate consistent with sigma_palm
sigma = sigmarate*sqrt(k) # Standard error of Brownian movement after k seconds
animalspd = getspeed(sigmarate,248)*1000 # average speed of animals in m/s, after 248 seconds
planespd/(animalspd/1000) # How much faster plane is going than average animal speed
#speed2sigmarate(.95/1000,248) # Bownian motion parameter that gives the observed speed over 248 seconds of .95 m/s
sigma.mult=5 # multiplier of sigma, to set bound for maximum distance moved in k seconds

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

seed = 12345
progbar = TRUE # set to FALSE if you don't want to generate a bar to track progress of simulations
# Uncomment the lines below if you want to do simulations with Palm simulator
#palmtmvt = system.time(testpalmvt <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
#                                          fix.N=TRUE,En=NULL,Nsim=3,writeout=TRUE,seed=seed,simethod="Palm",
#                                          control.opt=control.opt))
#harvestsim(testpalmvt$file)
mletmvt = system.time(testmlemvt <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
                                          fix.N=TRUE,En=NULL,Nsim=3,writeout=TRUE,seed=seed,simethod="MLE",
                                          control.opt=control.opt,adj.mvt=TRUE,ft.normal=FALSE,sim.ft.normal=TRUE,
                                          progbar=progbar))
harvestsim(testmlemvt$file)
# --------------- END of A few simulations to chck stuff seems to be working -------------------



# --------------- Set up and run the full set of simulations -------------------
Nsim = 150 # Number of simulations to do

sigmarates = c(0.5, 0.95, 1.5)/1000
kappas = c(0.2, 0.5, 0.8)*tau
ks = c(10, 20, 50, 80)

fns = c(rep("",length(sigmarates)*length(kappas)*length(ks))) # filenames
start.a=1; start.k=1; start.s=1
end.a=length(kappas); end.k=length(ks); end.s=length(sigmarates)
simnum = 0
#simethod = "Palm"
simethod = "MLE"
Ltype = "FixedL"
fix.N=FALSE  # Allows the abundance to vary betwen simulations (as Poission with rate D.2D*2*b*L)
En = NULL # Let E[n] be determined by L and D.2D.
if(fix.N) Ntype = "FixedN" else Ntype = "RandomN"
if(!is.null(En)) Ltype = "RandomL"

startime=date()
for(na in start.a:end.a) {
  for(nk in start.k:end.k) {
    for(ns in start.s:end.s) {
      sigma = sigmarate/(sqrt(2)/sqrt(ks[nk]))
      b <- w + sigma.mult*sigma
      simnum = simnum+1
      fns[simnum] = dosim(D.2D,L,w,b,sigmarates[ns],ks[nk],planespd,kappas[na],tau,p=p,movement=movement,
                          fix.N=fix.N,En=En,Nsim=Nsim,writeout=TRUE,seed=seed,simethod=simethod,
                          control.opt=control.opt,adj.mvt=TRUE,ft.normal=FALSE,sim.ft.normal=TRUE,
                          progbar=progbar)
    }
  }
}
endtime=date()

nscenarios = length(kappas)*length(ks)*length(sigmarates)
NAs = rep(NA,nscenarios)
simtab =  data.frame(Nsim=NAs,gamma=NAs,k=NAs,speed=NAs,D=NAs,
                     pc.bias.mle=NAs,pc.cv.mle=NAs,cover.mle=NAs,
                     pc.bias.palm=NAs,pc.cv.palm=NAs,
                     Dhat.cor=NAs,
                     n1=NAs,n2=NAs,m=NAs,
                     kappa=NAs,
                     sigmarate=NAs,
                     sehat.Dhat=NAs,
                     nbadD=NAs,nbadse=NAs)
ns = 0
#for(i in 1:length(gammas)) {
for(i in 1:length(kappas)) {
    for(j in 1:length(ks)) {
    for(l in 1:length(sigmarates)) {
      ns = ns+1
      simtab[ns,] = harvestsim(fns[[ns]])
    }
  }
}

simID = paste("-D_",signif(D.2D,3),"-",Ltype,"-",Ntype,"-simethod_",simethod,"-Nsim_",Nsim,sep="")
writeRDS(fns,file=paste("./inst/results/filenames",simID,".Rds",sep=""))
writeRDS(simtab,file=paste("./inst/results/simtab",simID,".Rds",sep=""))
