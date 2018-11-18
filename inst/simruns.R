library(twoplane)
library(palm)

#names(sigmas)=c("Hiby","Westgate","mle","palm")

Ec = 100 # mean dive cycle length
alphas = c(10,20,50,80)/100
ks = c(10,20,50,80)
animalspeeds = c(0.65, 0.95, 1.5)/1000 # mean speed in km/sec
# Convert using fact that E(U)=sqrt(2)*gamma(1)/gamma(0.5), if U~Chi(1)
# (See here: https://math.stackexchange.com/questions/1059938/whats-the-expectation-of-square-root-of-chi-square-variable)
sigmarates = animalspeeds/(sqrt(2)*gamma(1)/gamma(0.5))
planeknots=100 # observer speed in knots
nm2km=1.852 # multiplier to convert nautical miles to kilometres
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec


w=0.125*2 # width of strip in km from porpoise data
D = 1.24
En = 100

Nsim=30

alpha = alphas[2]
k = ks[4]
sigmarate = sigmarates[2]
#sigmarate = 0.01/(sqrt(2)*gamma(1)/gamma(0.5))
sigma.mult=5
dmax.km=sigma.mult*sigmarate*sqrt(k) # max dist apart (in km) observations could be considered duplicates; 
dmax.time=(dmax.km/planespd)   # max time apart (in seconds) observations could be considered duplicates
tw=w/planespd # width of searched stip, in plane seconds
btw = tw+2*dmax.time # width of stip with buffer, in plane seconds
tw/btw


seed = 1
# Do one scenaio with a few simulations to check it works:
tm   = system.time(testsim   <- dosim(alpha,Ec,k,w,sigmarate,planespd,D,En=En,sigma.mult=sigma.mult,
                                      seed=seed,Nsim=30,writeout=FALSE,iomvt=FALSE))
tmvt = system.time(testsimvt <- dosim(alpha,Ec,k,w,sigmarate,planespd,D,En=En,sigma.mult=sigma.mult,
                                      seed=seed,Nsim=30,writeout=FALSE,iomvt=TRUE))
harvestsim(alpha,k,sigmarate,D,En=En,Nsim=10,simresults=testsim)
harvestsim(alpha,k,sigmarate,D,En=En,Nsim=10,simresults=testsimvt)

# Then do a bunch
start.a=1; start.k=1; start.s=1
end.a=4; end.k=4; end.s=3
startime=date()
for(na in start.a:end.a) {
  for(nk in start.k:end.k) {
    for(ns in start.s:end.s) {
      dosim(alphas[na],Ec,ks[nk],w,sigmarates[ns],planespd,D,En=En,Nsim,seed=12345)
    }
  }
}
endtime=date()

nscenarios = length(alphas)*length(ks)*length(sigmarates)
NAs = rep(NA,nscenarios)
simtab =  data.frame(Nsim=NAs,alpha=NAs,k=NAs,speed=NAs,D=NAs,
                     pc.bias.mle=NAs,pc.cv.mle=NAs,cover.mle=NAs,
                     pc.bias.palm=NAs,pc.cv.palm=NAs,
                     n1=NAs,n2=NAs,m=NAs,
                     E1=NAs,
                     sigmarate=NAs,
                     sehat.Dhat=NAs,
                     nbadD=NAs,nbadse=NAs)
ns = 0
#for(i in 1:length(alphas)) {
for(i in 3:length(alphas)) {
    for(j in 1:length(ks)) {
    for(l in 1:length(sigmarates)) {
      ns = ns+1
      simtab[ns,] = harvestsim(alphas[i],ks[j],sigmarates[l],D,En,Nsim)
    }
  }
}

saveRDS(simtab,file=paste("./inst/results/simtab_En",En,"_",Nsim,".Rds",sep=""))

make3Dplot = function(simtab,speed,response,zlab=NULL,main="",zlim=NULL,...) {
  responsename = names(simtab)[which(names(simtab)==response)]
  tab3d = simtab[simtab$speed==speed,c("alpha","k",responsename)]
  alpha = unique(tab3d$alpha)
  k = unique(tab3d$k)
  z = matrix(tab3d[,3],nrow=length(x),byrow=TRUE)
  if(is.null(zlab)) zlab = response
  if(is.null(zlim)) zlim = range(z)
  hist3D(alpha,k,z,xlab=expression(alpha),ylab=expression(k),zlab=zlab,main=main,zlim=zlim,clim=zlim,...)
}

library(plot3D)

quartz(h=6,w=12)
par(mfrow=c(2,3))
zlim = range(simtab$pc.bias.mle)
make3Dplot(simtab,0.65,"pc.bias.mle","%Bias","speed=0.65",zlim=zlim)
make3Dplot(simtab,0.95,"pc.bias.mle","%Bias","speed=0.95",zlim=zlim)
make3Dplot(simtab,1.5,"pc.bias.mle","%Bias","speed=1.5",zlim=zlim)
zlim = range(simtab$pc.cv.mle)
make3Dplot(simtab,0.65,"pc.cv.mle","%CV","speed=0.65",zlim=zlim)
make3Dplot(simtab,0.95,"pc.cv.mle","%CV","speed=0.95",zlim=zlim)
make3Dplot(simtab,1.5,"pc.cv.mle","%CV","speed=1.5",zlim=zlim)


quartz(h=12,w=8)
par(mfrow=c(3,2))
zlimBias = range(simtab$pc.bias.mle)
zlimCV = range(simtab$pc.cv.mle)
make3Dplot(simtab,0.65,"pc.bias.mle","%Bias","speed=0.65",zlim=zlimBias)
make3Dplot(simtab,0.65,"pc.cv.mle","%CV","speed=0.65",zlim=zlimCV)
make3Dplot(simtab,0.95,"pc.bias.mle","%Bias","speed=0.95",zlim=zlimBias)
make3Dplot(simtab,0.95,"pc.cv.mle","%CV","speed=0.95",zlim=zlimCV)
make3Dplot(simtab,1.5,"pc.bias.mle","%Bias","speed=1.5",zlim=zlimBias)
make3Dplot(simtab,1.5,"pc.cv.mle","%CV","speed=1.5",zlim=zlimCV)


# Coverage prob

quartz(h=4,w=12)
par(mfrow=c(1,3))
simtab$covererr = simtab$cover.mle-0.95
zlimCover = range(simtab$covererr)
make3Dplot(simtab,0.65,"covererr","Coverage error","speed=0.65",zlim=zlimCover)
make3Dplot(simtab,0.95,"covererr","Coverage error","speed=0.95",zlim=zlimCover)
make3Dplot(simtab,1.5,"covererr","Coverage error","speed=1.5",zlim=zlimCover)



blims = range(simtab$pc.bias.palm,simtab$pc.bias.mle)
cvlims = range(simtab$pc.cv.mle,simtab$pc.cv.palm)
alpha10 = which(simtab$alpha==0.1)
lag80 = which(simtab$k==80)

quartz(h=4,w=8)
par(mfrow=c(1,2))
plot(simtab$pc.bias.mle,simtab$pc.bias.palm,xlab="%Bias MLE",ylab="%Bias Palm",main="",pch=19,cex=0.5)
segments(blims[1],blims[1],blims[2],blims[2])
points(simtab$pc.bias.mle[lag80],simtab$pc.bias.palm[lag80])
points(simtab$pc.bias.mle[alpha10],simtab$pc.bias.palm[alpha10],pch=4)
plot(simtab$pc.cv.mle,simtab$pc.cv.palm,xlab="%CV MLE",ylab="%CV Palm",main="",pch=19,cex=0.5)
segments(cvlims[1],cvlims[1],cvlims[2],cvlims[2])
points(simtab$pc.cv.mle[lag80],simtab$pc.cv.palm[lag80])
points(simtab$pc.cv.mle[alpha10],simtab$pc.cv.palm[alpha10],pch=4)



dosim = function(alpha,Ec,k,w,sigmarate,planespd,D,En=100,Nsim=100,writeout=TRUE,seed=1,iomvt=FALSE,sigma.mult=5) {

  p.up = alpha # proportion of time up
  E1 = alpha*Ec
  p=c(1, 1) # definitely see if available in searched strip
  dmax.t=sigma.mult*(sigmarate*sqrt(k))/planespd # max time apart (in seconds) observations could be considered duplicates
  dmax.d = dmax.t*planespd # max distance apart observations could be considered duplicates
  ps = p.t(E1,Ec,p,sigmarate,k,dmax.t,planespd,w/2,io=iomvt) # capture history probabilities
  p. = sum(ps) # prob detect
  
  b = w
  if(iomvt) b = w + 2*dmax.d
#  L=1100  # length in km from porpoise data
#  L=En/(D*w*alpha) # set L to get desired sample size
  L=En/(D*b*p.) # set L to get desired sample size
  
  p.see.up=c(1,1) # Prob see if up and in
  
  N=D*(L*b)
  
  Dstrip=N/(L*b) # density in number per sq km
  Dstrip.t=D*(planespd^2) # density in planespd units
  D.line.t=Dstrip.t*b/planespd # density in planespd along LINE units (1-dimensional)
  control.opt=list(trace=0,maxit=1000)
  
  estsio=ests=data.frame(Dhat=rep(0,Nsim),E1=rep(0,Nsim),E2=rep(0,Nsim),sigma=rep(0,Nsim),
                         n1=rep(0,Nsim),n2=rep(0,Nsim),m=rep(0,Nsim),mu_c=rep(0,Nsim),
                         se=rep(0,Nsim),inci=rep(0,Nsim))
  estsnspp=estsna=data.frame(Dhat=rep(0,Nsim),n1=rep(0,Nsim),n2=rep(0,Nsim))
  
  ests.kd=ests
  plot.sample=FALSE
  plot.displacement=FALSE
  plot.cuts=FALSE
  checkdists=FALSE
  fromfile=FALSE
  segiotime=segtime=nspptime=rep(NA,Nsim)
  
  estimate=c("D","sigma","E1") # parameters to estimate
  true=list(D=Dstrip,sigma=sigmarate,E1=E1) # parameters to use in simulation
  

  set.seed(seed) # initialise random number sequence (for repeatability)
  skip=c()
  startime=date()
  for(sim in 1:Nsim) {
    sdat=sim.2plane(N,L,w,sigmarate,k,planespd,p.up,Ec,p=p.see.up,
                    sigma.mult=sigma.mult,movement=list(forward=TRUE,sideways=iomvt))
    
    # fit accounting for leakage of animals in and out of strip
    segiotime[sim]=system.time(fitio<-segfit(sdat,D.line.t,E1=E1,Ec=Ec,sigmarate=sigmarate,planespd=planespd,p=c(1,1),sigma.mult=sigma.mult,
                                             control.opt=control.opt,method="BFGS",estimate=estimate,set.parscale=TRUE,
                                             io=iomvt,Dbound=NULL,hessian=TRUE))[3]
    
    # Palm
    pdat = format4palm(sdat,planespd,Ec,sigmarate,sigma.mult=sigma.mult)
    nspptime[sim]=system.time(palmfit<-fit.twocamera(pdat$points,pdat$cameras,pdat$d,pdat$w,pdat$b,pdat$l,pdat$tau,pdat$R,trace=FALSE))[3]
    est.palm=coef(palmfit)
    #  Convert 'activity centre' sigma of Palm into our sigma by *sqrt(2)  and convert to sigmarate.
    est.palm[3]=est.palm[3]*sqrt(2)/sqrt(sdat$k)
    #  
    estsio$n1[sim]=length(sdat$s1)
    estsio$n2[sim]=length(sdat$s2)
    estsnspp$n1[sim]=length(sdat$s1)
    estsnspp$n2[sim]=length(sdat$s2)
    ests$m[sim]=sdat$n11 # record number of actual duplicates
    estsio$m[sim]=sdat$n11 # record number of actual duplicates
    estsnspp$m[sim]=sdat$n11 # record number of actual duplicates
    
    infmat=try(solve(fitio$hessian),silent=TRUE)
    if(!inherits(infmat, "try-error")) {
      intest=logn.seci(log(fitio$D),sqrt(infmat[1,1]))
      estsio$se[sim]=intest$se/(b*planespd)
      estsio$inci[sim]=(intest$lower/(b*planespd)<=Dstrip & Dstrip<=intest$upper/(b*planespd))
    } else skip=c(skip,sim)
    estsio$Dhat[sim]=fitio$D/(b*planespd)
    estsio$E1[sim]=fitio$E[1]
    estsio$E2[sim]=fitio$E[2]
    estsio$sigma[sim]=fitio$sigmarate
    estsio$mu_c[sim]=fitio$mu_c
    
    estsnspp$Dhat[sim]=est.palm[1]
    estsnspp$sigma[sim]=est.palm[3]
    
    #    # Petersen estimator for known recaptures, within strip:
    #    ests.kd$Dhat[sim] = ((sdat$n1+1)*(sdat$n2+1)/(sdat$n11+1) - 1)/(L/planespd)
    
    #    if(sim==1) cat("\nCounter: \n")
    #    if(sim %% 50 == 0) cat(sim,"\n")
    #    else if(sim %% 10 == 0) cat(sim)
    #    else if(sim %% 5 == 0) cat("+")
    #    else cat("-")
  }  # End sim loop
  
  results = list(mle=estsio,palm=estsnspp)
  dir = "./inst/results/"
  fn = paste("sim-alpha_",alpha,"-k_",k,"-sigmarate_",signif(sigmarate,3),"-D_",D,"-En_",En,"-Nsim_",Nsim,".Rds",sep="")
  dirfn = paste(dir,fn,sep="")
  if(writeout) saveRDS(results,file=dirfn)
  else return(list(fn=fn,sim=results))

}  


harvestsim = function(alpha,k,sigmarate,D,En,Nsim,badcut=10,simresults=NULL) {
  
  if(is.null(simresults)) {
    fn = paste("./inst/results/sim-alpha_",alpha,"-k_",k,"-sigmarate_",signif(sigmarate,3),"-D_",D,"-En_",En,"-Nsim_",Nsim,".Rds",sep="")
    sim = readRDS(fn)
  } else {
    fn=simresults$fn
    sim=simresults$sim
  }
  
  alpha = as.numeric(strsplit(strsplit(fn,"alpha_")[[1]][2],"-")[[1]][1])
  k = as.numeric(strsplit(strsplit(fn,"k_")[[1]][2],"-")[[1]][1])
  sigmarate = as.numeric(strsplit(strsplit(fn,"sigmarate_")[[1]][2],"-")[[1]][1])
  D = as.numeric(strsplit(strsplit(fn,"D_")[[1]][2],"-")[[1]][1])
  En = as.numeric(strsplit(strsplit(fn,"En_")[[1]][2],"-")[[1]][1])
  Nsim = as.numeric(strsplit(strsplit(fn,"Nsim_")[[1]][2],".Rds")[[1]][1])
  
  badD = abs(sim[[1]]$Dhat-D)>=badcut | abs(sim[[2]]$Dhat)>=badcut 
  nbadD = sum(badD)
  bad = abs(sim[[1]]$Dhat-D)>=badcut | abs(sim[[2]]$Dhat)>=badcut | is.na(sim[[1]]$se)
  nbad = sum(bad)
  nbadse = nbad - nbadD
  
  pc.bias.mle = 100*(mean(sim[[1]]$Dhat[!bad])-D)/D
  pc.cv.mle = 100*sqrt(var(sim[[1]]$Dhat[!bad]))/mean(sim[[1]]$Dhat[!bad])
  cover.mle = sum(sim[[1]]$inci[!bad])/Nsim
  
  pc.bias.palm = 100*(mean(sim[[2]]$Dhat[!bad])-D)/D
  pc.cv.palm = 100*sqrt(var(sim[[2]]$Dhat[!bad]))/mean(sim[[1]]$Dhat[!bad])
  
  n1 = mean(sim[[1]]$n1[!bad])
  n2 = mean(sim[[1]]$n2[!bad])
  m = mean(sim[[1]]$m[!bad])
  
  E1 = mean(sim[[1]]$E1[!bad])
  
  sigmrate = mean(sim[[1]]$sigma[!bad])

  sehat.Dhat = mean(sim[[1]]$se[!bad])
  
  out = data.frame(Nsim=Nsim,alpha=alpha,k=k,speed=getspeed(sigmarate),D=D,
                   pc.bias.mle=pc.bias.mle,pc.cv.mle=pc.cv.mle,cover.mle=cover.mle,
                   pc.bias.palm=pc.bias.palm,pc.cv.palm=pc.cv.palm,
                   n1=n1,n2=n2,m=m,
                   E1=E1,
                   sigmarate=sigmarate,
                   sehat.Dhat=sehat.Dhat,
                   nbadD=nbadD,nbadse=nbadse)
  
  return(out)
}  

getspeed = function(sigmarate) return(1000*sigmarate*(sqrt(2)*gamma(1)/gamma(0.5)))
