library(twoplane)
library(palm)

#names(sigmas)=c("Hiby","Westgate","mle","palm")
nm2km=1.852 # multiplier to convert nautical miles to kilometres
planeknots=100 # observer speed in knots   CHECK THIS
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec

D = D.2D <- 1.05
## Time between cameras (seconds).
k = l <- 20
## Mean dive-cycle duration (seconds).
tau <- 100
## Mean duration of surface phase (seconds).
kappa <- 84
## Animal movement (roughly based on 0.95 m per s).
sigma <- gamma(0.5)*sqrt(((0.95*l/1000)^2)/2)/2
sigmarate = sigma*sqrt(2)/sqrt(k)
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

seed = 123
# Do one scenaio with a few simulations to check it works:
palmtmvt = system.time(testpalmvt <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
                                          fix.N=TRUE,En=NULL,Nsim=3,writeout=TRUE,seed=seed,simethod="Palm",
                                          control.opt=control.opt))
                       
harvestsim(testpalmvt$file)
mletmvt = system.time(testmlemvt <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
                                          fix.N=TRUE,En=NULL,Nsim=3,writeout=TRUE,seed=seed,simethod="MLE",
                                          control.opt=control.opt))
harvestsim(testmlemvt$file)


sigmarates = c(0.65, 0.95, 1.5)/1000
kappas = c(0.2, 0.5, 0.8)*tau
ks = c(10, 20, 50, 80)

# Then do a bunch
fns = c(rep("",length(sigmarates)*length(kappas)*length(ks)))
Nsim = 30
start.a=1; start.k=1; start.s=1
end.a=length(kappas); end.k=length(ks); end.s=length(sigmarates)
simnum = 0
startime=date()
for(na in start.a:end.a) {
  for(nk in start.k:end.k) {
    for(ns in start.s:end.s) {
      sigma = sigmarate/(sqrt(2)/sqrt(ks[nk]))
      b <- w + sigma.mult*sigma
      simnum = simnum+1
      fns[simnum] = dosim(D.2D,L,w,b,sigmarates[ns],ks[nk],planespd,kappas[na],tau,p=p,movement=movement,
            fix.N=TRUE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="MLE",
            control.opt=control.opt)
    }
  }
}
endtime=date()

nscenarios = length(gammas)*length(ks)*length(sigmarates)
NAs = rep(NA,nscenarios)
simtab =  data.frame(Nsim=NAs,gamma=NAs,k=NAs,speed=NAs,D=NAs,
                     pc.bias.mle=NAs,pc.cv.mle=NAs,cover.mle=NAs,
                     pc.bias.palm=NAs,pc.cv.palm=NAs,
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
      simtab[ns,] = harvestsim(fns[ns])
    }
  }
}

saveRDS(simtab,file=paste("./inst/results/simtab_En",En,"_",Nsim,".Rds",sep=""))


require(plot3D)

make3Dplot = function(simtab,speed,response,zlab=NULL,main="",zlim=NULL,decimals=2,...) {
  responsename = names(simtab)[which(names(simtab)==response)]
  tab3d = simtab[signif(simtab$speed,decimals)==signif(speed,decimals),c("gamma","k",responsename)]
  gamma = unique(tab3d$gamma)
  k = unique(tab3d$k)
  z = matrix(tab3d[,3],nrow=length(gamma),byrow=TRUE)
  if(is.null(zlab)) zlab = response
  if(is.null(zlim)) zlim = range(z)
  hist3D(gamma,k,z,xlab=expression(gamma),ylab=expression(k),zlab=zlab,main=main,zlim=zlim,clim=zlim,...)
}

decimals = 2
speeds = unique(signif(simtab$speed,decimals))

pdf(h=9,w=6,file="./inst/figs/biascv.pdf")
par(mfrow=c(3,2))
zlimBias = range(simtab$pc.bias.mle)
zlimCV = range(simtab$pc.cv.mle)
make3Dplot(simtab,speeds[1],"pc.bias.mle","%Bias",paste("speed=",speeds[1],sep=""),zlim=zlimBias)
make3Dplot(simtab,speeds[1],"pc.cv.mle","%CV",paste("speed=",speeds[1],sep=""),zlim=zlimCV)
make3Dplot(simtab,speeds[2],"pc.bias.mle","%Bias",paste("speed=",speeds[2],sep=""),zlim=zlimBias)
make3Dplot(simtab,speeds[2],"pc.cv.mle","%CV",paste("speed=",speeds[2],sep=""),zlim=zlimCV)
make3Dplot(simtab,speeds[3],"pc.bias.mle","%Bias",paste("speed=",speeds[3],sep=""),zlim=zlimBias)
make3Dplot(simtab,speeds[3],"pc.cv.mle","%CV",paste("speed=",speeds[3],sep=""),zlim=zlimCV)
dev.off()


# Coverage prob

quartz(h=4,w=12)
par(mfrow=c(1,3))
simtab$covererr = simtab$cover.mle-0.95
zlimCover = range(simtab$covererr)
make3Dplot(simtab,0.65,"covererr","Coverage error","speed=0.65",zlim=zlimCover)
make3Dplot(simtab,0.95,"covererr","Coverage error","speed=0.95",zlim=zlimCover)
make3Dplot(simtab,1.5,"covererr","Coverage error","speed=1.5",zlim=zlimCover)


scenario = as.character(1:length(simtab$pc.bias.mle))
xylim.bias = c(0,max(abs(simtab$pc.bias.mle),abs(simtab$pc.bias.palm)))
xylim.cv = c(0,max(abs(simtab$pc.cv.mle),abs(simtab$pc.cv.palm)))

pdf(h=4,w=8,file="./inst/figs/mlepalm.pdf")
par(mfrow=c(1,2))
plot(abs(simtab$pc.bias.palm),abs(simtab$pc.bias.mle),xlab="Abs %Bias Palm",ylab="Abs %Bias MLE",main="",
     pch=scenario,col="white",cex=0.5,xlim=xylim.bias,ylim=xylim.bias)
text(abs(simtab$pc.bias.palm),abs(simtab$pc.bias.mle),labels=scenario,cex=0.5)
lines(xylim.bias,xylim.bias)
plot(abs(simtab$pc.cv.palm),abs(simtab$pc.cv.mle),xlab="%CV Palm",ylab="%CV MLE",main="",
     pch=scenario,col="white",cex=0.,xlim=xylim.cv,ylim=xylim.cv)
text(abs(simtab$pc.cv.palm),abs(simtab$pc.cv.mle),labels=scenario,cex=0.5)
lines(xylim.cv,xylim.cv)
dev.off()

simtab[4:6,c("gamma","k","speed","n1","n2","m")]
simtab[7:9,c("gamma","k","speed","n1","n2","m")]


pc.down = 1-pc.up

lim1 = calc.p.01(ks[1],pc.up[1]*tau,tau,p=c(1,1),sigmarate,planespd)
lim2 = calc.p.01(ks[np.k],pc.up[length(pc.up)]*tau,tau,p=c(1,1),sigmarate,planespd)
ylim = range((lim1/pc.down[1]*100)-100,(lim2/pc.down[length(pc.up)]*100)-100)
for(i in 1:np.k) p.k[i]=calc.p.01(ks[i],pc.up[1]*tau,tau,p=c(1,1),sigmarate,planespd)
bias = (p.k/pc.down[1]*100)-100
plot(ks,bias,xlim=c(left,tau),xlab="lag as % of dive cycle length",cex.lab=0.75,
     ylab=expression(100%*%(P(X[2]==1~"|"~X[1]==0)-alpha)/alpha),type="l",col="gray",ylim=ylim)
text(left,bias[1],paste("available" ,round(pc.up[1]*100),"%",sep=""),cex=0.5,pos=4)
for(j in 2:length(pc.up)) {
  for(i in 1:np.k) p.k[i]=calc.p.01(ks[i],pc.up[j]*tau,tau,p=c(1,1),sigmarate,planespd)
  bias = (p.k/pc.down[j]*100)-100
  lines(ks,bias,col="gray") 
  if(j>length(pc.up)-3) text(left,bias[1],paste("available" ,round(pc.up[j]*100),"%",sep=""),cex=0.5,pos=4)
}
for(i in 1:np.k) p.k[i]=calc.p.01(ks[i],mu1,tau,p=c(1,1),sigmarate,planespd)
bias.mle = ((p.k/(tau-mu1)*100)-1)*100
lines(ks,bias.mle)

