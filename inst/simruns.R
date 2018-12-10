library(twoplane)
library(palm)

nm2km=1.852 # multiplier to convert nautical miles to kilometres
planeknots=100 # observer speed in knots   CHECK THIS
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec

D = D.2D <- 1.24
## Time between cameras (seconds).
k = l <- 20
## Mean dive-cycle duration (seconds).
tau <- 110
## Mean duration of surface phase (seconds).
kappa <- 80
## Animal movement (roughly based on 0.95 m per s).
#sigma <- gamma(0.5)*sqrt(((0.95*l/1000)^2)/2)/2
#sigmarate = sigma*sqrt(2)/sqrt(k)
sigma_palm = 0.13 # estimated Palm-type sigma (in km) from porpoise data, with lag 248 seconds
# Calcluate Brownian sigmarate consistent with this displacement sigma over 248 seconds:
# Variance of Palm-type movement is half that of LCE-type movement, and Brownian variance is proportional to lag
sigmarate = sigma_palm*2/sqrt(248) # Brownian sigmarate consistent with sigma_palm
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

seed = 12345
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
fns = c(rep("",length(sigmarates)*length(kappas)*length(ks))) # filenames
Nsim = 150
start.a=1; start.k=1; start.s=1
end.a=length(kappas); end.k=length(ks); end.s=length(sigmarates)
simnum = 0
simethod = "Palm"
simethod = "MLE"
Ltype = "RandomL"
fix.N=TRUE
Ntype = "FixedN"
if(!fix.N) Ntype = "RandomN"
En = NULL
if(!is.null(En)) {
  if(!is.null(L)) warning("Ignoring L because En was specified; L has been calculated to give En with given D.2D.")
  L=En/(D.2D*b*p.) # set L to get desired sample size
  Ltype = "RandomL"
}

startime=date()
for(na in start.a:end.a) {
  for(nk in start.k:end.k) {
    for(ns in start.s:end.s) {
      sigma = sigmarate/(sqrt(2)/sqrt(ks[nk]))
      b <- w + sigma.mult*sigma
      simnum = simnum+1
#      fns[simnum] = dosim(D.2D,L,w,b,sigmarates[ns],ks[nk],planespd,kappas[na],tau,p=p,movement=movement,
#                          fix.N=TRUE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="MLE",
#                          control.opt=control.opt)
      fns[simnum] = dosim(D.2D,L,w,b,sigmarates[ns],ks[nk],planespd,kappas[na],tau,p=p,movement=movement,
                          fix.N=fix.N,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod=simethod,
                          control.opt=control.opt)
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
saveRDS(fns,file=paste("./inst/results/filenames",simID,".Rds",sep=""))
saveRDS(simtab,file=paste("./inst/results/simtab",simID,".Rds",sep=""))


boxplotsim(fns,method="mle",stat="D.est")


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

# Bias plot with sign
scenario = as.character(1:length(simtab$pc.bias.mle))
xylim.bias = range(simtab$pc.bias.mle,simtab$pc.bias.palm)
#quartz(h=3,w=9)
pdf(h=3,w=9,file="./inst/figs/mlepalmbias.pdf")
par(mfrow=c(1,3))
plot((simtab$pc.bias.palm),(simtab$pc.bias.mle),xlab="%Bias Palm",ylab="%Bias MLE",main="",
     pch=scenario,col="white",cex=0.5,xlim=xylim.bias,ylim=xylim.bias)
lines(xylim.bias,xylim.bias,col="gray")
text((simtab$pc.bias.palm),(simtab$pc.bias.mle),labels=scenario,cex=0.5)
hist(simtab$pc.bias.mle,nclass=15,xlab="%Bias MLE",main="")
hist(simtab$pc.bias.palm,nclass=15,xlab="%Bias Palm",main="")
dev.off()

# Mean percentage bias:
mean(simtab$pc.bias.mle)
mean(simtab$pc.bias.palm)
# Correlation between estimators
hist(simtab$Dhat.cor)
range(simtab$Dhat.cor)

n1 = simtab$n1
shift = 8
xlim = c(min(n1),max(shift+n1))
plot(n1,simtab$pc.bias.mle,ylim=range(simtab$pc.bias.mle,simtab$pc.bias.palm),
     xlim=xlim,xlab=expression(bar(n)[1]),ylab="% Bias")
points(shift+n1,simtab$pc.bias.palm,pch=3,col="gray")
lines(range(simtab$n1),c(0,0),lty=2)

quartz(h=4,w=8)
par(mfrow=c(1,2))
plot(n1,simtab$pc.cv.mle,ylim=range(simtab$pc.cv.mle,simtab$pc.cv.palm),
     xlim=xlim,xlab=expression(bar(n)[1]),ylab="% CV",main="(a)")
points(shift+n1,simtab$pc.cv.palm,pch=3,col="gray")
lines(range(simtab$n1),c(0,0),lty=2)

cvdiff = 100*(simtab$pc.cv.palm - simtab$pc.cv.mle)/simtab$pc.cv.mle
plot(n1,cvdiff,ylim=range(cvdiff),
     xlim=xlim,xlab=expression(bar(n)[1]),ylab="relative % difference in CV",main="(n)")
lines(range(simtab$n1),c(0,0),lty=2)

# Abs bias and CV plot
scenario = as.character(1:length(simtab$pc.bias.mle))

#badone = which(simtab$pc.bias.mle>10)
#xylim.abias = c(0,max(abs(simtab$pc.bias.mle[-badone]),abs(simtab$pc.bias.palm[-badone])))
#xylim.cv = c(0,max(abs(simtab$pc.cv.mle[-badone]),abs(simtab$pc.cv.palm[-badone])))

xylim.abias = c(0,max(abs(simtab$pc.bias.mle),abs(simtab$pc.bias.palm)))
xylim.cv = c(0,max(abs(simtab$pc.cv.mle),abs(simtab$pc.cv.palm)))

#quartz(h=4,w=8)
pdf(h=4,w=8,file="./inst/figs/mlepalm.pdf")
par(mfrow=c(1,2))
plot(abs(simtab$pc.bias.palm),abs(simtab$pc.bias.mle),xlab="Abs %Bias Palm",ylab="Abs %Bias MLE",main="",
     pch=scenario,col="white",cex=0.5,xlim=xylim.abias,ylim=xylim.abias)
lines(xylim.abias,xylim.abias,col="gray")
text(abs(simtab$pc.bias.palm),abs(simtab$pc.bias.mle),labels=scenario,cex=0.5)
plot(abs(simtab$pc.cv.palm),abs(simtab$pc.cv.mle),xlab="%CV Palm",ylab="%CV MLE",main="",
     pch=scenario,col="white",cex=0.,xlim=xylim.cv,ylim=xylim.cv)
lines(xylim.cv,xylim.cv,col="gray")
text(abs(simtab$pc.cv.palm),abs(simtab$pc.cv.mle),labels=scenario,cex=0.5)
dev.off()

quartz(h=3,w=9)
par(mfrow=c(1,3))
pccvdiff = 100*(simtab$pc.cv.palm-simtab$pc.cv.mle)/simtab$pc.cv.mle
plot(simtab$pc.cv.mle,pccvdiff)
lines(c(0,max(simtab$pc.cv.mle)),c(0,0),lty=2)
plot(simtab$n1,pccvdiff)
difflm = lm(pccvdiff~n1,data.frame(pccvdiff=pccvdiff,n1=simtab$n1))
predf = data.frame(n1=seq(min(simtab$n1),max(simtab$n1),length=100))
diffy = predict(difflm,type="response",newdata=predf)
lines(predf$n1,diffy)
coefficients(difflm)
lines(c(0,max(simtab$pc.cv.mle)),c(0,0),lty=2)
hist(pccvdiff)
mean(pccvdiff)







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

