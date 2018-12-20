Nsim = 150 # Number of simulations done

D.2D = 1.24
simethod = "MLE"
Ltype = "FixedL"
Ntype =  "RandomN"
# Construct simulation set ID from the above, and read the associated filenames and table of results
simID = paste("-D_",signif(D.2D,3),"-",Ltype,"-",Ntype,"-simethod_",simethod,"-Nsim_",Nsim,sep="")
fns = readRDS(paste("./inst/results/filenames",simID,".Rds",sep=""))
simtab = readRDS(paste("./inst/results/simtab",simID,".Rds",sep=""))


boxplotsim(fns,method="mle",stat="D.est")

# Commented out code no longer used in paper
#require(plot3D)

#make3Dplot = function(simtab,speed,response,zlab=NULL,main="",zlim=NULL,decimals=2,...) {
#  responsename = names(simtab)[which(names(simtab)==response)]
#  tab3d = simtab[signif(simtab$speed,decimals)==signif(speed,decimals),c("gamma","k",responsename)]
#  gamma = unique(tab3d$gamma)
#  k = unique(tab3d$k)
#  z = matrix(tab3d[,3],nrow=length(gamma),byrow=TRUE)
#  if(is.null(zlab)) zlab = response
#  if(is.null(zlim)) zlim = range(z)
#  hist3D(gamma,k,z,xlab=expression(gamma),ylab=expression(k),zlab=zlab,main=main,zlim=zlim,clim=zlim,...)
#}

#decimals = 2
#speeds = unique(signif(simtab$speed,decimals))

#pdf(h=9,w=6,file="./inst/figs/biascv.pdf")
#par(mfrow=c(3,2))
#zlimBias = range(simtab$pc.bias.mle)
#zlimCV = range(simtab$pc.cv.mle)
#make3Dplot(simtab,speeds[1],"pc.bias.mle","%Bias",paste("speed=",speeds[1],sep=""),zlim=zlimBias)
#make3Dplot(simtab,speeds[1],"pc.cv.mle","%CV",paste("speed=",speeds[1],sep=""),zlim=zlimCV)
#make3Dplot(simtab,speeds[2],"pc.bias.mle","%Bias",paste("speed=",speeds[2],sep=""),zlim=zlimBias)
#make3Dplot(simtab,speeds[2],"pc.cv.mle","%CV",paste("speed=",speeds[2],sep=""),zlim=zlimCV)
#make3Dplot(simtab,speeds[3],"pc.bias.mle","%Bias",paste("speed=",speeds[3],sep=""),zlim=zlimBias)
#make3Dplot(simtab,speeds[3],"pc.cv.mle","%CV",paste("speed=",speeds[3],sep=""),zlim=zlimCV)
#dev.off()


# Coverage prob

#quartz(h=4,w=12)
#par(mfrow=c(1,3))
#simtab$covererr = simtab$cover.mle-0.95
#zlimCover = range(simtab$covererr)
#make3Dplot(simtab,0.5,"covererr","Coverage error","speed=0.5",zlim=zlimCover)
#make3Dplot(simtab,0.95,"covererr","Coverage error","speed=0.95",zlim=zlimCover)
#make3Dplot(simtab,1.5,"covererr","Coverage error","speed=1.5",zlim=zlimCover)

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

