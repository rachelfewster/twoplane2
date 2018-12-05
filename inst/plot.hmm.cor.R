rm(list=ls())

w=0.125;sigma.mult=5;tau=110; k=20

mle <- readRDS("./inst/results/survey.mle.Rds") 
sigmarate=mle$sigmarate["est"]
gama=mle$gamma["est"]
nm2km=1.852 # multiplier to convert nautical miles to kilometres
planeknots=100 # observer speed in knots   CHECK THIS
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec

nts = 300
ts = seq(0.01,480,length=nts) 
b = w + sigma.mult*sigmarate*sqrt(max(ts))
gamas = c(round(mle$gamma["est"],2),(2:9)/10)
ngamas = length(gamas)
hmm.cor = hmm.cor0 = p.c.1 = p.c.1nomvt = matrix(rep(NA,ngamas*nts),nrow=ngamas)
for(i in 1:ngamas) {
  #  pcbias[i,] = calc.pcbias(ts,gamas[i],tau,sigmarate,w,sigma.mult)
  p.c.1[i,] = p.cond.1det(ts,gamas[i],tau,sigmarate,w,b=b)
  p.c.1nomvt[i,] = p.cond.1det(ts,gamas[i],tau,sigmarate,w,b=b,io=FALSE)
  #  hmm.cor[i,] = hmmcor(ts,gamas[i],tau,sigmarate,planespd,sigma.mult=sigma.mult,io=TRUE,p=c(1,1))
  hmm.cor[i,] = hmmcor(ts,gamas[i],tau,sigmarate,planespd,b=b,io=TRUE,p=c(1,1))
  #  hmm.cor0[i,] = hmmcor(ts,gamas[i],tau,sigmarate,planespd,sigma.mult=sigma.mult,io=FALSE,p=c(1,1))
  hmm.cor0[i,] = hmmcor(ts,gamas[i],tau,sigmarate,planespd,b=b,io=FALSE,p=c(1,1))
}
# Here to plot correlation
for(i in 1:ngamas) {
  if(i==1) {
    plot(ts,hmm.cor[1,],ylim=c(0,1),type="l",ylab="Correlatoin",xlab="Lag as % of dive cycle length",lwd=2)
    for(j in 1:ngamas) lines(ts,hmm.cor0[j,],lty=2,col="gray")
  }
  else lines(ts,hmm.cor[i,],lty=1,col="black")
}
lines(ts,hmm.cor[1,],lty=1,col="black",lwd=2)

#lines(range(ts),rep(0,2),lty=2)
## Here for case where animals do not move
#for(i in 1:ngamas) {
#  if(i==1) plot(ts,hmm.cor0[1,],ylim=c(0,1),type="l",ylab="Correlatoin",xlab="Lag as % of dive cycle length",lwd=2)
#  else lines(ts,hmm.cor0[i,],lty=i)
#  lines(range(ts),c(0,0))
#}
# Here to plot conditional detection prob
ratios = p.c.1*0
for(i in 1:ngamas) ratios[i,] = p.c.1[i,]/(gamas[i]*w/b)
for(i in 1:ngamas) ratios[i,] = p.c.1[i,]/(gamas[i])
for(i in 1:ngamas) ratios[i,] = p.c.1[i,]/p.c.1nomvt[i,]
for(i in 1:ngamas) {
  if(i==1) plot(ts,ratios[1,],ylim=range(1,ratios),xlim=c(0,100),type="l",ylab=expression(PRatio(X[2]==1~"|"~X[1]==1)),xlab="Lag (seconds)",lwd=2)
  else lines(ts,ratios[i,],lty=1)
}
lines(range(ts),rep(1,2),lty=2)

for(i in 1:ngamas) {
  if(i==1) plot(ts,p.c.1[i,],ylim=range(1,p.c.1,p.c.1nomvt),
                type="l",ylab=expression(P(X[2]==1~"|"~X[1]==1)),xlab="Lag (seconds)",lwd=2)
  else lines(ts,p.c.1[i,],lty=1)
}
for(i in 1:ngamas) {
  lines(ts,p.c.1nomvt[i,],lty=1,col="red")
}
