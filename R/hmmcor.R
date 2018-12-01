rm(list=ls())

w=0.125;sigma.mult=8;tau=100; k=20

mle <- readRDS("./inst/results/survey.mle.Rds") 
sigmarate=mle$sigmarate["est"]
gama=mle$gamma["est"]
nm2km=1.852 # multiplier to convert nautical miles to kilometres
planeknots=100 # observer speed in knots   CHECK THIS
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec

nts = 300
ts = seq(0.01,100,length=nts) 
idbn = c(1,0)
gamas = c(round(mle$gamma["est"],2),(1:9)/10)
ngamas = length(gamas)
hmm.cor = hmm.cor0 = pcbias = matrix(rep(NA,ngamas*nts),nrow=ngamas)
for(i in 1:ngamas) {
  kappa = gamas[i] * tau
  Qmat=matrix(c(-1/kappa,1/(tau-kappa),1/kappa,-1/(tau-kappa)),nrow=2)
  pcbias[i,] = 100*(p.omega.t(ts,idbn,p1=1,p2=1,Qmat,omega=11,IO=NULL)-gamas[i])/gamas[i]
  hmm.cor[i,] = hmmcor(ts,gamas[i],tau,sigmarate,planespd,sigma.mult,io=TRUE,p=c(1,1))
  hmm.cor0[i,] = hmmcor(ts,gamas[i],tau,sigmarate,planespd,sigma.mult,io=FALSE,p=c(1,1))
}
for(i in 1:ngamas) {
  if(i==1) plot(ts,hmm.cor[1,],ylim=c(0,1),type="l",ylab="Correlatoin",xlab="Lag as % of dive cycle length",lwd=2)
  else lines(ts,hmm.cor[i,],ylim=c(0,max(pcbias)),lty=i)
  lines(range(ts),c(0,0))
}
for(i in 1:ngamas) {
  if(i==1) plot(ts,hmm.cor0[1,],ylim=c(0,1),type="l",ylab="Correlatoin",xlab="Lag as % of dive cycle length",lwd=2)
  else lines(ts,hmm.cor0[i,],ylim=c(0,max(pcbias)),lty=i)
  lines(range(ts),c(0,0))
}
for(i in 1:ngamas) {
  if(i==1) plot(ts,pcbias[1,],ylim=c(0,max(pcbias)),type="l",ylab="%bias",xlab="Lag as % of dive cycle length",lwd=2)
  else lines(ts,pcbias[i,],ylim=c(0,max(pcbias)),lty=i)
}


E1=kappa
Ec=tau
p=c(1,1)
k=20
sigmarate=mle$sigmarate["est"]
gama=mle$gamma["est"]
nm2km=1.852 # multiplier to convert nautical miles to kilometres
planeknots=100 # observer speed in knots   CHECK THIS
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec
sigma.mult = 8

hmmcor = function(ts,gamma,tau,sigmarate,planespd,sigma.mult=5,io=TRUE,p=c(1,1)) {
  nts = length(ts)
  hmm.cor = rep(NA,nts)
  for(i in 1:nts) {
    dmax.t = sigma.mult*sigmarate*sqrt(ts[i])/planespd
    if(io) {
      b = w + dmax.t*planespd
      idbn4 = c(gamma*w/b, gamma*(1-w/b), (1-gamma)*w/b, (1-gamma)*(1-w/b))
      p4 = p.t(gamma*tau, tau, p, sigmarate, ts[i], dmax.t, planespd, halfw.dist=w, io=TRUE, idbn=idbn4)
      pnames = names(p4)
      p4 = c(as.numeric(p4), 1-sum(as.numeric(p4)))
      names(p4) = c(pnames,"ch00")
      px4 = gamma*w/b
      cov4 = px4^2*p4["ch00"] - px4*(1-px4)*p4["ch10"] - (1-px4)*px4*p4["ch01"] + (1-px4)^2*p4["ch11"]
      vx = px4*(1-px4)
      hmm.cor[i] = cov4/vx
    } else {
      b = w
      idbn2 = c(gamma, (1-gamma))
      p2 = p.t(gamma*tau, tau, p, sigmarate, ts[i], dmax.t, planespd, halfw.dist=w, io=FALSE, idbn=idbn2)
      pnames = names(p2)
      p2 = c(as.numeric(p2), 1-sum(as.numeric(p2)))
      names(p2) = c(pnames,"ch00")
      px2 = gamma
      cov2 = px2^2*p2["ch00"] - px2*(1-px2)*p2["ch10"] - (1-px2)*px2*p2["ch01"] + (1-px2)^2*p2["ch11"]
      vx = px2*(1-px2)
      hmm.cor[i] = cov2/vx
    }
  }
  return(hmm.cor)
}



hmmcor(k,gama,tau,sigmarate,planespd,sigma.mult=5,io=FALSE,p=c(1,1))
hmmcor(k,gama,tau,sigmarate,planespd,sigma.mult=5,io=TRUE,p=c(1,1))
hmmcor(k,gama,tau,0,planespd,sigma.mult=5,io=TRUE,p=c(1,1))


hmmcor(ts,gama,tau,sigmarate,planespd,sigma.mult=5,io=FALSE,p=c(1,1))


