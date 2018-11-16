format4palm=function(sdat,planespd,Ec,sigma,sigma.mult=5){
  pts=c(sdat$s1,sdat$s2)*planespd
  cameras=c(rep(1,length(sdat$s1)),rep(2,length(sdat$s2)))
  d=sdat$tL*planespd
  w=0.5*sdat$tw*planespd # sdat$tw is width in plane-seconds; w must be half-width in distance
  b=w+sigma.mult*sigma*sqrt(sdat$k)
  l=sdat$k
  tau=Ec
  R=2*b
  
  return(list(points=pts,cameras=cameras,d=d,w=w,b=b,l=l,tau=tau,R=R))
}


twoplane.fit=function(sdat,planespd,Ec,sigma,sigma.mult=6,trace=FALSE,all=FALSE){
  pts=c(sdat$s1,sdat$s2)*planespd
  cameras=c(rep(1,length(sdat$s1)),rep(2,length(sdat$s2)))
  
  l=sdat$tL*planespd
  w=0.5*sdat$tw*planespd # sdat$tw is width in plane-seconds; w must be half-width in distance
  b=w+sigma.mult*sigma*sqrt(sdat$k)
  t=sdat$k
  C=Ec
  R=2*b
  
  est=fit.twocamera(pts,cameras,l,w,b,t,C,R,trace=trace)
  
  params=coef(est)
  
  #  Convert 'activity centre' sigma of Palm into our sigma by *sqrt(2)  and convert to sigmarate.
  params[3]=params[3]*sqrt(2)/sqrt(sdat$k)
  
  #cat(paste("Estimate from Palm:", params[1]," real value:", Dstrip))
  if(all) return(est) else return(params)
}
