twoplane.fit=function(sdat,tau,R,all=TRUE){
  pts=matrix(c(sdat$y1,sdat$y2),ncol=1)
  cameras=c(rep(1,length(sdat$y1)),rep(2,length(sdat$y2)))
  d=sdat$L
  w=sdat$w
  b=sdat$b
  l=sdat$k

  est = fit.twocamera(points=pts, cameras=cameras, d=d, w=w, b=b, l=l, tau=tau, R=R)
  
  params=coef(est)
  
  if(all) return(est) else return(params)
}

Palm2mleSimData = function(palmsimdata) {
  ys = palmsimdata$points
  obs = palmsimdata$sibling.list$cameras
  y1 = ys[obs==1]
  y2 = ys[obs==2]
  sdat = list(y1=y1,y2=y2,k=k,L=L,w=w,b=b)
  return(sdat)
}