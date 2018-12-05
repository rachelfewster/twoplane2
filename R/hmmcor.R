
hmmcor = function(ts,gamma,tau,sigmarate,planespd,sigma.mult=NULL,b=NULL,io=TRUE,p=c(1,1)) {
  if(is.null(sigma.mult) & is.null(b)) stop("Need to specify b or sigma.mult")
  if(!is.null(b)) {
    if(!is.null(sigma.mult)) {
      warning("b and sigma.mult specified; using b and ignoring sigma.mult.")
      sigma.mult = NULL
    }
    dmax = b - w
  }
  nts = length(ts)
  hmm.cor = rep(NA,nts)
  for(i in 1:nts) {
    if(io) {
      if(!is.null(sigma.mult)) {
        dmax = sigma.mult*sigmarate*sqrt(ts[i])
        b = w + dmax
        dmax.t = dmax/planespd
      } else {
        dmax.t = (b-w)/planespd
      }
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


p.cond.1det = function(ts,gamma,tau,sigmarate,w,sigma.mult=NULL,b=NULL,io=TRUE) {
  if(is.null(sigma.mult) & is.null(b)) stop("Need to specify b or sigma.mult")
  if(!is.null(b)) {
    if(!is.null(sigma.mult) & io) {
      warning("b and sigma.mult specified; using b and ignoring sigma.mult.")
      sigma.mult = NULL
    }
    dmax = b - w
  }
  nts = length(ts)
  cp = rep(NA,nts)
  kappa = gamma * tau
  Qmat=matrix(c(-1/kappa,1/(tau-kappa),1/kappa,-1/(tau-kappa)),nrow=2)
  for(i in 1:nts) {
    if(!is.null(sigma.mult) & io) {
      dmax = sigma.mult*sigmarate*sqrt(ts[i])
      b = w + dmax
    }
    if(io) {
      idbn1. = c(1,0,0,0)
      TPM = make.inout.tpm(sigma=sigmarate*sqrt(ts[i]),dmax=dmax,w=w)
      cp[i] = p.omega.t(ts[i],idbn1.,p1=1,p2=1,Qmat,omega=11,IO=TPM)
    } else {
      idbn1. = c(1,0)
      cp[i] = p.omega.t(ts[i],idbn1.,p1=1,p2=1,Qmat,omega=11,IO=NULL)
      
    }
  }
  return(cp)
}


