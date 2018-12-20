# Dummy function while working on forward movement correction code:
f.d=function(delta,k,theta3,dmax) return(1)

#' @title PDF of Brownian hitting time
#'
#' @description
#'  Calculates the PDF of the time to hit a barrier for a particle moving with 1D 
#'  Brownian motion with movement rate paramter \code{exp(theta3)} and drift 
#'  speed \code{spd} towards the barrier, when the barrier is initially at a
#'  distance \code{k*spd} away.
#'  
#' @param tt The hitting time
#' @param k The hitting time if the particle did not move.
#' @param spd The speed at which the barrier drifts.
#' @param theta3 The log of the Brownian motion rate parameter.
#' 
#' @examples 
dft = function(ft,k,spd,theta3) { 
  sigmarate=exp(theta3)
  return(spd*k*exp(-spd^2*(k-ft)^2/(2*sigmarate^2*ft))/sqrt(2*pi*sigmarate^2*ft^3))
}

#' @title PDF of Brownian hitting distance
#'
#' @description
#'  Calculates the PDF of the distance a particle has moved when it hits a 
#'  barrier, for a particle moving with 1D Brownian motion with movement rate 
#'  paramter \code{exp(theta3)} and drift speed \code{spd} towards the barrier, 
#'  when the barrier is initially at a distance \code{k*spd} away.
#'  
#' @param fd The hitting time
#' @param k The hitting time if the particle did not move.
#' @param spd The speed at which the barrier drifts.
#' @param theta3 The log of the Brownian motion rate parameter.
#' 
#' @examples 
dfd = function(fd,k,spd,theta3) dft(fd/spd,k,spd,theta3)*spd

#' @title Random variable for Brownian hitting time
#'
#' @description
#'  Generaes a random variable for the time to hit a barrier for a particle moving with 1D 
#'  Brownian motion with movement rate paramter \code{exp(theta3)} and drift 
#'  speed \code{spd} towards the barrier, when the barrier is initially at a
#'  distance \code{k*spd} away.
#'  
#' @param ft The hitting time
#' @param k The hitting time if the particle did not move.
#' @param spd The speed at which the barrier drifts.
#' @param theta3 The log of the Brownian motion rate parameter.
#' 
#' @examples 
rft = function(n,k,spd,theta3,prop.mult=1.5) {
  sigmarate=exp(theta3)
  sigma.t = sqrt(k)*sigmarate/planespd
  done=FALSE
  t1 = NULL
  while(!done) {
    t0 = rnorm(10*n,mean=k,sd=sigma.t*prop.mult) # proposal dbn (normal with bigger sigma)
    t0 = t0[t0>0] # chuck the negative times
    t1 = c(t1,t0[runif(length(t0))<=dft(t0,k,planespd,theta3)]) # acceptance step
    if(length(t1)>=n) done = TRUE
  }
  ft = t1[1:n]
  return(ft)
}




#' @title Calculates observable capture histories
#'
#' @description
#'  Calculates an observable capture history \code{omega}, given an initial state
#'  vector,for animals available according to Markov process and allowed to move 
#'  in and out of searched strip according to Brownian motion. 
#'  
#' @param t The time(s) bewtween 2 observers passing
#' @param idbn The initial state distribution.
#' If no movement assumed (\code{IO} is NULL), this is  (p(up), p(down));
#' If movement is assumed, it is  (p(in & up), p(out & up), p(in & down), p(out & down)).
#' @param p1 Probability that obs 1 detects animal that is available and in searched strip
#' @param p2 Probability that obs 2 detects animal that is available and in searched strip
#' @param Qmat Transition intensity matrix for up-down Markov process
#' @param IO The in-out transition probability matrix. Must either be 2x2 or 2x2xlength(t), 
#' with IO[,,i] being the in-out transition probability matrix for time t[i].
#' 
#' @examples 
p.omega.t=function(t,idbn,p1,p2,Qmat,omega,IO=NULL){
  nt=length(t)
  p=rep(NA,nt)
  if(is.null(IO)) {
    if(omega==10){
      P1=diag(c(p1,0))
      P2c=diag(c(1-p2,1))
      for(i in 1:nt) p[i]=idbn%*%P1%*%expm(t[i]*Qmat)%*%P2c%*%c(1,1)
    }else if(omega==01){
      P1c=diag(c(1-p1,1))
      P2=diag(c(p2,0))
      for(i in 1:nt) p[i]=idbn%*%P1c%*%expm(t[i]*Qmat)%*%P2%*%c(1,1)    
    }else if(omega==11){
      P1=diag(c(p1,0))
      P2=diag(c(p2,0))
      for(i in 1:nt) p[i]=idbn%*%P1%*%expm(t[i]*Qmat)%*%P2%*%c(1,1)
    }else stop(paste("Invalid omega:",omega))
  } else { # Here if you have an IO transition pobability matrix
    d = c(omega%/%10,omega%%10)
    P1 = diag(c(p1^d[1]*(1-p1)^(1-d[1]),rep((1-d[1]),3)))
    P2 = diag(c(p2^d[2]*(1-p2)^(1-d[2]),rep((1-d[2]),3)))
    if(length(dim(IO))==2) { # here if have only one TPM matrix
      # warning("You have multiple times but only one IO transition probability matrix.")
      for(i in 1:nt) {
        TPM = kronecker(expm(t[i]*Qmat),IO)
        p[i]=idbn%*%P1%*%TPM%*%P2%*%rep(1,4)
      }
    } else {
      if(length(dim(IO))==3) {
        for(i in 1:nt) {
          TPM = kronecker(expm(t[i]*Qmat),IO[,,i]) # Use ith IO transition probability matrix
          p[i]=idbn%*%P1%*%TPM%*%P2%*%rep(1,4)
        }
        } else {
        stop("Your IO transition probability matrix has to be 2D or 3D (wih appropriate dimension lengths).")
      }
    }
  }
  return(round(p,10))
}



#' @title Calculates the probabilitis of the observable capture histories
#'
#' @description
#'  Calculates the probabilitis of the three observable capture histories for animals available 
#'  according to Markov process and allowed to move in and out of searched
#'  strip according to Brownian motion. 
#'  
#' @param E1 expected time available per dive cycle
#' @param Ec expected dive cycle length
#' @param p vector of probs that each obs detects available animal in searched strip
#' @param sigmarate Brownian movement rate parameter
#' @param k time lag between observers
#' @param dmax.t maximum distance can move (in plane-seconds)
#' @param planespd speed observers move along transect
#' @param halfw.dist half-width of searched strip
#' @param adj.mvt if TRUE, uses Brownian hitting time PDF to average over in-out TPMs.
#' @param nts number of integration points for averaging over in-out TPMs.
#' @param ft.normal if TRUE uses normal to approximate Brownian hitting times, else 
#' uses exact expression for Brownian hitting times.
#' 
#' @examples 
p.t = function(E1,Ec,p,sigmarate,k,dmax.t,planespd,halfw.dist=NULL,adj.mvt=FALSE,io=TRUE,
               idbn=NULL,nts=200,ft.normal=TRUE) {
  Qmat=matrix(c(-1/E1,1/(Ec-E1),1/E1,-1/(Ec-E1)),nrow=2)
  # deal with in-out movement:
#  TPM = make.inout.tpm(sigma=sigmarate*sqrt(k)*planespd,dmax=dmax.t*planespd,w=halfw.dist) # in-out transition probability matrix
  if(io) {
    if(is.null(halfw.dist)) stop("Need halfw.dist if io=TRUE.")
    if(adj.mvt) {
      ts = seq(k-dmax.t,k+dmax.t,length=nts) # integration points
      dts = diff(ts)
      dts = c(dts[1]/2,dts[2:(nts-1)],dts[nts-1]/2)
      TPM = make.inout.tpm.array(sigma=sigmarate*sqrt(ts),dmax=dmax.t*planespd,w=halfw.dist) # array with multiple in-out transition probability matrices
    } else {
      TPM = make.inout.tpm(sigma=sigmarate*sqrt(k),dmax=dmax.t*planespd,w=halfw.dist) # in-out transition probability matrix
    }
    p.in = halfw.dist/(halfw.dist+dmax.t*planespd) # unconditional probability is in strip, given within dmax.t*planespd of strip
    if(!is.null(idbn)) {
      if(length(idbn)!=4) stop ("idbn must be of length 4")
    } else {
      idbn = c((E1/Ec)*p.in,(E1/Ec)*(1-p.in),(1-E1/Ec)*p.in,(1-E1/Ec)*(1-p.in))
    }
  } else {
    TPM = NULL
    if(!is.null(idbn)) {
      if(length(idbn)!=2) stop ("idbn must be of length 2")
    } else {
      idbn = c((E1/Ec),(1-(E1/Ec)))
    }
  }
  if(adj.mvt & io) {
    if(ft.normal) {
      p01.k=sum(p.omega.t(ts,idbn,p[1],p[2],Qmat,omega=01,TPM)*dnorm(ts,mean=k,sd=sqrt(k)*sigmarate/planespd)*dts)
      p10.k=sum(p.omega.t(ts,idbn,p[1],p[2],Qmat,omega=10,TPM)*dnorm(ts,mean=k,sd=sqrt(k)*sigmarate/planespd)*dts)
      p11.k=sum(p.omega.t(ts,idbn,p[1],p[2],Qmat,omega=11,TPM)*dnorm(ts,mean=k,sd=sqrt(k)*sigmarate/planespd)*dts)
    } else {
      p01.k=sum(p.omega.t(ts,idbn,p[1],p[2],Qmat,omega=01,TPM)*dft(ts,k,planespd,log(sigmarate))*dts)
      p10.k=sum(p.omega.t(ts,idbn,p[1],p[2],Qmat,omega=10,TPM)*dft(ts,k,planespd,log(sigmarate))*dts)
      p11.k=sum(p.omega.t(ts,idbn,p[1],p[2],Qmat,omega=11,TPM)*dft(ts,k,planespd,log(sigmarate))*dts)
    }
  } else {
    p01.k=p.omega.t(k,idbn,p[1],p[2],Qmat,omega=01,TPM)
    p10.k=p.omega.t(k,idbn,p[1],p[2],Qmat,omega=10,TPM)
    p11.k=p.omega.t(k,idbn,p[1],p[2],Qmat,omega=11,TPM)
  }
  return(data.frame(ch01=p01.k,ch10=p10.k,ch11=p11.k))
}



#' @title Prob goes from outside strip to inside strip
#'
#' @description
#'  Returns the transition probability from outside strip to inside strip. 
#'  Starting coordinate outside strip is assumed to be a uniform random variable 
#'  on \code{(w,w+dmax)}, where \code{w} is half-width of strip, \code{dmax} is maximum 
#'  distance from strip to consider, and cooridnate zero is centre of strip. Position 
#'  after movement is assumed to be normal with mean equal to starting coordinate 
#'  and standared deviation \code{sigma}. 

#'  
#'  Deals only with positive initial coordinate as problem is symmetric about zero.
#'  Integrates CDF of normal distribution only inside strip, from initial coordinate 
#'  \code{w} to \code{w+dmax}. Integration is done by trapezoid rule. 
#'  
#' @param sigma standard deviation of nomal location distribution
#' @param dmax maximum distance can move
#' @param w half-width of strip
#' @param nx number of intervals to use for approximate integral.
#' 
#' @examples 
#' sigma = 3
#' dmax = 9*sigma
#' nx=30
#' dx=dmax/nx
#' w = 1
#' xs=seq(w,dmax+w,dx)
#' ws = rep(w,length(xs))
#' F.rs = F.r.norm(ws,sigma,xtrunc=dmax,mean=xs) - F.r.norm(-ws,sigma,xtrunc=dmax,mean=xs)
#' plot(xs,F.rs,type="l")
#' p.o2i(sigma,dmax=sigma,w,nx=30)
#p.o2i = function(sigma,dmax,w,nx=500) {
#  dx=dmax/nx
#  xs=seq(w,dmax+w,dx)
#  ws = rep(w,length(xs))
#  xs = xs + dx/100000 # so that all starting points are outside strip
#  F.rs = F.r.norm(ws,sigma,xtrunc=dmax,mean=xs) - F.r.norm(-ws,sigma,xtrunc=dmax,mean=xs)
#  p = (F.rs[1]/2 + sum(F.rs[2:nx]) + F.rs[nx+1]/2)*dx/dmax
#  return(p)
#}

p.o2i = function(sigma,dmax,w,nx=500) {
  if(sigma==0) sigma = 1e-15
  dx=dmax/nx
  xs=seq(w,dmax+w,length=(nx+1))
  ws = rep(w,length(xs))
#  xs = xs + dx/100000 # so that all starting points are outside strip
  F.rs = pnorm(ws-xs,mean=0,sd=sigma) - pnorm(-ws-xs,mean=0,sd=sigma)
  p = (F.rs[1]/2 + sum(F.rs[2:nx]) + F.rs[nx+1]/2)*dx/(dmax)
  return(p)
}


#' @title Prob goes from inside strip to outside strip
#'
#' @description
#'  Returns the transition probability from inside strip to outside strip. 
#'  Starting coordinate outside strip is assumed to be a uniform random variable 
#'  on \code{(-w,w)}, where \code{w} is half-width of strip and cooridnate zero 
#'  is centre of strip. Position after movement is assumed to be normal with 
#'  mean equal to starting coordinate and standared deviation \code{sigma}. 
#'  
#'  Deals only with positive initial coordinate as problem is symmetric about zero. 
#'  Integrates CDF of normal distribution only inside strip, from initial coordinate 
#'  \code{0} to \code{w}. Integration is done by trapezoid rule. 
#'  
#' @param sigma standard deviation of nomal location distribution
#' @param dmax maximum distance can move
#' @param w half-width of strip
#' @param nx number of intervals to use for approximate integral.
#' 
#' @examples 
#' w = 1
#' sigma = 0.1
#' dmax = 5*sigma
#' p.i2o(sigma,dmax=dmax,w)
#p.i2o = function(sigma,dmax,w,nx=500) return(1-p.i2i(sigma,dmax,w,nx))
p.i2o = function(sigma,dmax,w,nx=500) {
  if(sigma==0) sigma = 1e-15
  xs=seq(0,w,length=(nx+1))
  dx = w/nx
  ws = rep(w,length(xs))
  F.rs = pnorm(xs-ws,mean=0,sd=sigma) + pnorm(-xs-ws,mean=0,sd=sigma)
  p = (F.rs[1]/2 + sum(F.rs[2:nx]) + F.rs[nx+1]/2)*dx/w
  return(p)
}

#' @title Prob goes from inside strip to inside strip
#'
#' @description
#'  Returns the transition probability from inside strip to inside strip. 
#'  Starting coordinate inside strip is assumed to be a uniform random variable 
#'  on \code{(-w,w)}, where \code{w} is half-width of strip and cooridnate zero 
#'  is centre of strip. Position after movement is assumed to be normal with 
#'  mean equal to starting coordinate and standared 
#'  deviation \code{sigma}. 
#'  
#'  Deals only with positive initial coordinate as problem is symmetric about zero. 
#'  Integrates CDF of normal distribution only inside strip, from initial coordinate 
#'  \code{0} to \code{w}. Integration is done by trapezoid rule. 
#'  
#' @param sigma standard deviation of nomal location distribution
#' @param dmax maximum distance can move
#' @param w half-width of strip
#' @param nx number of intervals to use for approximate integral.
#' 
#' @examples 
#' w = 1
#' sigma = 0.1
#' dmax = 5*sigma
#' nx=500
#' dx=w/nx
#' xs=seq(0,w,dx)
#' ws = rep(w,length(xs))
#' F.rs = F.r.norm(ws,sigma,xtrunc=dmax,mean=xs) - F.r.norm(-ws,sigma,xtrunc=dmax,mean=xs)
#' plot(xs,F.rs,type="l") # for debugging - remove when confident is debugged
#' p.i2i(sigma,dmax=dmax,w)
#' 
#p.i2i = function(sigma,dmax,w,nx=500) {
#  dx=w/nx
#  xs=seq(0,w,dx)
#  xs = xs - dx/100000 # so that all starting points are inside strip
#  ws = rep(w,length(xs))
#  F.rs = F.r.norm(ws,sigma,xtrunc=dmax,mean=xs) - F.r.norm(-ws,sigma,xtrunc=dmax,mean=xs)
#  p = (F.rs[1]/2 + sum(F.rs[2:nx]) + F.rs[nx+1]/2)*dx/w
#  return(p)
#}
p.i2i = function(sigma,dmax,w,nx=500) return(1-p.i2o(sigma,dmax,w,nx))


#' @title Prob goes from outside strip to outside strip
#'
#' @description
#'  Returns the transition probability from outside strip to outside strip. 
#'  Starting coordinate outside strip is assumed to be a uniform random variable 
#'  on \code{(w,w+dmax)}, where \code{w} is half-width of strip, \code{dmax} is maximum 
#'  distance from strip to consider, and cooridnate zero is centre of strip. Position 
#'  after movement is assumed to be normal with mean equal to starting coordinate 
#'  and standared deviation \code{sigma}. 

#'  
#'  Deals only with positive initial coordinate as problem is symmetric about zero.
#'  Integrates CDF of normal distribution only outside strip, from initial coordinate 
#'  \code{w} to \code{w+dmax}. Integration is done by trapezoid rule. 
#'  
#' @param sigma standard deviation of nomal location distribution
#' @param dmax maximum distance can move
#' @param w half-width of strip
#' @param nx number of intervals to use for approximate integral.
#' 
#' @examples 
#' w = 1
#' sigma = 0.1
#' dmax = 5*sigma
#' p.o2o(sigma,dmax=dmax,w)
p.o2o = function(sigma,dmax,w,nx=500) return(1-p.o2i(sigma,dmax,w,nx))


#' @title Make in-out transition probability matrix
#'
#' @description
#'  Makes the transition probability matrix for in-out movement from strip. Position 
#'  after movement is assumed to be normal with mean equal to #'  starting 
#'  coordinate and standared deviation \code{sigma}. 
#'  
#' @param sigma standard deviation of nomal location distribution
#' @param dmax maximum distance can move
#' @param w half-width of strip
#' @param nx number of intervals to use for approximate integral.
#' 
#' @examples 
#' w = 1
#' sigma = 0.1
#' dmax = 5*sigma
#' tpm = make.inout.tpm(sigma,dmax,w)
#' apply(tpm,1,sum) # check that rows sum to 1
#' 
#TPM = make.inout.tpm(sigma=exp(theta_new[3])*sqrt(k),dmax=dmax*planespd,w=w) 
make.inout.tpm=function(sigma,dmax,w,nx=500){
  if(sigma==0) {
  inout = diag(c(1,1))  
  } else {
    inout=matrix(c(
      p.i2i(sigma,dmax,w,nx),
      p.o2i(sigma,dmax,w,nx),
      p.i2o(sigma,dmax,w,nx),
      p.o2o(sigma,dmax,w,nx)
    ),ncol=2)
  }
  colnames(inout)=row.names(inout)=c("in","out")
  return(inout)
}


# *************
make.inout.tpm.array=function(sigma,dmax,w,nx=500){
  ntpms = length(sigma)
  inout = array(rep(NA,4*ntpms),dim=c(2,2,ntpms),
                dimnames=list(c("in","out"),c("in","out"),as.character(1:ntpms)))
  for(i in 1:ntpms) {
    if(sigma[i]==0) {
      inout[,,i] = diag(c(1,1))  
    } else {
      inout[,,i]=matrix(c(
        p.i2i(sigma[i],dmax,w,nx),
        p.o2i(sigma[i],dmax,w,nx),
        p.i2o(sigma[i],dmax,w,nx),
        p.o2o(sigma[i],dmax,w,nx)
      ),ncol=2)
    }
  }
  return(inout)
}




