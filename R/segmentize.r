segmentize=function(s1,s2,dmax){
  ss=c(s1,s2)
  ord=order(ss)
  ds=ss[ord]
  diffs=diff(ds)
  cutafter=which(diffs>dmax)
  cuts=(ds[cutafter]+ds[cutafter+1])/2
  return(cuts)
}

segmentsummary = function(dat,planespd,cutstretch=1) {
  s1 = dat$y1/planespd
  s2 = dat$y2/planespd
  tL = dat$L/planespd
  tw = dat$w/planespd
  tb = dat$b/planespd
  dmax.t = tb-tw # max dist animal can move (in plane seconds)
  cuts=sort(unique(c(0,tL,segmentize(s1,s2,dmax.t*cutstretch))))
  nseg=length(cuts)-1
  ns1 = ns2 = rep(0,nseg)
  for(i in 1:nseg) {
    ss1=s1[cuts[i]<s1 & s1<=cuts[i+1]]
    ss2=s2[cuts[i]<s2 & s2<=cuts[i+1]]
    ns1[i] = length(ss1)
    ns2[i] = length(ss2)
  }
  sumdat = data.frame(segment=1:nseg,n1=ns1,n2=ns2)
  return(sumdat)
}

segmentplot = function(dat,planespd,cutstretch=1) {
  sumdat = segmentsummary(dat,planespd,cutstretch)
  big = sumdat$n1
  small = sumdat$n2
  nseg = length(big)
  for(i in 1:nseg) {
    if(big[i]<small[i]) {
      big[i] = sumdat$n2[i]
      small[i] = sumdat$n1[i]
    }
  }
  pair = rep("",nseg)
  for(i in 1:nseg) pair[i] = paste(big[i],"-",small[i],sep="")
  
  barplot(table(pair))
  
}

# Copy of nch; forgot I wrote nch and rewrote it here; need to tidy up!
numinsegment = function(n1,n2=n1) {
  if(n1<n2) {
    temp = n1
    n1 = n2
    n2 = n1
  }
  nn = 1
  for(m in 1:n2) {
    nn = nn + choose(n2,m)*exp(lgamma(n1+1)-lgamma(n1-m+1))
  }
  return(nn)
}

#' @title Calculates number of possible capture histories
#'
#' @description
#'  Calculates number of possible capture histories, given the number of 
#'  detections by each observer.
#'  
#' @param n1 The number of observer 1 detections.
#' @param n2 The number of observer 2 detections.
nch = function(n1,n2=n1) {
  n = 1
  nbig = n1; nsmall = n2
  if(n2>n1) {nbig=n2; nsmall=n1}
  for(m in 1:nsmall) {
    n = n + choose(nsmall,m) * exp(lgamma(nbig+1)-lgamma(nbig-m+1))
  }
  return(n)
}


