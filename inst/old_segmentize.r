segmentize=function(s1,s2,dmax){
  #--------------------------------------------------------------------------------------------------
  # Returns endpoints of all segments containing only one set of linked possible pairings
  #
  # Depends: mergeranges(), pairdists()
  #
  #--------------------------------------------------------------------------------------------------
  # cat("sigma,k,dmax: ",sigma,", ",k,", ",dmax,"\n")
  ds=pairdists(s1,s2) # distances between s1 and s2 detections
  db=ds<=dmax # logical indicating which distances possible pairs
  joined1=which(apply(db,1,sum)>0) # indices of s1 detections that have at least 1 possible pairing
  njoined1=length(joined1) # number of s1 detections that are joined
  ranges=matrix(rep(NA,njoined1*2),ncol=2)
  for(i in 1:njoined1) {
    s2s=which(db[joined1[i],])
    ranges[i,]=range(s1[joined1[i]],s2[s2s]) # ranges of locations of possible pairs
  }
  pairanges=mergeranges(ranges)
  nranges=length(pairanges[,1])
  cuts=sort(c(0,s1,s2)) + dmax/10
  outs=c()
  for(i in 1:nranges) {
    out=which(pairanges[i,1]<=cuts & cuts<=pairanges[i,2])
    if(!is.null(out)) outs=c(outs,out)
  }
  cuts=cuts[-outs]
  #print(paste("cuts: ",cuts))
  return(cuts)
}  


mergeranges=function(lrange){
  #--------------------------------------------------------------------------------------------------
  # For each detection that has a detection by the other platform close enough to be a possible pair, 
  # returns min dist (col 1) and max dist (col 2) over all such pairs.
  #
  # Depends: overlap()
  #
  #--------------------------------------------------------------------------------------------------
  n=dim(lrange)[1] # number of ranges in starting set (i.e. number of rows in lrange)
  inds=1:n
  left=rep(TRUE,n) # rows that have not been merged with another row
  for(i in 1:n){
    if(left[i]) {
      leaveout=NULL # index for rows out of consideration for this i
      for(j in inds[left][inds[left]>i]){
        if(overlap(lrange[i,],lrange[j,])) { # if ranges overlap
          lrange[i,]=range(lrange[i,],lrange[j,]) # extend range i to include range j
          if(is.null(leaveout)) leaveout=j
          else leaveout=unique(c(leaveout,j)) # mark range j for removal from set of unmerged ranges
          #        cat("i= ",i);cat("j= ",j);cat(paste("leaveout: ",leaveout)); cat("\n")
        }
      }
    }
    if(!is.null(leaveout)) left[leaveout]=FALSE # remove those ranges that were merged with range i from set of ranges left
    #    cat("i= ",i);cat(paste("left: ",left)); cat("\n")
  }
  return(lrange[left,,drop=FALSE])
}


overlap=function(r1,r2){
  #--------------------------------------------------------------------------------------------------
  # Returns true if the ranges r1 and r1 overlap; else returns FALSE.
  #--------------------------------------------------------------------------------------------------
  in1=r1[1]<=r2[1] & r2[1]<=r1[2]
  in2=r2[1]<=r1[1] & r1[1]<=r2[2]
  return(in1 || in2)
}

pairdists=function(s1,s2){
  #--------------------------------------------------------------------------------------------------
  # Calculate distances between all detection locations in s1 and those in s2 
  #--------------------------------------------------------------------------------------------------
  n1=length(s1)
  n2=length(s2)
  n=n1+n2
  ss=matrix(c(c(s1,s2),rep(0,n)),ncol=2)
  dmat=as.matrix(dist(ss))[1:n1,(n1+1):(n1+n2)]  
}