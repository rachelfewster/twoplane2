#' @examples 
#' sigma.beta=0.24
#' mu.beta= -1.6
#' logn.seci(mu.beta,sigma.beta)
#' 
logn.seci=function(mu,sd) {
  se = sqrt((exp(sd^2)-1)*exp(2*mu + sd^2))
  bnd = exp(mu+c(-1,1)*1.96*sd)
  return(list(se=se,lower=bnd[1],upper=bnd[2]))
}