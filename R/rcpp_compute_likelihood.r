

rcpp_compute_likelihood <- function(s1, s2, dmax, theta_1, theta_2, theta_3, theta_4, p1, p2, k, tL, planespd, p10k, p01k, p11k){
  x=.Call( "compute_likelihood", s1, s2, dmax, theta_1, theta_2, theta_3, theta_4, p1, p2, k, tL, planespd, p10k, p01k, p11k, PACKAGE="twoplane")
  return(x)
}

#  Testing, access to likelihood function with known recaptures.
calculate_Ljk <- function(s1, s2, dmax, theta_1, theta_2, theta_3, theta_4, p1, p2, k, tL, planespd, p10k, p01k, p11k, varstate){
  x=.Call("calculate_Ljk", s1, s2, dmax, theta_1, theta_2, theta_3, theta_4, p1, p2, k, tL, planespd, p10k, p01k, p11k, varstate, PACKAGE="twoplane")
  return(x)
}
