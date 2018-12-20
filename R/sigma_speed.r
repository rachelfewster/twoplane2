#' @title Average speed of Brownian motion particle
#'
#' @description
#'  Calculates the average speed E(displacement/time) of a particle doing 
#'  Brownian motion with movement rate paramter \code{sigmarate}, over a time
#'  period \code{lag}. 
#'  Uses the fact that the expected value of a Chi random variable is sqrt(2)/gamma(0.5).
#'  
#' @param sigmarate The Brownian motion with movement rate paramter
#' @param lag The time over which the average speed is required.
getspeed = function(sigmarate,lag=1) return(sigmarate*sqrt(2)/(gamma(0.5)*sqrt(lag)))

#' @title Brownian motion sigma from average speed and time
#'
#' @description
#'  Calculates the Brownian motion with movement rate paramter from the average 
#'  speed E(displacement/time) of a particle doing Brownian motion over a time
#'  period \code{lag}. 
#'  Uses the fact that the expected value of a Chi random variable is sqrt(2)/gamma(0.5).
#'  NB: Assumes that speed is in km/sec and returns Brownian motion rate paramter in m/sec.
#'  
#' @param sigmarate The Brownian motion with movement rate paramter
#' @param lag The time over which the average speed is required.
speed2sigmarate = function(speed,lag) return(speed*gamma(0.5)*sqrt(lag)/sqrt(2))


#' @title Converts between Palm sigma and LCE sigma
#'
#' @description
#'  Converts between LCE sigma and Palm sigma.
#'  
#' @param sigmarate The Brownian motion with movement rate paramter
#' @param lag The time lag.
sigmarate2sigmapalm = function(sigmarate,lag) return(sigmarate*sqrt(lag)/2)

#' @title Converts between Palm sigma and LCE sigma
#'
#' @description
#'  Converts between LCE sigma and Palm sigma.
#'  
#' @param sigmarate The Palm distibution sigma
#' @param lag The time lag.
sigmapalm2sigmarate = function(sigma,lag) return(sigma*2/sqrt(lag))

