% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pomega.r
\name{dft}
\alias{dft}
\title{PDF of Brownian hitting time}
\usage{
dft(ft, k, spd, theta3)
}
\arguments{
\item{k}{The hitting time if the particle did not move.}

\item{spd}{The speed at which the barrier drifts.}

\item{theta3}{The log of the Brownian motion rate parameter.}

\item{tt}{The hitting time}
}
\description{
Calculates the PDF of the time to hit a barrier for a particle moving with 1D 
 Brownian motion with movement rate paramter \code{exp(theta3)} and drift 
 speed \code{spd} towards the barrier, when the barrier is initially at a
 distance \code{k*spd} away.
}
