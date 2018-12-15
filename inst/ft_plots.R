# test that dft() works:
# ---------------------
planeknots=100 # observer speed in knots
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec
k=248 # time separation of observers in seconds
#sigmaknots=7.5 
#sigmarate=sigmaknots*nm2km/60/60 # std dev of movement distance per second (in km/sec)
animalspd =0.3*planespd
sigmarate = speed2sigmarate(animalspd,k)
sigma.mult = 5
sigma.t = sqrt(k)*sigmarate/planespd
fts = seq(max(0,k-sigma.mult*sigma.t),k+sigma.mult*sigma.t,length=500)
plot(fts,dft(fts,k,planespd,log(sigmarate)),type="l",ylab=expression(f[T](t)),xlab="t")
lines(fts,dnorm(fts,mean=k,sd=sigma.t),lty=2,col="blue")

# Check f_T(t) simulator rft:
randts = rft(100000,k,planespd,log(sigmarate))
hist(randts,probability=TRUE)
lines(fts,dft(fts,k,planespd,log(sigmarate)))
lines(fts,dnorm(fts,mean=k,sd=sigma.t),lty=2,col="blue")
# Conclusion: rft has slightly too high variance
sd(randts)
#[1] 4.834763
sigma.t
#[1] 5.831533