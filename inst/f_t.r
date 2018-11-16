f.t = function(t,sigma=10,k=100,v=100) {
  v*k*exp(-v^2*(k-t)^2/(2*sigma^2*t))/sqrt(2*pi*sigma^2*t^3)
}

sigma = 10
v = 100
k=200
dt = 8
t=seq(k-dt,k+dt,length=200)
t=seq(k-dt,k+dt,length=200)
correct = f.t(t,sigma=sigma,v=v,k=k)
simple = dnorm(t,mean=k,sd=sigma/sqrt(v))
quartz()
plot(t,correct,type="l",ylim=range(correct,simple),ylab=expression(f[T](t)),lwd=2)
lines(t,simple,type="l",lty=2,lwd=1.5,col="gray")

