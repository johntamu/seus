# some cleaned up versions of the scripts from a different file
# import 3808-4
# Look at Disk 2 of #3808-4
View(D2_3808_4)
b4.d2 <- D2_3808_4$`Sr/Ca ALV-3808 #4 disk 1 C2 T1`[1:697]
# remember that analyses start from the outside, so need to reverse
b4.d2.rev <- rev(b4.d2)
# make it a time series
b4.ts <- ts(b4.d2.rev, start=c(1938), end=c(2002), frequency=12)
plot(b4.ts)
plot(diff(b4.ts))
b4.ssa <- ssa(b4.ts, L = 36)
plot(b4.ssa)
plot(b4.ssa, type="vectors")
plot(b4.ssa, type="paired")
plot(wcor(b4.ssa))
recon.b4 <- reconstruct(b4.ssa, groups=list(2, 3))
plot(recon.b4)
trend <- recon.b4$F3
plot(trend)
res.trend <- residuals(recon.b4)
spec.pgram(res.trend, detrend = FALSE, log = "no")
plot(res.trend)

t.fft <- fft(npi)
plot.frequency.spectrum(t.fft)

plot.frequency.spectrum <- function(X.k, xlimits=c(0,length(X.k))) {
  plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))
  
  # TODO: why this scaling is necessary?
  plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2] 
  
  plot(plot.data, t="h", lwd=2, main="", 
       xlab="Frequency (Hz)", ylab="Strength", 
       xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
}

plot.frequency.spectrum(t.fft)

get.trajectory <- function(X.k,ts,acq.freq) {
  
  N   <- length(ts)
  i   <- complex(real = 0, imaginary = 1)
  x.n <- rep(0,N)           # create vector to keep the trajectory
  ks  <- 0:(length(X.k)-1)
  
  for(n in 0:(N-1)) {       # compute each time point x_n based on freqs X.k
    x.n[n+1] <- sum(X.k * exp(i*2*pi*ks*n/N)) / N
  }
  
  x.n * acq.freq 
}

get.trajectory(t.fft, npi.ts)

plot.harmonic <- function(Xk, i, ts, acq.freq, color="red") {
  Xk.h <- rep(0,length(Xk))
  Xk.h[i+1] <- Xk[i+1] # i-th harmonic
  harmonic.trajectory <- get.trajectory(Xk.h, ts, acq.freq=acq.freq)
  points(ts, harmonic.trajectory, type="l", col=color)
}
 plot.harmonic(t.fft, 1)
 
 
library(spectral)
spectrum(b3.d2)
spectrum(b4.d2)

b3.1000 <- b3.d2*1000
b3.anom <- b3.d2-mean(b3.d2)
View(b3.1000)
spectrum(b3.anom)

periodogram(b3.anom)

spec.ar(b3.anom)
spec.ar(b4.ts)

plot(reconstruct(b4.cd.ssa, "series", groups=c(1,2)))
b4.recon1 <- reconstruct(b4.cd.ssa, "series", groups=c(1,2))
plot(b4.recon1$F1)
abline(v=c(1927, 1977), col="red", lty=5)
abline(v=1947, col="blue")
abline(h=0.003139172, lty=5)