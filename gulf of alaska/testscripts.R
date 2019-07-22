x.spec <- spectrum(data,log="no",span=10,plot=FALSE)
spx <- x.spec$freq/1
spy <- 2*x.spec$spec
plot(spy~spx,xlab="frequency",ylab="spectral density",type="l")

raw.spec <- spec.pgram(data, taper = 0, plot = FALSE, demean = TRUE, na.action = na.omit)
plot(raw.spec, log = "no")
