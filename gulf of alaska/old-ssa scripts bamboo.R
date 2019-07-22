# Test scripts doing SSA with R

plot(LakeHuron, main = 'Lake Huron Water Levels', xlab = 'Lake Level (ft)')
plot(diff(LakeHuron), main = 'Lake Huron Water Levels: First Difference',
     xlab = 'Lake Level (ft)')

plot(diff(coral2))

str(LakeHuron)
str(coral2)
cc <- ts(coral2$`Sr/Ca 3808-3D2-Alpha T1`, start=c(2002), end=c(1927), frequency=12)
cx <- rev(coral2$`Sr/Ca 3808-3D2-Alpha T1`)
cx.t <- ts(cx, start=c(1925), end=c(2002), frequency=12)
plot(cx.t)
plot(diff(cx.t))
acf(diff(cx.t))
acf(diff(LakeHuron))

# For the LakeHuron dataset, we create trajectory matrices of 10 (m = 10) and 25 (m = 25) year windows
# note: for the matrix, (n - m + 1) X m, so:
traj10 <- matrix(nrow = length(diff(LakeHuron)) - 10 + 1, ncol = 10)
traj25 <- matrix(nrow = length(diff(LakeHuron)) - 25 + 1, ncol = 25)
# nrow is the "length" (or number of observations) of diff(LakeHuron)
View(traj10)
# As you can see, right now these matrices are empty of values, so we need to fill them
# We can use a "for loop" to fill these with values
for (i in 1:nrow(traj10)) {
  # literally, for every possible value in each row of traj10, from beginning to end ('nrow')
  for (j in 1:ncol(traj10)) {
    # same as above, but with columns
    traj10[i, j] <- diff(LakeHuron)[i + j - 1]
    # fill the empty cells with data from diff(LakeHuron)
  }
}

# now use View(traj10)
View(traj10)

S.traj10 <- (t(traj10) * 1/sqrt(nrow(traj10))) %*% (traj10 * 1/sqrt(nrow(traj10)))
S.traj25 <- (t(traj25) * 1/sqrt(nrow(traj25))) %*% (traj25 * 1/sqrt(nrow(traj25)))
# you can always use View() to see the data you are actually manipulating
# t() transposes your data
# %*% is a matrix multiplication operator, can use ?'%*%' to see what it does

# Spectral decomposition of the lagged covariance matrix (columns are eigenvectors)
S.traj10.eigen <- eigen(S.traj10, symmetric = T)$vectors

# Alternatively,
# Singular Value Decomposition (SVD) of the trajectory matrix
S.traj10.by.svd <- svd(traj10)

S.traj10.eigen <- eigen(S.traj10, symmetric = T)$vectors
dat.traj10 <- as.data.frame(S.traj10.eigen)
colnames(dat.traj10) <- 1:ncol(dat.traj10)
dat.traj10$time <- 1:ncol(dat.traj10)

library(dplyr)
library(tidyr)
library(ggplot2)
dat.traj10 %>%
  gather(key = 'eigenvector', value = 'value', -time) %>%
  mutate(eigenvector = ordered(eigenvector, levels = 1:1000)) %>%
  ggplot(mapping = aes(x = time, y = value)) +
  geom_line(size = 0.8) +
  facet_wrap(~ eigenvector) +
  theme_linedraw()

# from here, if we want to do a smoothing of the time series, we do a reconstruction of the original signal
# using a subset of the components we just produced
library(Rssa)
plot(co2)
s <- ssa(co2, L=120)
plot(s)
# Reconstruct the series, grouping elementary series.
r <- reconstruct(s, groups = list(Trend = c(1, 4), Season1 = c(2,3), Season2 = c(5, 6)))
plot(r)
# 'groups' argument might contain duplicate entries as well
q <- reconstruct(s, groups = list(1, 1:4, 1:6))
plot(q)
w <- reconstruct(s, groups=list(3, 35:37))
plot(w)
plot(s)
plot(s, type="vectors")
plot(s, type="paired")
plot(wcor(s))


plot(cx.t)

cx.ssa <- ssa(cx.t, L=456)
plot(cx.ssa)
plot(cx.ssa, type="vectors")
plot(cx.ssa, type="paired")
plot(wcor(cx.ssa))

#cx.r <- reconstruct(cx.ssa, groups=list(X1 = c(3, 4), X2 = c(5,6), X3 = c(9, 10)))
cx.r <- reconstruct(cx.ssa, groups=list(2, 3:4, 5:6, 9:10))
plot(cx.r)

cx.res <- cx.r$F1
trend <- cx.res
plot(trend)

res.trend <- residuals(cx.r)
spec.pgram(res.trend, detrend = FALSE, log = "no")

# CO2
res2 <- reconstruct(s, groups = list(1))
trend2 <- res2$F1
res.trend2 <- residuals(res2)
spec.pgram(res.trend2, detrend = FALSE, log = "no")

plot(co2)
plot(res.trend2)
plot(trend2)

plot(b4.ts)
plot(ma(b3.ts, order = 12))