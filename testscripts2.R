library(lattice)
library(latticeExtra)
library(dplyr)

View(df.bulk)

df.bulk %>%
  filter(coral.id == "jack-4684-bc-unk" | coral.id == "jack-4907-bc1-d1" |
           coral.id == "jack-4907-bc1-d3" | coral.id == "sav-4902-bc1-unk" |
           coral.id == "stet-4904-bc1-d2") -> bulk2

xyplot() # plan is to create an xyplot of each of these corals to see what the "overall" record looks like

View(corals)
plot(TP ~ bp, corals)
max(corals$TP)
min(corals$TP)
mean(corals$TP)
sd(corals$TP)

corals %>%
  filter(Sample.ID2 == "Jacksonville-4684") -> t.coral
corals %>%
  filter(Sample.ID2 == "Jacksonville-4907" | Sample.ID2 == "Savannah Banks-4902") -> t.coral2
corals %>%
  filter(Sample.ID2 == "Savannah Banks-4902") -> t.coral3

mean(t.coral$TP)
sd(t.coral$TP)
mean(t.coral2$TP)
sd(t.coral2$TP)

max(corals$Sum.V)
min(corals$Sum.V)
mean(corals$Sum.V)
sd(corals$Sum.V)

plot(Sum.V ~ bp, corals)
corals <- droplevels(corals)
xyplot(Sum.V ~ bp, data = corals,
       group = Sample.ID2, auto.key = TRUE)

# Plotting D14C vs distance for each of the corals

par(mfrow = c(2,3))
plot(D14C ~ Distance.microns, r.jack, type = "o", xlab = "Distance (microns)", main = "jack4907")
plot(D14C ~ Distance.microns, r.jack2, type = "o", xlab = "Distance (microns)", main = "jack4684")
plot(D14C ~ Distance.microns, r.sav, type = "o", xlab = "Distance (microns)", main = 'sav4902')
plot(D14C ~ Distance.microns, r.stet, type = "o", xlab = "Distance (microns)", main = 'stetson4904 second')
plot(D14C ~ Distance.microns, r.stet2, type = "o", xlab = "Distance (microns)", main = 'stetson4904 first')

par(mfrow = c(1,2))
plot(D14C ~ Distance.microns, r.jack2, type = "o", xlab = "Distance from edge (um)", ylab = expression({Delta}^14*"C"))
plot(D14C ~ Distance.microns, r.stet, type = "o", xlab = "Distance from edge (um)", ylab = expression({Delta}^14*"C"))

par(mfrow = c(1,2))
plot(D14C ~ Distance.microns, r.jack, type = "o", xlab = "Distance (microns)", ylab = expression({Delta}^14*"C"))
plot(D14C ~ Distance.microns, r.sav, type = "o", xlab = "Distance (microns)", ylab = expression({Delta}^14*"C"))


View(t.stet)
ttt <- detrend(t.stet$d15n)
ttc <- detrend(t.stet$d13c)
t.stet$n.anom <- ttt$Y
t.stet$c.anom <- ttc$Y

plot(n.anom ~ bp, t.stet, type="l")
plot(c.anom ~ bp, t.stet, type="l")

par(mfrow = c(2,1))
plot(d15n ~ bp, t.stet, type="l")
plot(d13c ~ bp, t.stet, type="l")



par(mfrow=c(2,1))
plot(d13c ~ linear.ad, t.stet,
                       type = "l",
                       cex = 0.5,
                       col = alpha("black", 0.3),
                       # axes = FALSE,
                       ylab = "",
                       yaxt="n",
                       xlab = "Years BP",
                       bty = "l",
                       xlim = c(2005, 1700)
                       # ylim = c(-16.5,-15)
)
  minor.tick(nx = nx)
  axis(side = 2)
  mtext(c, side = 2, line=2.5)
  lines(rollmean(d13c, 5, fill = NA, align = "center") ~ linear.ad, t.stet,
        col = "black")

plot(d15n ~ linear.ad, t.stet,
       type = "l",
       cex = 0.5,
       col = alpha("black", 0.3),
       # axes = FALSE,
       ylab = "",
       yaxt="n",
       xlab = "Years BP",
       bty = "l",
       xlim = c(2005, 1700)
       # ylim = c(-16.5,-15)
  )
minor.tick(nx = nx)
axis(side = 2)
mtext(n, side = 2, line=2.5)
lines(rollmean(d15n, 5, fill = NA, align = "center") ~ linear.ad, t.stet,
        col = "black")

obj2 <- xyplot(rollmean(d15n, 5, fill = NA, align = "center") ~ linear.ad, 
               data = t.stet, 
               # key=list(space="right",
               #          lines=list(col=c("#4d4d4d","#878787"), lty=c(1,1), cex = 0.1),
               #          text=list(c("A"," B"))),
               type = "l", lwd = 1.25, ylab = n,  xlim = c(2100, 400), col = "#D55E00", xlab = x)
obj1 <- xyplot(rollmean(d13c, 5, fill = NA, align = "center") ~ linear.ad, data = t.stet, type= "l", lwd = 1.25, ylab = c, col = "#0072B2", xlim = c(2100, 400))
doubleYScale(obj2, obj1, add.ylab2 = TRUE, style1= NULL, style2 = NULL)

obj2 <- xyplot(rollmean(d15n, 5, fill = NA, align = "center") ~ bp, 
               data = t.jack, 
               # key=list(space="right",
               #          lines=list(col=c("#4d4d4d","#878787"), lty=c(1,1), cex = 0.1),
               #          text=list(c("A"," B"))),
               type = "l", lwd = 1.25, ylab = n, col = "#D55E00", xlab = "year BP")
obj1 <- xyplot(rollmean(d13c, 5, fill = NA, align = "center") ~ bp, data = t.jack, type= "l", lwd = 1.25, ylab = c, col = "#0072B2")
doubleYScale(obj2, obj1, add.ylab2 = TRUE, style1= NULL, style2 = NULL)


# Linear regression on subset of Stetson peels
t.stet %>%
  filter(linear.ad > 1699 & linear.ad < 1801) -> t.df

summary(lm(d13c ~ d15n, t.df))
plot(d13c ~ d15n, t.df)


t.bulk <- df.bulk
t.bulk %>%
  filter(coral.id == "jack-4684-bc-unk" | coral.id == "jack-4907-bc1-d3" | coral.id == "sav-4902-bc1-unk" |
           coral.id == "stet-4904-bc1-d2") -> t.bulk

t.bulk %>%
  filter(coral.id == "jack-4907-bc1-d3" |
           coral.id == "stet-4904-bc1-d2") -> t.bulk2

droplevels(t.bulk$coral.id)
xyplot(rollmean(d15n, 3, fill = NA, align = "center") ~ linear.ad,
       group = coral.id,
       data = t.bulk,
       type = "l",
       ylab = n,
       xlab = x,
       # auto.key = TRUE,
       par.settings = list(superpose.line = list(lwd = 1.5, col = c('#33a02c', '#b2df8a', '#1f78b4', '#a6cee3'))))


f <- ggplot(t.bulk, aes(x=linear.ad, y=rollmean(d15n, 3, fill = NA, align = "center")))
f + geom_line(aes(color = coral.id)) +
                xlab(x) +
                ylab(n) +
  scale_color_manual(values=c("#999999", '#009E73', '#0072B2', '#D55E00'), 
                     labels = c("Jacksonville-4684 BC1", "Jacksonville-4907 BC1", "Savannah-4902 BC1", "Stetson-4904 BC1"), 
                     name = NULL) +
  scale_x_continuous(breaks = c(-1000, -500, 0, 500, 1000, 1500, 2000)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

d <- ggplot(t.bulk2, aes(x=linear.ad, y=rollmean(d13c, 3, fill = NA, align = "center")))
d + geom_line(aes(color = coral.id), size = 1) +
  xlab(x) +
  ylab(c) +
  scale_color_manual(values=c('#1f78b4', '#1b9e77'), labels = c("A", "B", "C", "D"), 
                     name = NULL) +
  scale_x_continuous(breaks = c(-1000, -500, 0, 500, 1000, 1500, 2000)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))
  
obj2 <- xyplot(rollmean(d15n, 5, fill = NA, align = "center") ~ bp, data = t.jack, type = "l", lwd = 2, ylab = n,  xlim = c(2005, 1500), col = "#56B4E9", xlab = x)
obj1 <- xyplot(rollmean(d13c, 5, fill = NA, align = "center") ~ linear.ad, data = t.stet, type= "l", lwd = 2, ylab = c, col = "#009E73", xlim = c(2005, 1500))
doubleYScale(obj2, obj1, add.ylab2 = TRUE, style1= NULL, style2 = NULL)




xlim1 <- c(2000, 3100)
par.set <- par(mfrow=c(2,1), oma = c(2, 2, 0, 0), # two rows of text at the outer left and bottom margin
               mar = c(1, 1, 0, 0)+0.3) # space for one row of text at ticks and to separate plots
# mgp = c(2, 1, 0))   # axis label at 2 rows distance, tick labels at 1 row

plot(d13c ~ bp, t.jack,
     type = "l",
     cex = 0.25,
     col = alpha("black", 0.3),
     # axes = FALSE,
     ylab = c,
     yaxt="l",
     xaxt = "n",
     las = 1,
     xlab = "Years BP",
     bty = "n",
     xlim = xlim1,
     ylim = c(-17,-14))
# minor.tick(nx = nx)
# axis(side = 2, las = 1)
# mtext(c, side = 2, line=2.5)
# mtext("Years BP", side = 1, line = 0.8, outer = TRUE, cex = 0.85, font = 2)
lines(rollmean(d13c, 5, na.pad = TRUE) ~ bp, t.jack,
      col = "black")

plot(d15n ~ bp, t.jack,
     type = "l",
     cex = 0.25,
     col = alpha("black", 0.3),
     # axes = FALSE,
     ylab = n,
     yaxt="l",
     las = 1,
     xlab = "Years BP",
     bty = "l",
     xlim = xlim1)
minor.tick(nx = nx)
# axis(side = 2, las = 1)
# mtext(c, side = 2, line=2.5)
mtext("Years BP", side = 1, line = 0.8, outer = TRUE, cex = 0.85, font = 2)
lines(rollmean(d15n, 5, na.pad = TRUE) ~ bp, t.jack,
      col = "black")
abline(v=2850, lty = "dotted")
abline(v=2250, lty = "dotted")

#4575b4

plot(d13c ~ linear.ad, t.jack,
     xlab = x,
     ylab = c,
     type = "l",
     las = 1,
     cex = 0.5,
     col = alpha("#4575b4", 0.3),
     xlim = c(c(-1100, 2005)),
     ylim = c(-17, -14.5))
lines(rollmean(d13c, 5, fill = NA, align = "center") ~ linear.ad, df.jack,
      col = "#4575b4", lwd = 1.75)
lines(d13c ~ linear.ad, t.stet, col = alpha("#CC79A7", 0.3))
lines(rollmean(d13c, 5, fill = NA, align = "center") ~ linear.ad, t.stet, col = "#CC79A7", lwd = 1.75)
lines(d13c ~ linear.ad, t.sav, col = alpha("#D55E00", 0.3))
lines(rollmean(d13c, 5, fill = NA, align = "center") ~ linear.ad, t.sav, col = "#D55E00", lwd = 1.75)

plot(d15n ~ linear.ad, t.stet,
     xlab = x,
     ylab = c,
     type = "l",
     las = 1,
     cex = 0.5,
     col = alpha("#4575b4", 0.3),
     xlim = c(1200, 2005))
lines(rollmean(d15n, 5, fill = NA, align = "center") ~ linear.ad, t.stet, col = "#4575b4", lwd = 1.75)
moy.melt %>%
  filter(variable == 'Red.Color.Intensity.Units') %>%
  na.omit() %>%
  plot(value ~ yrAD, .,
       type = "l",
       bty = "l",
       col = alpha("black", 0.3),
       xlab = "Years BP",
       ylab = "Red Intensity",
       xlim = c(1300, 2005))
       
moy.melt %>%
  filter(variable == 'Red.Color.Intensity.Units') %>%
  na.omit() %>%
  lines(rollmean(value, 100, na.pad = TRUE) ~ yrAD, ., col = "black", xlim = c(1300, 2005))

trouet.melt %>%
  filter(variable == 'NAOms') %>%
  na.omit() %>%
  plot(value ~ yrAD, .,
       type = "l",
       bty = "l",
       col = alpha("black", 0.3),
       xlab = "Years",
       # main = "Trouet NAO",
       xlim = c(1200, 2005))
trouet.melt %>%
  filter(variable == 'NAOms') %>%
  na.omit() %>%
  lines(rollmean(value, 100, na.pad = TRUE) ~ yrAD, ., col = "black", xlim = c(1200, 2005))

trouet.melt %>%
  filter(variable == 'NAOms') %>%
  na.omit() %>%
  plot(value ~ yrAD, .,
       type = "l",
       bty = "l",
       col = alpha("black", 0.3),
       xlab = "Years",
       # main = "Trouet NAO",
       xlim = c(1200, 2005))
trouet.melt %>%
  filter(variable == 'NAOms') %>%
  na.omit() %>%
  lines(rollmean(value, 100, na.pad = TRUE) ~ yrAD, ., col = "black", xlim = c(1200, 2005))

# Test data frame to create an age model with a single linear growth rate

