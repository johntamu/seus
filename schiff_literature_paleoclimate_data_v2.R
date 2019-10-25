#' -----------------------------------------------------
#' Topic: Compiled paleoclimate data from other sources
#' John Schiff
#' 08-14-2019
#' -----------------------------------------------------
#' 
#' This script file is meant for importing, analyzing,
#' and combining paleoclimate data with the isotope data
#' from the black corals and bamboo corals. Different
#' datasets are imported and organized by the study
#' utilizing it.
#' 
#' Other datasets can be utilized. I got these from the
#' NOAA paleoclimate database.
#' 

library(zoo)
library(ggplot2)
library(dplyr)
library(grid)
library(reshape2)
library(pryr)

# Calculate yrs BP for each of the corals
t.stet <- df.stet
t.sav <- df.sav
t.jack <- df.jack
t.jack4684 <- df.jack4684

t.stet$bp <- 1950 - t.stet$linear.ad
t.sav$bp <- 1950 - t.sav$linear.ad
t.jack$bp <- 1950 - t.jack$linear.ad
t.jack4684$bp <- 1950 - t.jack4684$linear.ad
ndata$bp <- 1950 - ndata$Year.CE

# Keiwgin 1996; Bermuda Rise
string1 <- '~/Documents/GitHub/data/paleoclimate_data/bc004a-tab.txt'
string2 <- '~/Documents/GitHub/data/paleoclimate_data/bc004d-tab.txt'
keigwin4a <- read.delim(string1, header = TRUE, comment.char = '#')
keigwin4d <- read.delim(string2, header = TRUE, comment.char = '#')

# Richey 2009; Gulf of Mexico
fisk <- read.csv('~/Documents/GitHub/data/paleoclimate_data/richey2009-fisk.csv', header = TRUE)
garrison <- read.csv('~/Documents/GitHub/data/paleoclimate_data/richey2009-garrison.csv', header = TRUE)

# Richey 2007; Gulf of Mexico
richey2007 <- read.csv('~/Documents/GitHub/data/paleoclimate_data/richey2007.csv', header = TRUE)

# Schmidt 2012; Gulf of Mexico
schmidt <- read.delim('~/Documents/GitHub/data/paleoclimate_data/schmidt2012.txt', header = TRUE, comment.char = '#')

# Saenger 2011; Carolina Slope
saenger.core1 <- '~/Documents/GitHub/data/paleoclimate_data/saenger2011_core1.csv'
saenger.core2 <- '~/Documents/GitHub/data/paleoclimate_data/paleoclimate_data/saenger2011_core2.csv'
core1 <- read.csv(saenger.core1)
core2 <- read.csv(saenger.core2)

# Saenger 2009; Bahamas
bahamas.sst <- read.csv('~/Documents/GitHub/data/paleoclimate_data/bahamas2009sst.csv')

# Lund and Curry 2006; Great Bahama Bank


# Trouet 2009; NAO z-score
trouet <- read.delim('~/Documents/GitHub/data/paleoclimate_data/nao-trouet2009.txt', header = TRUE, comment.char = '#', sep = "")
trouet$z <- (trouet$NAOms - mean(trouet$NAOms)) / sd(trouet$NAOms)
# ***********************
# Build a new data frame
#
#
# ***********************

schmidt$age <- schmidt$age_calkaBP*1000

schmidt.melt <- melt(schmidt, "age")
richey2007.melt <- melt(richey2007, "Cal.yr.B.P.")
fisk.melt <- melt(fisk, "yrBP")
core1.melt <- melt(core1, "Year.AD.")
keigwin1.melt <- melt(keigwin4a, "yrBP")
keigwin2.melt <- melt(keigwin4d, "yrBP")
bahamas.melt <- melt(bahamas.sst, "Year..A.D.")
trouet.melt <- melt(trouet, "Year")

colnames(schmidt.melt)[colnames(schmidt.melt)=="age"] <- "yrBP"
colnames(richey2007.melt)[colnames(richey2007.melt)=="Cal.yr.B.P."] <- "yrBP"
colnames(fisk.melt)[colnames(fisk.melt)=="yrBP"] <- "yrBP"
colnames(core1.melt)[colnames(core1.melt)=="Year.AD."] <- "yrAD"
colnames(keigwin1.melt)[colnames(keigwin1.melt)=="yrBP"] <- "yrBP"
colnames(keigwin2.melt)[colnames(keigwin2.melt)=="yrBP"] <- "yrBP"
colnames(bahamas.melt)[colnames(bahamas.melt)=="Year..A.D."] <- "yrAD"
colnames(trouet.melt)[colnames(trouet.melt)=="Year"] <- "yrAD"

core1.melt$yrBP <- 1950 - core1.melt$yrAD
bahamas.melt$yrBP <- 1950 - bahamas.melt$yrAD
trouet.melt$yrBP <- 1950 - trouet.melt$yrAD

#' ------------------------
#' Now let's plot them
#' all together
#' ------------------------

bp <- "Years before 1950 CE"
sst <- "SST (deg C)"
lim <- c(600,2000)
{
plotrichey <- richey2007.melt %>%
  filter(variable == 'SST') %>%
  na.omit() %>%
  ggplot(aes(x = yrBP, y = rollmean(value, 5, na.pad = TRUE))) +
  geom_line(color = "#bd0026") +
  geom_line(aes(x = yrBP, y = value), color = alpha("black", 0.3)) +
    annotate(geom = "text", x = 200, y = 26.25, label = "GOM sediment core (Richey et al, 2007)", size = 3) +
    theme_bw() +
    theme(axis.text.y   = element_text(size=10, color = "black"),
          axis.text.x   = element_blank(),
          axis.title.y  = element_text(size=10),
          axis.title.x  = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_line(size = 0.25, color = "black"),
          axis.ticks.y = element_line(size = 0.25, color = "black"),
          panel.background = element_rect(size = 0.0, color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab(bp) +
  xlab(bp) +
  ylab(sst) +
  xlim(lim)
plotfisk <- fisk.melt %>%
  filter(variable == "SST") %>%
  na.omit() %>%
  ggplot(aes(x = yrBP, y = value), size = 0.5, alpha = 0.75) +
  geom_line(color = "black", size = 0.75) +
  annotate(geom = "text", x = 200, y = 26.25, label = "GOM sediment core (Richey et al, 2009)", size = 3) +
  theme_bw() +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_text(size=10, color = "black"),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10, color = "black"),
        axis.line.x = element_blank(),
        axis.ticks.x = element_line(size = 0.25, color = "black"),
        axis.ticks.y = element_line(size = 0.25, color = "black"),
        panel.background = element_rect(size = 0.0, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(bp) +
  ylab(sst) +
  xlim(lim)
stet.n <- ggplot(data=t.stet, aes(x = bp, y = rollmean(d15n, 15, na.pad = TRUE))) +
  geom_line(color = "black") +
  geom_line(aes(bp, y = d15n), color = alpha("black", 0.3)) +
  # geom_line(data=t.jack4684, aes(x=bp, y=rollmean(d15n, 5, na.pad = TRUE)), color = "red") +
  # geom_line(data=t.jack4684, aes(bp, y = d15n), color = alpha("red", 0.3)) +
  # geom_point() +
  theme_bw() +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_blank(),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_line(size = 0.25, color = "black"),
        axis.ticks.y = element_line(size = 0.25, color = "black"),
        panel.background = element_rect(size = 0.0, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(bp) +
  ylab(n) +
  xlim(lim) +
  ylim(7, 10)
stet.c <- ggplot(data=t.stet, aes(x = bp, y = rollmean(d13c, 15, na.pad = TRUE))) +
  geom_line(color = "black") +
  geom_line(aes(bp, y = d13c), color = alpha("black", 0.3)) +
  # geom_point() +
  theme_bw() +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_blank(),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_line(size = 0.25, color = "black"),
        axis.ticks.y = element_line(size = 0.25, color = "black"),
        panel.background = element_rect(size = 0.0, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(bp) +
  ylab(c) +
  xlim(lim) +
  ylim(-16, -15)
jack4684.n <- ggplot(data=t.jack4684, aes(x = bp, y = rollmean(d15n, 5, na.pad = TRUE))) +
  geom_line(color = "black") +
  geom_line(aes(bp, y = d15n), color = alpha("black", 0.3)) +
  theme_classic() +
  xlab(bp) +
  ylab(n) +
  xlim(lim)
plot.jackn <- ggplot(data=t.jack, aes(x = bp, y = rollmean(d15n, 3, na.pad = TRUE))) +
  geom_line(color = "blue") +
  theme_classic() +
  xlab(bp) +
  ylab(n) +
  xlim(lim)
plot.jackc <- ggplot(data=t.jack, aes(x = bp, y = rollmean(d13c, 3, na.pad = TRUE))) +
  geom_line(color = "blue") +
  theme_classic() +
  xlab(bp) +
  ylab(c) +
  xlim(lim)
plot.savn <- ggplot(data=t.sav, aes(x = bp, y = d15n)) +
  geom_line(color = "blue") +
  theme_classic() +
  xlab(bp) +
  ylab(n) +
  xlim(lim)
plot9 <- ndata %>%
  filter(Region == 'SEUS') %>%
  ggplot(data = ., aes(bp, y = Sum.V)) +
  geom_point(color = "blue") +
  theme_classic() +
  xlab(bp) +
  ylab(n) +
  xlim(lim)
plotk1 <- keigwin1.melt %>%
  filter(variable == 'carb.') %>%
  ggplot(aes(x = yrBP, y = value)) +
  geom_point(color = "red") +
  theme_classic() +
  xlab(bp) +
  xlim(lim)
plotschmidt <- schmidt.melt %>%
  filter(variable == 'sst') %>%
  ggplot(aes(x = yrBP*1000, y = value)) +
  geom_point(color = "red") +
  geom_line(color = "red") +
  theme_classic() +
  ylab(sst) +
  xlab(bp) +
  xlim(lim)
plotsaenger <- core1.melt %>%
  # filter(variable == 'SST.Anand.') %>%
  filter(variable == 'SSS') %>%
  # filter(variable == 'd18Oc') %>%
  # filter(variable == 'd18Osw') %>%
  # filter(variable == 'd18Ow') %>%
  # filter(variable == 'SST') %>%
  ggplot(aes(x = yrBP, y = value)) +
  # geom_point(color = "red") +
  geom_line(color = "red") +
  theme_classic() +
  ylab(sst) +
  xlab(bp) +
  xlim(lim)
plotbahamas <- bahamas.melt %>%
  filter(variable == 'SST.anomaly') %>%
  ggplot(aes(x = yrBP, y = rollmean(value, 20, na.pad = TRUE))) +
  geom_line(color = "#bd0026") +
  geom_line(aes(x = yrBP, y = value), color = alpha("black", 0.3)) +
  theme_classic() +
  ylab("SST anomaly (deg C)") +
  xlab(bp) +
  xlim(lim)
plottrouet <- trouet.melt %>%
  filter(variable == 'NAOms') %>%
  ggplot(aes(x = yrBP, y = rollmean(value, 50, na.pad = TRUE))) +
  geom_line(color = "#bd0026") +
  geom_line(aes(x = yrBP, y = value), color = alpha("black", 0.3)) +
  theme_classic() +
  ylab("NAO z-score") +
  xlab(bp) +
  xlim(lim)
}

grid.newpage()
grid.draw(rbind(ggplotGrob(stet.n), ggplotGrob(plot.jackn), ggplotGrob(plotfisk), size = "last")) #  ggplotGrob(stet.c)

#' -------------------------------------------------
#' Comparing bulk record data to
#' paleoclimate data. This is cumbersome
#' but essentially: Determine whether paleo-record 
#' is on yrs BP or year AD/CE timescale. Then go
#' from there.
#' -------------------------------------------------
#' 

# Redoing the above but with base R code
ma1 <- 3
ma2 <- 5
ma3 <- 10

par.set <- par(pty = "s", mfrow=c(2,1), oma = c(2, 2, 0, 0), # two rows of text at the outer left and bottom margin
               mar = c(1, 1, 0, 0)+0.3) # space for one row of text at ticks and to separate plots
# mgp = c(2, 1, 0))   # axis label at 2 rows distance, tick labels at 1 row

pryr.jackn %<a-% {plot(d15n ~ bp, t.jack,
                    type = "l",
                    cex = 0.5,
                    col = alpha("black", 0.3),
                    axes = FALSE,
                    ylab = "",
                    yaxt="l",
                    # xlab = "Years BP",
                    xlim = c(1750, 3000))
  axis(side = 2)
  mtext(n, side = 2, line=2.5)
  lines(rollmean(d15n, 5, na.pad = TRUE) ~ bp, t.jack,
        col = "black")
}
pryr.jackc %<a-% {plot(d13c ~ bp, t.jack,
     type = "l",
     cex = 0.25,
     col = alpha("black", 0.3),
     # axes = FALSE,
     ylab = "",
     yaxt="n",
     # xlab = "Years BP",
     bty = "l",
     xlim = c(1750, 3000),
     ylim = c(-17,-14))
axis(side = 2)
mtext(c, side = 2, line=2.5)
mtext("Years BP", side = 1, line = 0.8, outer = TRUE, cex = 0.85, font = 2)
lines(rollmean(d13c, 10, na.pad = TRUE) ~ bp, t.jack,
      col = "black")
}
pryr.savn %<a-% {plot(d15n ~ bp, t.sav,
                       type = "l",
                       cex = 0.5,
                       col = alpha("black", 0.3),
                       axes = FALSE,
                       xlab = "Years BP",
                       ylab = "",
                       yaxt="n",
                       # xlab = "",
                       ylim = c(6,11))
  axis(side = 2)
  mtext(n, side = 2, line=2.5)
  lines(rollmean(d15n, ma1, na.pad = TRUE) ~ bp, t.sav,
        col = "black")
}
pryr.stetn %<a-% {plot(d15n ~ linear.ad, t.stet,
     type = "l",
     cex = 0.5,
     col = alpha("black", 0.3),
     # axes = FALSE,
     ylab = "",
     yaxt="n", xlim = c(1800, 2005), ylim = c(7.5, 10.5))
axis(side = 2)
mtext(n, side = 2, line=2.5)
lines(rollmean(d15n, 5, na.pad = TRUE) ~ linear.ad, t.stet,
      col = "black")
}
pryr.stetc %<a-% {plot(d13c ~ linear.ad, t.stet,
                       type = "l",
                       cex = 0.5,
                       col = alpha("black", 0.3),
                       # axes = FALSE,
                       ylab = "",
                       yaxt="n",
                       xlab = "Years BP",
                       bty = "n", xlim = c(1800, 2005))
  axis(side = 2)
  mtext(c, side = 2, line=2.5)
  lines(rollmean(d13c, 20, fill = NA, align = "center") ~ linear.ad, t.stet,
        col = "black")
}

pryr.jackn
pryr.jackc
pryr.savn
pryr.stetn
pryr.stetc
segments(x0 = 256, x1 = 456, y0 = 7.5, y1 = 7.5, lty = "dashed")

# *****************************
# Paleoclimate reconstructions 
# in base R graphics
#
#
#
# *****************************

pryr.richey %<a-% {richey2007.melt %>%
    filter(variable == 'SST') %>%
    na.omit() %>%
    plot(value ~ yrBP, .,
         type = "l",
         bty = "l",
         col = alpha("black", 0.3))
  richey2007.melt %>%
    filter(variable == 'SST') %>%
    na.omit() %>%
  lines(rollmean(value, 5, na.pad = TRUE) ~ yrBP, ., col = "black")
  
}
pryr.schmidt %<a-% {schmidt.melt %>%
    filter(variable == 'sst') %>%
    na.omit() %>%
    plot(value ~ yrBP, .,
         type = "l",
         bty = "l",
         col = alpha("black", 0.3))
  schmidt.melt %>%
    filter(variable == 'sst') %>%
    na.omit() %>%
    lines(rollmean(value, 3, na.pad = TRUE) ~ yrBP, ., col = "black")
  
}
pryr.fisk %<a-% {fisk.melt %>%
    filter(variable == 'SST') %>%
    na.omit() %>%
    plot(value ~ yrBP, .,
         type = "l",
         bty = "l",
         col = alpha("black", 0.3))
  fisk.melt %>%
    filter(variable == 'SST') %>%
    na.omit() %>%
    lines(rollmean(value, 3, na.pad = TRUE) ~ yrBP, ., col = "black")
  
}
pryr.keigwin1 %<a-% {keigwin1.melt %>%
  filter(variable == 'carb.') %>% # carb. or d18Og.rub
  na.omit() %>%
  plot(value ~ yrBP, .,
       type = "o",
       bty = "l",
       col = alpha("black", 0.3),
       xlab = "Years BP")
# keigwin1.melt %>%
#   filter(variable == 'carb.') %>%
#   na.omit() %>%
#   lines(rollmean(value, 1, na.pad = TRUE) ~ yrBP, ., col = "black")
}
pryr.core1 %<a-% {core1.melt %>%
    filter(variable == 'd18Ow') %>% # Lots of proxies in this one: SST.Anand., d18oC, d18Osw, d18Ow
    na.omit() %>%
    plot(value ~ yrBP, .,
         type = "o",
         bty = "l",
         col = alpha("black", 0.3),
         xlab = "Years BP")
  core1.melt %>%
    filter(variable == 'd18Oc') %>%
    na.omit() %>%
    lines(rollmean(value, 1, na.pad = TRUE) ~ yrBP, ., col = "black")
}
pryr.bahamas %<a-% {bahamas.melt %>%
    filter(variable == 'SST.anomaly') %>%
    na.omit() %>%
    plot(value ~ yrBP, .,
         type = "l",
         bty = "l",
         col = alpha("black", 0.3),
         xlab = "Years BP")
  bahamas.melt %>%
    filter(variable == 'SST.anomaly') %>%
    na.omit() %>%
    lines(rollmean(value, 5, na.pad = TRUE) ~ yrBP, ., col = "black")
}
pryr.trouet %<a-% {trouet.melt %>%
    filter(variable == 'NAOms') %>%
    na.omit() %>%
    plot(value ~ yrBP, .,
         type = "l",
         bty = "l",
         col = alpha("black", 0.3),
         xlab = "Years BP")
  trouet.melt %>%
    filter(variable == 'NAOms') %>%
    na.omit() %>%
    lines(rollmean(value, 10, na.pad = TRUE) ~ yrBP, ., col = "black")
}


# ************
# Run them
#
#
# ************
pryr.richey
pryr.schmidt
pryr.fisk
pryr.keigwin1
pryr.core1
pryr.bahamas
pryr.trouet

# Test plots
{
plot(d15n ~ linear.ad, t.jack,
     type = "l",
     cex = 0.5,
     col = alpha("black", 0.0),
     # axes = FALSE,
     ylab = "",
     xlab = "Year CE",
     yaxt="l",
     ylim = c(7,11),
     # xlab = "Years BP",
     xlim = c(250,2005)
     )
axis(side = 2)
mtext(n, side = 2, line=2.5)
lines(rollmean(d15n, 1, na.pad = TRUE) ~ linear.ad, t.jack,
      col = "#1d91c0")
lines(rollmean(d15n, 3, na.pad = TRUE) ~ linear.ad, t.sav,
      col = "#225ea8")
lines(rollmean(d15n, 3, na.pad = TRUE) ~ linear.ad, t.stet,
      col = "#253494")
# lines(rollmean(d15n, 3, na.pad = TRUE) ~ linear.ad, t.jack4684,
#       col = "#081d58")

nfix <- 0
t.stet$nfix <- (1 - ((t.stet$d15n - nfix)/(nitrate - nfix)))

plot(nfix ~ bp, t.stet,
     type = "l",
     cex = 0.5,
     col = alpha("black", 0.3))
}

pdf('~/Documents/GitHub/rstudio/paleo_richey2007_GOM.pdf')
pryr.richey
dev.off()

pdf('~/Documents/GitHub/rstudio/paleo_schmidt2012_GOM.pdf')
pryr.schmidt
dev.off()

pdf('~/Documents/GitHub/rstudio/paleo_richey2009_GOM.pdf')
pryr.fisk
dev.off()

pdf('~/Documents/GitHub/rstudio/paleo_keigwin.pdf')
pryr.keigwin1
dev.off()

pdf('~/Documents/GitHub/rstudio/paleo_saenger2011_carolina_slope.pdf')
pryr.core1
dev.off()

pdf('~/Documents/GitHub/rstudio/paleo_saenger2009_bahamas_recon_sst.pdf')
pryr.bahamas
dev.off()

pdf('~/Documents/GitHub/rstudio/paleo_trouet_recon_NAOms.pdf')
pryr.trouet
dev.off()