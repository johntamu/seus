#' ------------------------------------------
#' Note: Make sure to have data loaded from
#' figures script file v2.
#' ------------------------------------------
#' 08/03/2019
#' 

# Required packages
library(maps)
library(ggplot2)
library(marmap)
library(oceanmap)
library(dplyr)
library(ggplot2)
library(lattice)
library(latticeExtra)
library(reshape2)
library(Hmisc)
library(grid)
library(zoo) 
library(forecast)

##################################
## Load df.bulk from timeseries ##
## script file                  ##
##################################
path1 <- '~/Documents/GitHub/data/schiff_bulk_years_09-04-2019.csv'
# path2 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff_bulk_years_08-14-2019.csv'
# path3 <- '/home/john/Desktop/data/schiff_bulk_years_08-14-2019.csv'

path <- path1

df.bulk <- read.csv(path, header = TRUE)
# View(df.bulk)

colnames(df.bulk)[names(df.bulk) == "distance..mm."] <- "distance" # Rename some columns for easier coding
colnames(df.bulk)[names(df.bulk) == "d15n.vs.air"] <- "d15n"
colnames(df.bulk)[names(df.bulk) == "d13c.vs.vpdb"] <- "d13c"

df.jack <- df.bulk %>% filter(coral.id == 'jack-4907-bc1-d3')
df.sav <- df.bulk %>% filter(coral.id == 'sav-4902-bc1-unk')
df.stet <- df.bulk %>% filter(coral.id == 'stet-4904-bc1-d2')
df.jack2 <- df.bulk %>% filter(coral.id == 'jack-4907-bc1-d1')
df.jack4684 <- df.bulk %>% filter(coral.id == 'jack-4684-bc-unk')
df.jack4686 <- df.bulk %>% filter(coral.id == 'jack-4686-bc-d1-t1')

x <- 'Year CE'
phe <- expression({delta}^15*"N"[" Phe"]*" (\u2030)")
n <- expression(delta^{15}*"N (\u2030)")
c <- expression(delta^{13}*"C (\u2030)")
eaa.neaa <- c('Phe', 'Thr', 'Ile', 'Leu', 'Val', 'Asx', 'Glx', 'Pro', 'Ala', 'Ser', 'Gly') # For Essential/Non-Essential ordering
tr.srcaa <- c('Glu', 'Asp', 'Ala', 'Ile', 'Leu', 'Pro', 'Val', 'Gly', 'Ser', 'Lys', 'Tyr', 'Phe', 'Thr') # For Trophic/Source AA ordering

#' --------------------------------
#' Summaries of any data frame you
#' are interested in
#' --------------------------------

summary(df.stet$d15n)




#' --------------------------------------------------------------------------------------
#' Linear age models (radiocarbon)
#' --------------------------------------------------------------------------------------

par(pty = "s", mfrow = c(2,2), mar=c(4,5,1,1))
plot(X14C.Age ~ Distance..um., r.jack,
     pch = 22, bg = '#ffeda0', col = "black", 
     ylab = expression({Delta}^14*"C Age (CRA)"), xlab = 'Distance from edge (um)')
# abline(mod1)
plot(X14C.Age ~ Distance..um., r.sav,
     pch = 23, bg = '#feb24c', col = "black", 
     ylab = expression({Delta}^14*"C Age (CRA)"), xlab = 'Distance from edge (um)')
# abline(mod2)
plot(X14C.Age ~ Distance..um., r.stet,
     pch = 21, bg = '#f03b20', col = "black",
     ylab = expression({Delta}^14*"C Age (CRA)"), xlab = 'Distance from edge (um)')
# abline(mod3)
r.jack$mm <- r.jack$Distance..um./1000
r.jack2$mm <- r.jack2$Distance..um./1000
r.stet$mm <- r.stet$Distance..um./1000
r.sav$mm <- r.sav$Distance..um./1000

par(pty = "s", mfrow = c(2,2), mar=c(4,5,1,1))
plot(mean ~ mm, r.jack,
     type = "o", bg = 'gray', col = "black", 
     ylab = "Years BP", xlab = 'Distance from edge (um)')
abline(jlm1.mm, lty = "dashed")
abline(jlm2.mm, lty = "dashed")
plot(X14C.Age ~ mm, whole,
     type = "o", bg = 'gray', col = "black", 
     ylab = "Years BP", xlab = 'Distance from edge (um)')
abline(lm.whole.mm, lty = "dashed")
plot(mean ~ mm, r.sav,
     type = "o", bg = 'gray', col = "black", 
     ylab = "Years BP", xlab = 'Distance from edge (um)')
abline(lm.sav.mm, lty = "dashed")
plot(mean ~ mm, r.stet,
     type = "o", bg = 'gray', col = "black",
     ylab = "Years BP", xlab = 'Distance from edge (um)')
abline(s1.mm, lty = "dashed")
abline(s2.mm)

#' Figure: 
#' --------------------------------------------------------------------------------------
#' Bomb spike age model for Jack-4684 BC1
#' --------------------------------------------------------------------------------------


#' Figure: 
#' --------------------------------------------------------------------------------------
#' Bulk data - overall
#' --------------------------------------------------------------------------------------

mod1 <- lm(d15n ~ d13c, df.jack)
mod2 <- lm(d15n ~ d13c, df.sav)
mod3 <- lm(d15n ~ d13c, df.stet)

summary(mod)
summary(mod2)
summary(mod3)

par(mfrow = c(1,3))
plot(d15n ~ d13c, df.jack,
     pch = 22, bg = '#ffeda0', col = "black", 
     ylab = n, xlab = c)
abline(mod1)
plot(d15n ~ d13c, df.sav,
     pch = 23, bg = '#feb24c', col = "black", 
     ylab = n, xlab = c)
abline(mod2)
plot(d15n ~ d13c, df.stet,
     pch = 21, bg = '#f03b20', col = "black",
     ylab = n, xlab = c)
abline(mod3)

#' Figure: 
#' --------------------------------------------------------------------------------------
#' Bulk data - testing reproducibility
#' --------------------------------------------------------------------------------------

# Jack-4907 bulk data vs distance from edge
# Link: https://stats.stackexchange.com/questions/4489/removing-borders-in-r-plots-for-achieving-tuftes-axis
# https://stackoverflow.com/questions/13239986/avoid-wasting-space-when-placing-multiple-aligned-plots-onto-one-page
p1.pryr %<a-% {plot(d15n ~ distance, df.jack,
     type = "l",
     cex = 0.5,
     col = alpha("black", 0.3),
     axes = FALSE,
     ylab = "",
     yaxt="n")
axis(side = 2)
mtext(n, side = 2, line=2.5)
lines(d15n ~ distance, df.jack2,
      col = alpha("black", 0.75), lwd = 1.5)
}
par(pty = "s", mfrow=c(2,1), oma = c(2, 2, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(1, 1, 0, 0)+0.1, # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0))   # axis label at 2 rows distance, tick labels at 1 row
plot(d15n ~ distance, df.jack,
     type = "l",
     cex = 0.5,
     col = alpha("black", 0.3),
     axes = FALSE,
     ylab = "",
     yaxt="n",
     xlab = "",
     ylim = c(6,11))
axis(side = 2)
mtext(n, side = 2, line=2.5)
lines(d15n ~ distance, df.jack2,
      col = alpha("black", 0.75), lwd = 1.5)

plot(d13c ~ distance, df.jack,
     xlab = "Distance (mm)",
     ylab = "",
     type = "l",
     cex = 0.5,
     col = alpha("black", 0.3),
     # axes = FALSE,
     bty = "n",
     yaxt = "n",
     ylim = c(-18,-14))
axis(side = 4)
mtext(c, side=4, line=2.5)
lines(d13c ~ distance, df.jack2,
      col = alpha("black", 0.75), lwd = 1.5)

par(pty = "s")
plot(d15n ~ linear.ad, df.jack,
     xlab = x,
     ylab = n,
     type = "l",
     cex = 0.5,
     col = alpha("black", 0.0))
lines(forecast::ma(df.jack$d15n, order = 3, centre = TRUE) ~ linear.ad, df.jack,
      col = alpha("black", 0.4), lwd = 1.5)
lines(forecast::ma(df.jack2$d15n, order = 3, centre = TRUE) ~ linear.ad, df.jack2,
      col = "#1f78b4", lwd = 1.5)
points(forecast::ma(df.jack2$d15n, order = 3, centre = TRUE) ~ linear.ad, df.jack2,
      col = "#1f78b4", lwd = 1.5)

par(pty = "s")
plot(d13c ~ linear.ad, df.jack,
     xlab = x,
     ylab = n,
     type = "l",
     cex = 0.5,
     col = alpha("black", 0.0))
lines(forecast::ma(df.jack$d13c, order = 3, centre = TRUE) ~ linear.ad, df.jack,
      col = alpha("black", 0.4), lwd = 1.5)
lines(forecast::ma(df.jack2$d13c, order = 2, centre = TRUE) ~ linear.ad, df.jack2,
      col = "#1f78b4", lwd = 1.5)
points(forecast::ma(df.jack2$d13c, order = 2, centre = TRUE) ~ linear.ad, df.jack2,
      col = "#1f78b4", lwd = 1.5)

# For now, use the figure already generated and copy paste that into the file

#' Figure: Recent history
#' --------------------------------------------------------------------------------------
#' Bulk data - Stet4904 and Jack4684 d15N
#' --------------------------------------------------------------------------------------

par(pty = "s")
plot(d15n ~ linear.ad, df.stet,
     xlab = x,
     ylab = n,
     type = "o",
     cex = 0.5,
     xlim = c(1500,2005),
     col = alpha("black", 0.9))
points(forecast::ma(df.stet$d15n, order = 1, centre = TRUE) ~ linear.ad, df.stet,
      col = "#33a02c", lwd = 1.5)
lines(forecast::ma(df.jack4684$d15n, order = 3, centre = TRUE) ~ linear.ad, df.jack4684,
      col = "#1f78b4", lwd = 1.5)

s1 <- ggplot(df.stet, aes(x=linear.ad, y=d15n))
s1 + geom_line(color = alpha("#009E73", 0.01)) +
  coord_fixed(ratio = 160) +
  geom_line(data=df.stet, aes(x=linear.ad, y=rollmean(d15n, 1, na.pad = TRUE)), color = "#009E73", size = 0.75) +
  geom_line(data=df.jack4684, aes(x=linear.ad, y=rollmean(d15n, 1, na.pad = TRUE)), color = "#D55E00", size = 0.75) +
  # theme_bw() +
  # theme_classic() +
  # geom_point() +
  # geom_point(shape = 21, size = 3, color = "black", fill = "red") +
  # theme_classic() +
  # xlim(1400, 2010) +
  # scale_x_continuous(breaks = seq(1000,2010,by=100), limits = c(1000, 2010)) +
  # scale_x_continuous(breaks = seq(1650, 1800, by = 25), limits = c(1650, 1800)) +
  theme(axis.text.y   = element_text(size=15, color = "black"),
        axis.text.x   = element_text(size=15, color = "black"),
        axis.title.y  = element_text(size=15),
        axis.title.x  = element_text(size=15),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        legend.position = "topright") +
  annotate("text",label="Stetson-4904 BC1", x = 1200, y = 10, color = "#009E73", size = 5, fontface = "bold", hjust = 0) +
  annotate("text",label="Jacksonville-4907 BC1", x = 1200, y = 9.75, color = "#D55E00", size = 5, fontface = "bold", hjust = 0) +
  xlim(1200, 2010) +
  xlab(x) +
  ylab(n)
ggsave('~/Documents/GitHub/rstudio/ggplot_figs/bulk_n_jack4684_stet4904.png', plot = last_plot(), dpi = 300, width = 6, height = 5.5)
dev.off()


#' Figure: Recent history, bulk
#' --------------------------------------------------------------------------------------
#' Bulk data - d15N and d13C plotted againt distance from the edge (mm)
#' --------------------------------------------------------------------------------------

# Nitrogen
par(pty = "s", mfrow = c(2,2), mar=c(4,5,1,1))
plot(d15n ~ distance, df.stet,
     xlab = 'Distance from edge (mm)',
     ylab = n,
     type = "l",
     cex = 0.75,
     col = alpha("#0072B2", 0.99))
# text(2,7, labels = "A", cex = 2.25)
# legend(7,2, box.lty = 0, legend = "Stetson-4904-BC1 Disk 1", bg = NULL)
plot(d15n ~ distance, df.jack,
     xlab = 'Distance from edge (mm)',
     ylab = n,
     type = "l",
     cex = 0.75,
     col = alpha("#0072B2", 0.99))
# text(2,7, labels = "A", cex = 2.25)
# legend("bottomleft", box.lty = 0, legend = "Jacksonville-4907-BC1 Disk 3", bg = NULL)
plot(d15n ~ distance, df.jack4684,
     xlab = 'Distance from edge (mm)',
     ylab = n,
     type = "l",
     cex = 0.75,
     col = alpha("#0072B2", 0.99))
# text(2,7, labels = "A", cex = 2.25)
# legend("bottomleft", box.lty = 0, legend = "Jacksonville-4684-BC1 Disk 1", bg = NULL)
plot(d15n ~ distance, df.sav,
     xlab = 'Distance from edge (mm)',
     ylab = n,
     type = "l",
     cex = 0.75,
     col = alpha("#0072B2", 0.99))
# text(2,7, labels = "A", cex = 2.25)
# legend("bottomleft", box.lty = 0, legend = "Savannah Banks-BC1 Base 1", bg = NULL)

# Carbon
par(pty = "s", mfrow = c(2,2), mar=c(4,5,1,1))
plot(d13c ~ distance, df.stet,
     xlab = 'Distance from edge (mm)',
     ylab = c,
     type = "l",
     cex = 0.75,
     col = alpha("#0072B2", 0.99))
# text(2,7, labels = "A", cex = 2.25)
# legend(7,2, box.lty = 0, legend = "Stetson-4904-BC1 Disk 1", bg = NULL)
plot(d13c ~ distance, df.jack,
     xlab = 'Distance from edge (mm)',
     ylab = c,
     type = "l",
     cex = 0.75,
     col = alpha("#0072B2", 0.99))
# text(2,7, labels = "A", cex = 2.25)
# legend("bottomleft", box.lty = 0, legend = "Jacksonville-4907-BC1 Disk 3", bg = NULL)
df.jack4684 %>%
  filter(d13c )
plot(d13c ~ distance, df.jack4684,
     xlab = 'Distance from edge (mm)',
     ylab = c,
     type = "l",
     cex = 0.75,
     col = alpha("#0072B2", 0.99))
# text(2,7, labels = "A", cex = 2.25)
# legend("bottomleft", box.lty = 0, legend = "Jacksonville-4684-BC1 Disk 1", bg = NULL)
plot(d13c ~ distance, df.sav,
     xlab = 'Distance from edge (mm)',
     ylab = c,
     type = "l",
     cex = 0.75,
     col = alpha("#0072B2", 0.99))
# text(2,7, labels = "A", cex = 2.25)
# legend("bottomleft", box.lty = 0, legend = "Savannah Banks-BC1 Base 1", bg = NULL)


#' Figure: Recent history, bulk
#' --------------------------------------------------------------------------------------
#' Bulk data - Stet4904
#' --------------------------------------------------------------------------------------

plot(d13c ~ linear.ad, df.stet,
     xlab = x,
     ylab = c,
     type = "l",
     cex = 0.5,
     # xlim = c(1500,2005),
     col = alpha("black", 0.3))
abline(v = 950, col = "black", lty = 'dashed')
abline(v = 1250, col = "black", lty = 'dashed')
lines(forecast::ma(df.stet$d13c, order = 3, centre = TRUE) ~ linear.ad, df.stet,
      col = "#33a02c", lwd = 1.5)

plot(d15n ~ linear.ad, df.stet,
     xlab = x,
     ylab = c,
     type = "l",
     cex = 0.5,
     col = alpha("black", 0.3))
abline(v = 950, col = "black", lty = 'dashed')
abline(v = 1250, col = "black", lty = 'dashed')
lines(forecast::ma(df.stet$d15n, order = 3, centre = TRUE) ~ linear.ad, df.stet,
      col = "#33a02c", lwd = 1.5)

obj2 <- xyplot(forecast::ma(df.stet$d15n, order = 3, centre = TRUE) ~ linear.ad, data = df.stet, type = "l", lwd = 1.5,
               ylab = n, xlab = x, col = '#081d58', ylim = c(5, 10))
obj1 <- xyplot(forecast::ma(df.stet$d13c, order = 3, centre = TRUE) ~ linear.ad, data = df.stet, type = "l", lwd = 1.5,
               xlab = x, ylab = c, col = '#225ea8', ylim = c(-18, -10))
doubleYScale(obj2, obj1, add.ylab2 = TRUE, style1= NULL, style2 = NULL)

splot1 <- df.stet %>%
  dplyr::select(linear.ad, d15n) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d15n), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.45) +
  geom_line(aes(y=rollmean(d15n, 3, na.pad = TRUE)), size = 0.45, color = '#225ea8') +
  
  geom_vline(xintercept = 950, lty = 'dotted') +
  geom_vline(xintercept = 1250, lty = 'dotted') +
  geom_vline(xintercept = 1550, lty = 'dotted') +
  geom_vline(xintercept = 1850, lty = 'dotted') +
  
  annotate("text", x = 1100, y = 7, label = 'Medieval Warming', size = 3.5) +
  annotate("text", x = 1700, y = 7, label = "Little Ice Age", size = 3.5) +
  # annotate("segment", x = 625, xend = 975, y = 9.75, yend = 9.75, color = "black") +
  # annotate("text", x = 800, y = 9.85, label = "error = 350 yrs", size = 3) +
  # annotate("pointrange", x = 800, y = 9.75, ymin = 9.75, ymax = 9.75) +
  
  ylab(n) +
  xlab(NULL) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), position = "top") +
  theme_bw() +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_text(size=10, color = "black"),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        # axis.line.x = element_blank(),
        axis.ticks.x = element_line(size = 0.25, color = "black"),
        axis.ticks.y = element_line(size = 0.25, color = "black"),
        panel.background = element_rect(size = 0.0, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
        # plot.background = element_rect(size =0.75))

splot2 <- df.stet %>%
  dplyr::select(linear.ad, d13c) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d13c), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.45) +
  geom_line(aes(y=rollmean(d13c, 3, na.pad = TRUE)), size = 0.45, color = '#081d58') +
  
  geom_vline(xintercept = 950, lty = 'dotted') +
  geom_vline(xintercept = 1250, lty = 'dotted') +
  geom_vline(xintercept = 1550, lty = 'dotted') +
  geom_vline(xintercept = 1850, lty = 'dotted') +
  
  ylab(c) +
  theme_bw() +
  xlab(x) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_text(size=10, color = "black"),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        axis.ticks.x = element_line(size = 0.25),
        axis.ticks.y = element_line(size = 0.25),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

grid.newpage()
grid.draw(rbind(ggplotGrob(splot1), ggplotGrob(splot2), size = "last"))
dev.copy(pdf,"~/Documents/GitHub/rstudio/ggplot_figs/bulk_n_c_stet4904.pdf", width = 6, height = 6)
dev.off()

#' Figure: Recent history, bulk
#' --------------------------------------------------------------------------------------
#' Stet4904 d13C vs d15N XY plot during certain periods
#' --------------------------------------------------------------------------------------

df.stet %>%
  filter(linear.ad < 1100) -> mwp

plot(d13c ~ d15n, mwp)
abline(lm(d13c ~ d15n, mwp))
summary(lm(d13c ~ d15n, mwp))

#' Figure: Recent history, bulk
#' --------------------------------------------------------------------------------------
#' Bulk data - Jack4686 vs distance
#' --------------------------------------------------------------------------------------

plot(d15n ~ distance, df.jack4686,
     xlab = "Distance from edge (mm)",
     ylab = n,
     type = "o",
     cex = 0.5)
# xlim = c(1500,2005),
# col = alpha("black", 0.3))

#' Figure: Ancient bulk values
#' --------------------------------------------------------------------------------------
#' Savannah-4902 by itself
#' --------------------------------------------------------------------------------------

savplot1 <- df.sav %>%
  dplyr::select(linear.ad, d15n) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d15n), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.75) +
  geom_line(aes(y=rollmean(d15n, 3, na.pad = TRUE)), size = 0.75, color = '#02818a') +
  # geom_point(color = '#31a354', shape = 21) +
  
  geom_vline(xintercept = 950, lty = 'dashed') +
  geom_vline(xintercept = 1250, lty = 'dashed') +
  geom_vline(xintercept = -250, lty = "dotted") +
  geom_vline(xintercept = 400, lty = "dotted") +
  
  annotate("text", x = 1100, y = 7, label = 'Medieval Warming', size = 2.5) +
  annotate("text", x = -600, y = 7, label = "Iron Age Cold Epoch", size = 4) +
  annotate("text", x = 75, y = 7, label = "Roman Warm Period", size = 3) +
  
  ylab(n) +
  xlab(NULL) +
  # xlim(-1200, 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(-300,1300)) +
  theme_classic() +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_blank(),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

savplot2 <- df.sav %>%
  dplyr::select(linear.ad, d13c) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d13c), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.75) +
  geom_line(aes(y=rollmean(d13c, 3, na.pad = TRUE)), size = 0.75, color = '#3690c0') +
  # geom_point(color = '#addd8e', shape = 21) +
  
  geom_vline(xintercept = 950, lty = 'dashed') +
  geom_vline(xintercept = 1250, lty = 'dashed') +
  geom_vline(xintercept = -250, lty = "dotted") +
  geom_vline(xintercept = 400, lty = "dotted") +
  
  ylab(c) +
  theme_classic() +
  xlab(x) +
  # xlim(-1200, 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(-300,1300)) +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_text(size=10, color = "black"),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

grid.newpage()
grid.draw(rbind(ggplotGrob(savplot1), ggplotGrob(savplot2), size = "last"))




#' Figure: Ancient bulk values
#' --------------------------------------------------------------------------------------
#' Jack-4907 by itself
#' --------------------------------------------------------------------------------------

plot(d13c ~ linear.ad, df.jack,
     xlab = x,
     ylab = c,
     type = "l",
     cex = 0.5,
     col = alpha("black", 0.3))
abline(v = -250, col = "black", lty = "dashed")
abline(v = 400, col = "black", lty = "dashed")
abline(v = -900, col = alpha("black", 0.75), lty = "longdash")
abline(v = -300, col = alpha("black", 0.75), lty = "longdash")
lines(forecast::ma(df.jack$d13c, order = 3, centre = TRUE) ~ linear.ad, df.jack,
      col = "#4575b4", lwd = 1.5)

plot(d15n ~ linear.ad, df.jack,
     xlab = x,
     ylab = n,
     type = "l",
     cex = 0.5,
     col = alpha("black", 0.3))
abline(v = -250, col = "black", lty = "dashed")
abline(v = 400, col = "black", lty = "dashed")
abline(v = -900, col = alpha("black", 0.75), lty = "longdash")
abline(v = -300, col = alpha("black", 0.75), lty = "longdash")
lines(forecast::ma(df.jack$d15n, order = 3, centre = TRUE) ~ linear.ad, df.jack,
      col = "#4575b4", lwd = 1.5)

obj2 <- xyplot(forecast::ma(df.jack$d15n, order = 1, centre = TRUE) ~ distance, data = df.jack, type = "l", lwd = 1.5,
               ylab = n, xlab = x, col = 'red')
obj1 <- xyplot(forecast::ma(df.jack$d13c, order = 1, centre = TRUE) ~ distance, data = df.jack, type = "l", lwd = 1.5,
               xlab = x, ylab = c, col = 'blue')
doubleYScale(obj2, obj1, add.ylab2 = TRUE, style1= NULL, style2 = NULL)

t.jack <- df.jack
t.jack$bp <- 1950 - t.jack$linear.ad

jplot1 <- t.jack %>%
  dplyr::select(bp, d15n) %>%
  na.omit() %>%
  ggplot(aes(x = bp, y = d15n), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.5) +
  geom_line(aes(y=rollmean(d15n, 3, na.pad = TRUE)), size = 0.5, color = '#b30000') +
  
  geom_vline(xintercept = 1000, lty = 'dashed') +
  geom_vline(xintercept = 700, lty = 'dashed') +
  geom_vline(xintercept = 2850, lty = 'longdash') +
  geom_vline(xintercept = 2250, lty = 'longdash') +
  geom_vline(xintercept = 2200, lty = "dotted") +
  geom_vline(xintercept = 1550, lty = "dotted") +
  
  annotate("text", x = 825, y = 7, label = 'MCA', size = 3.4) +
  annotate("text", x = 2550, y = 7, label = "Iron Age Cold Epoch", size = 3.4) +
  annotate("text", x = 1875, y = 7, label = "Roman Warm Period", size = 3.4) +
  
  ylab(n) +
  xlab(NULL) +
  # xlim(-1200, 0) +
  # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_classic() +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_blank(),
        axis.line.x = element_blank(),
        # axis.ticks.x = element_line(size = 0.25, color = "black"),
        axis.ticks.y = element_line(size = 0.25, color = "black"),
        axis.line = element_line(size = 0.25, color = "black"),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# plot.background = element_rect(size =0.75))

jplot2 <- t.jack %>%
  dplyr::select(bp, d13c) %>%
  na.omit() %>%
  ggplot(aes(x = bp, y = d13c), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.5) +
  geom_line(aes(y=rollmean(d13c, 3, na.pad = TRUE)), size = 0.5, color = '#d7301f') +
  
  geom_vline(xintercept = 1000, lty = 'dashed') +
  geom_vline(xintercept = 700, lty = 'dashed') +
  geom_vline(xintercept = 2850, lty = 'longdash') +
  geom_vline(xintercept = 2250, lty = 'longdash') +
  geom_vline(xintercept = 2200, lty = "dotted") +
  geom_vline(xintercept = 1550, lty = "dotted") +

  ylab(c) +
  theme_classic() +
  xlab("Years BP") +
  # xlim(-1200, 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_text(size=10, color = "black"),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        axis.line = element_line(size = 0.25, color = "black"), 
        # axis.ticks.x = element_line(size = 0.25, color = "black"),
        axis.ticks.y = element_line(size = 0.25, color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

grid.newpage()
grid.arrange(jplot1, jplot2, ncol = 1)
g <- arrangeGrob(jplot1, jplot2, ncol = 1)
ggsave('~/Documents/GitHub/rstudio/ggplot_figs/bulk_n_c_jack4907.png', g, dpi = 300) # This works! 9-4-2019
# grid.draw(rbind(ggplotGrob(jplot1), ggplotGrob(jplot2), ggplotGrob(paleoplot1), size = "last"))

jplot1 <- t.jack %>%
  dplyr::select(linear.ad, d15n) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d15n), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.65) +
  geom_line(aes(y=rollmean(d15n, 3, na.pad = TRUE)), size = 0.65, color = '#cc4c02') +
  
  geom_vline(xintercept = 950, lty = 'dashed') +
  geom_vline(xintercept = 1250, lty = 'dashed') +
  geom_vline(xintercept = -900, lty = 'longdash') +
  geom_vline(xintercept = -300, lty = 'longdash') +
  geom_vline(xintercept = -250, lty = "dotted") +
  geom_vline(xintercept = 400, lty = "dotted") +

  annotate("text", x = 1100, y = 7, label = 'Medieval Warming', size = 2) +
  annotate("text", x = -600, y = 7, label = "Iron Age Cold Epoch", size = 4) +
  annotate("text", x = 75, y = 7, label = "Roman Warm Period", size = 2) +
  
  ylab(n) +
  xlab(NULL) +
  # xlim(-1200, 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_classic() +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_blank(),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

jplot2 <- t.jack %>%
  dplyr::select(linear.ad, d13c) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d13c), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.75) +
  geom_line(aes(y=rollmean(d13c, 3, na.pad = TRUE)), size = 0.85, color = '#fe9929') +
  
  geom_vline(xintercept = 950, lty = 'dashed') +
  geom_vline(xintercept = 1250, lty = 'dashed') +
  geom_vline(xintercept = -900, lty = 'longdash') +
  geom_vline(xintercept = -300, lty = 'longdash') +
  geom_vline(xintercept = -250, lty = "dotted") +
  geom_vline(xintercept = 400, lty = "dotted") +

  ylab(c) +
  theme_classic() +
  xlab("Year CE") +
  # xlim(-1200, 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_text(size=10, color = "black"),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

grid.newpage()
grid.draw(rbind(ggplotGrob(jplot1), ggplotGrob(jplot2), size = "last"))

#' Figure: Ancient bulk values
#' --------------------------------------------------------------------------------------
#' Jack-4907 vs Sav-4902
#' --------------------------------------------------------------------------------------

plot(d15n ~ linear.ad, df.jack, # With Base R
     xlab = x,
     ylab = n,
     type = "l",
     cex = 0.5,
     xlim = c(-1100,250),
     col = alpha("black", 0.3))
abline(v = 950, col = "black", lty = 'dashed')
abline(v = 1250, col = "black", lty = 'dashed')
abline(v = 400, col = alpha("black", 0.75), lty = "longdash")
abline(v = -900, col = "black", lty = 'dashed')
abline(v = -300, col = "black", lty = 'dashed')
lines(forecast::ma(df.jack$d15n, order = 15, centre = TRUE) ~ linear.ad, df.jack,
      col = "#4575b4", lwd = 1.5)
lines(forecast::ma(df.sav$d15n, order = 3, centre = TRUE) ~ linear.ad, df.sav,
      col = "#f46d43", lwd = 1.5)
lines(forecast::ma(df.stet$d15n, order = 5, centre = TRUE) ~ linear.ad, df.stet,
      col = "#f46d43", lwd = 1.5)

p1 <- ggplot() + # With ggplot2
  geom_line(data=df.jack, aes(x = linear.ad, y = d15n), color = "gray", alpha = 0.0, size = 0.5) +
  geom_line(data=df.jack, aes(x = linear.ad,
                              y = rollmean(d15n, 1.5, na.pad = TRUE)), color = "#1f78b4", alpha = 0.99, size = 0.75) +
  # geom_line(data=df.sav, aes(x = linear.ad, y = d15n), color = "gray", alpha = 0.0, size = 0.5) +
  # geom_line(data=df.sav, aes(x = linear.ad,
  #                            y = rollmean(d15n, 12, na.pad = TRUE)), color = "#33a02c", alpha = 0.99, size = 0.75) +
  geom_line(data=df.stet, aes(x = linear.ad, y = d15n), color = "gray", alpha = 0.0, size = 0.5) +
  geom_line(data=df.stet, aes(x = linear.ad,
                             y = rollmean(d15n, 15, na.pad = TRUE)), color = "#ff7f00", alpha = 0.99, size = 0.75) +
  
  annotate("text", x = 1100, y = 7, label = 'Medieval Warming', size = 3) +
  annotate("text", x = 75, y = 7, label = "Roman Warm Period", size = 3) +
  
  geom_vline(xintercept = 950, lty = 'dotted') +
  geom_vline(xintercept = 1250, lty = 'dotted') +
  # geom_vline(xintercept = -900, lty = 'longdash') +
  # geom_vline(xintercept = -300, lty = 'longdash') +
  geom_vline(xintercept = -250, lty = "dotted") +
  geom_vline(xintercept = 400, lty = "dotted") +
  
  ylab(n) +
  theme_classic() +
  xlab(x) +
  # xlim(-1200, 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(-300,1400)) +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_text(size=10, color = "black"),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p1

plot(d13c ~ linear.ad, df.sav,
     xlab = x,
     ylab = c,
     type = "l",
     cex = 0.5,
     ylim = c(-17,-14.50),
     xlim = c(600,1250),
     col = alpha("black", 0.0))
abline(v = 950, col = "black", lty = 'dashed')
abline(v = 1250, col = "black", lty = 'dashed')
abline(v = 400, col = alpha("black", 0.75), lty = "longdash")
lines(forecast::ma(df.jack$d13c, order = 1, centre = TRUE) ~ linear.ad, df.jack,
      col = alpha("#41b6c4", 0.99), lwd = 1.5)
lines(forecast::ma(df.sav$d13c, order = 1, centre = TRUE) ~ linear.ad, df.sav,
      col = alpha("#c7e9b4", 1.0), lwd = 1.5)
lines(forecast::ma(df.stet$d13c, order = 1, centre = TRUE) ~ linear.ad, df.stet,
      col = "#253494", lwd = 1.5)

plot(d15n ~ linear.ad, df.sav,
     xlab = x,
     ylab = n,
     type = "o",
     cex = 0.5,
     xlim = c(500,1500),
     col = alpha("black", 0.0))
abline(v = 950, col = "black", lty = 'dashed')
abline(v = 1250, col = "black", lty = 'dashed')
lines(d15n ~ linear.ad, df.jack,
      col = "#4575b4", lwd = 1.5, type = "o")
lines(d15n ~ linear.ad, df.sav,
      col = "#f46d43", lwd = 1.5, type = "l")

#' -----------------------------------------------------
#' Compound-specific figures below
#' 
#' 
#' 
#' 
#' 
#' -----------------------------------------------------

################################
##                            ##
##                            ##
## N-CSIAA Data visualization ##
##                            ##
##                            ##
################################
ndata <- read.csv("~/Google Drive/projects/rproj/seus/cleaned_ndata.csv") # Contains SEUS black coral data, GOM black corals, and some POM data from elsewhere

######################
## Overall AA plots ##
######################
# Helpful link: https://stackoverflow.com/questions/21982987/mean-per-group-in-a-data-frame

bwtheme <- standard.theme("pdf", color=FALSE)

myColours <- brewer.pal(9,'Greys')

my.settings <- list(
  superpose.symbol=list(fill=myColours[2:5],col = "black", border="transparent", pch=c(21,23,22,24,25,8)),
  strip.background=list(col=myColours[6]),
  strip.border=list(col='black'))

ndata %>%
  melt(., "Type") %>%
  filter(variable %in% tr.srcaa) %>% # Success! 03-04-2019
  xyplot(value ~ variable,
         .,
         panel = function(x,  y, ...){
           panel.xyplot(x, y, ...)
           panel.abline(h = 9, lty = 1)
           panel.abline(h = 9.25, lty = 2)
           panel.abline(h = 8.75, lty = 2)
         },
         cex=1.5,
         ylim=c(-20, 40), # ylim and xlim are "first class" parameters and don't need to be in scales=list()
         group = Type,
         xlab = NULL,
         ylab = expression({delta}^15*"N (\u2030)"),
         auto.key=list(columns=2, cex = 0.75),
         par.settings = my.settings)


ndata.melt <- melt(ndata[1:45,], "Sample.ID2")
ndata.melt$variable1 <- factor(ndata.melt$variable, 
                               levels=tr.srcaa)

bwtheme <- standard.theme("pdf", color = FALSE) # Figure for black and white overall plot, using myColors above
xyplot(value ~ variable1,
       ndata.melt,
       panel = function(x,  y, ...){
         panel.xyplot(x, y, ...)
         panel.abline(h = 9, lty = 1)
         panel.abline(h = 9.25, lty = 2)
         panel.abline(h = 8.75, lty = 2)
       },
       cex=1.5,
       ylim=c(-20, 40), # ylim and xlim are "first class" parameters and don't need to be in scales=list()
       group = Sample.ID2,
       xlab = NULL,
       ylab = expression({delta}^15*"N (\u2030)"),
       auto.key=list(columns=2, cex = 0.75),
       par.settings = my.settings)

######################
## Trophic Dynamics ## 
######################
ndata %>%
  filter(Type == "Black Coral - GOM" | Type == "Black Coral - SEUS") -> bc.ndata

p <- ggplot(bc.ndata, aes(x=Sum.V, y=TP, shape = Region, fill = Region), color = "black")
p + geom_point(size = 3) +
  facet_wrap(~Region, ncol = 1) +
  scale_shape_manual(values = c(22, 23)) +
  scale_color_manual(values = c("#D55E00", "#E69F00")) +
  ylab("Trophic Position (Glu - Phe)") +
  xlab(expression(paste(Sigma,"V"))) +
  theme_bw() +
  theme(axis.text=element_text(size=11, color = "black"))


########################
## SumV through time  ##
########################
ndata %>%
  filter(Region == "SEUS") %>%
  droplevels(.) %>%
  xyplot(Sum.V ~ Year.CE,
         data = .,
         group = Sample.ID2,
         cex = 1.5,
         # xlim = c(1500,2010),
         par.settings = my.settings,
         auto.key = list(columns=c(2), cex = 0.95),
         xlab = x,
         ylab = expression(paste(Sigma,"V")))

########################
## TP through time    ##
########################
ndata %>%
  # filter(Region == "SEUS") %>%
  # droplevels(.) %>%
  xyplot(TP ~ Year.CE,
         data = .,
         group = Sample.ID2,
         cex = 1.5,
         # xlim = c(750, 1500),
         par.settings = my.settings,
         auto.key = list(columns=c(2), cex = 0.95),
         xlab = x,
         ylab = "Trophic Position (Glu - Phe)")

######################
## Phe through time ##
######################

ndata %>% # Showing Phe through time, but only for select specimens (modify code as needed)
  filter(Phe > 1) %>%
  filter(Phe < 13) %>%
  filter(Sample.ID2 == "Jacksonville-4907" | Sample.ID2 == "Savannah Banks-4902" | Sample.ID2 == "Jacksonville-4684") -> t.ndata

par(mfrow=c(1,2))
xyplot(Phe ~ Year.CE,
       data = t.ndata,
       groups = Sample.ID2,
       # xlim =c(1500,1900),
       xlab = x,
       # xlim = c(-250,1500),
       ylab = phe,
       cex = 3,
       ylim = c(4, 14))


################################
##                            ##
##                            ##
## C-CSIAA Data visualization ##
##                            ##
##                            ##
################################
seus_carbon <- read.csv("C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff c-csiaa stetson_cleaned.csv")
seus_carbon <- read.csv("~/Google Drive/projects/rproj/seus/data/schiff c-csiaa stetson_cleaned.csv")

######################
## Overall AA Plots ##
######################
# Helpful link: https://stackoverflow.com/questions/21982987/mean-per-group-in-a-data-frame