#'
#' John Schiff
#' Created: 02/27/2019
#' 

library(dplyr)
library(ggplot2)
library(lattice)
library(latticeExtra)
library(reshape2)
library(Hmisc)
library(grid)

##################################
## Load df.bulk from timeseries ##
## script file                  ##
##################################
path1 <- '~/Documents/GitHub/rstudio/data/schiff_bulk_years_08-04-2019.csv'
path2 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff_bulk_years_08-04-2019.csv'
path3 <- '/home/john/Desktop/data/schiff_bulk_years_08-04-2019.csv'
  
path <- read.csv(path1)

df.bulk <- path
View(df.bulk)

colnames(df.bulk)[names(df.bulk) == "distance..mm."] <- "distance" # Rename some columns for easier coding
colnames(df.bulk)[names(df.bulk) == "d15n.vs.air"] <- "d15n"
colnames(df.bulk)[names(df.bulk) == "d13c.vs.vpdb"] <- "d13c"

# write.csv(df.bulk, 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff_bulk_years_03-04-2019.csv')

df.jack <- df.bulk %>% filter(coral.id == 'jack-4907-bc1-d3')
df.sav <- df.bulk %>% filter(coral.id == 'sav-4902-bc1-unk')
df.stet <- df.bulk %>% filter(coral.id == 'stet-4904-bc1-d2')
df.jack2 <- df.bulk %>% filter(coral.id == 'jack-4907-bc1-d1')
df.jack4684 <- df.bulk %>% filter(coral.id == 'jack-4684-bc-unk')
df.jack4686 <- df.bulk %>% filter(coral.id == 'jack-4686-bc-d1-t1')

#' --------------------------------------------------------------- 07/12/2019
#' Adding the THIRD growth rate for the Stetson coral.
#' This accounts for 1) the split growth rates
#' and also 2) the apparently faster growth rate identified in 
#' the iodine data.
#'
#' I SHOULD be doing this in the age model/time series scripts, but I would rather not
#' reload everything right now.
#' ---------------------------------------------------------------

# df.stet$linear.ad3 <- ... combined ages from the two growth rates
df.stet$linear.ad4 <- stet.linear.ad4



#'
#' Commonly used strings for axis labels
#' 
x <- 'Year CE'
phe <- expression({delta}^15*"N"[" Phe"]*" (\u2030)")
n <- expression(delta^{15}*"N (\u2030)")
c <- expression(delta^{13}*"C (\u2030)")
eaa.neaa <- c('Phe', 'Thr', 'Ile', 'Leu', 'Val', 'Asx', 'Glx', 'Pro', 'Ala', 'Ser', 'Gly') # For Essential/Non-Essential ordering
tr.srcaa <- c('Glu', 'Asp', 'Ala', 'Ile', 'Leu', 'Pro', 'Val', 'Gly', 'Ser', 'Lys', 'Tyr', 'Phe', 'Thr') # For Trophic/Source AA ordering

# Below is a general function to generate a continuous color palette
# e.g., I use it when comparing bulk d13c vs. bulk d15n and coloring by year
rbPal <- colorRampPalette(c('red','blue')) # '#fed976, #feb24c, #fd8d3c, #fc4e2a, #e31a1c, #bd0026, #800026'
df.stet$Col <- rbPal(10)[as.numeric(cut(df.stet$linear.ad, breaks = 10))]
t.df$Col <- rbPal(10)[as.numeric(cut(t.df$linear.ad, breaks = 10))]

##################################
# Some lattice par.settings code #
##################################
bwtheme <- standard.theme("pdf", color=FALSE)

myColours <- brewer.pal(9,'Greys')

my.settings <- list(
  superpose.symbol=list(fill=myColours[2:5],col = "black", border="transparent", pch=c(21,23,22,24,25,8)),
  strip.background=list(col=myColours[6]),
  strip.border=list(col='black'))

myColors <- brewer.pal(9, "YlGnBu")
myColors3 <- brewer.pal(9, "BuGn")

my.settings2 <- list(
  superpose.symbol=list(pch = c(21,22,23,24,25,8),
                        col = "black",
                        fill=myColors[1:8]),
  strip.background=list(col=myColors[8]),
  strip.border=list(col='black'))

my.settings3 <- list(
  superpose.symbol=list(pch = c(21,22,23,24,25,8),
                        col = "black",
                        fill=myColors3[1:8]),
  strip.background=list(col=myColors3[8]),
  strip.border=list(col='black'))

bw_theme <- trellis.par.get()
bw_theme$box.dot$pch <- "|"
bw_theme$box.rectangle$col <- "black"
bw_theme$box.rectangle$lwd <- 0.5
bw_theme$box.rectangle$fill <- "grey90"
bw_theme$box.umbrella$lty <- 1
bw_theme$box.umbrella$col <- "black"
bw_theme$plot.symbol$col <- "grey40"
bw_theme$plot.symbol$pch <- "*"
bw_theme$plot.symbol$cex <- 1.5
bw_theme$strip.background$col <- "grey80"

################
## Bulk plots ##
################

# Test plot
plot(d15n ~ linear.ad, df.jack, type = "l", col = alpha("black", 0.3), lwd = 1.5, ylim = c(6.5, 12), ylab = n)
lines(d15n.5pt ~ linear.ad, df.jack, col = "blue", lwd = 2)

# Test plot
par(pty = "s")
plot(d15n ~ d13c, df.stet, col = df.stet$Col, cex = 1.5, pch = 20)

##############################################
## Linear regressions between d13C and d15N ## 
## for all colonies/specimens               ##
##############################################

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

#################
# Stetson Banks #
#################

# Stetson Banks bulk d15N
# par(pty = "s", mfrow=c(2,1))
par(pty = "s")
plot(d15n ~ linear.ad2, df.stet,
     xlab = x,
     ylab = n,
     type = "l",
     cex = 0.5,
     # xlim = c(1200,1300),
     col = alpha("black", 0.3))
lines(d15n.3pt ~ linear.ad2, df.stet,
      col = "#33a02c", lwd = 3)
# lines(d15n ~ linear.ad, df.jack4684,
#       col = alpha("black", 0.3))
lines(d15n.3pt ~ linear.ad, df.jack4684,
      col = "#1f78b4", lwd = 3)


plot(d13c ~ linear.ad2, df.stet,
     xlab = x,
     ylab = c,
     type = "l",
     xlim=c(1200,1300),
     col = alpha("black", 0.4))
lines(d13c.5pt ~ linear.ad2, df.stet,
      col = "red", lwd = 2)


# with ggplot, focus on LIA and forward
s1 <- ggplot(df.stet, aes(x=linear.ad, y=d15n))
s1 + geom_line(color = alpha("#009E73", 0.01)) +
  coord_fixed(ratio = 160) +
  geom_line(data=df.stet, aes(x=linear.ad, y=d15n), color = "#009E73", size = 0.75) +
  geom_line(data=df.jack4684, aes(x=linear.ad, y=d15n), color = alpha("#D55E00", 0.01)) +
  geom_line(data=df.jack4684, aes(x=linear.ad, y=d15n), color = "#D55E00", size = 0.75) +
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

# Stetson Banks bulk d13C

s2 <- ggplot(df.stet, aes(linear.ad, d13c))
s2 + geom_line() +
  geom_point(shape = 21, fill = "gray", color = "black", size = 3) +
  # geom_ribbon(aes(ymin = d13c + 0.3, ymax = d13c - 0.3), alpha = 0.2, color = "black", size = 0.25) +
  # geom_point(shape = 21, size = 3, color = "black", fill = "red") +
  # theme_classic() +
  # xlim(1400, 2010) +
  # scale_x_continuous(breaks = seq(500,2010,by=200), limits = c(500, 2010)) +
  # scale_x_continuous(breaks = seq(1300,2010,by=100), limits = c(1300, 2010)) +
  # scale_x_continuous(breaks = seq(1000,2010,by=100), limits = c(1000, 2010)) +
  # scale_x_continuous(breaks = seq(1000,1500,by=100), limits = c(1000, 1500)) +
  # scale_x_continuous(breaks = seq(1650, 1800, by = 25), limits = c(1650, 1800)) +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_text(size=10, color = "black"),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75)) +
  xlab(x) +
  ylab(c)

par(pty = "s")
plot(d13c.corrected ~ linear.ad2, df.stet,
     xlab = "Calendar Years (C.E.)",
     ylab = expression("Suess-corrected "*delta^{13}*"C (\u2030)"),
     type = "l",
     col = alpha("black", 0.4))
lines(forecast::ma(d13c.corrected, order = 5, centre = TRUE) ~ linear.ad2, df.stet,
      col = "red", lwd = 2)

#' -------------------------------------------------------------
#' Figure comparing Stetson d13C and d15N together for ancient
#' biogeochemical regimes
#' -------------------------------------------------------------

# df.stet %>% 
#   select(linear.ad4, d15n, d13c) -> t.stet
# t.stet %>%
#   melt(., "linear.ad4") -> t.stet
# 
# t.stet$variable <- factor(t.stet$variable, labels = c(expression(paste(delta^{15}*'N (\u2030)')),
#                                                        expression(paste(delta^{13}*'C (\u2030)'))))
# 
# ggplot(data = t.stet, aes(x = linear.ad4, y = value, shape = 21)) +
#   # annotate("rect",xmin=1600,xmax=1850,
#   #          ymin=-Inf,ymax=Inf,fill="#9ecae1",color=NA,size=0.25,alpha=0.25) +
#   # annotate("rect",xmin=950,xmax=1250,
#   #          ymin=-Inf,ymax=Inf,fill="#e9a3c9",color=NA,size=0.25,alpha=0.25) +
#   geom_point(aes(color = "black", fill = factor(variable))) +
#   facet_wrap(~ variable, scales = "free_y",
#              strip.position = "left",
#              nrow = 4,
#              # labeller = as_labeller(proxynames)) +
#              labeller = label_parsed) +
#   xlab(x) +
#   ylab(NULL) +
#   xlim(500, 2010) +
#   theme_classic() +
#   theme(strip.background = element_blank(), strip.placement = "outside", legend.position = "none")
#   # geom_vline(xintercept=600, color="black", linetype="dashed")
# # ggsave("stetson_sst_cores.pdf", width=8, height=9)

#' -------------------------------------------------------------
#' Trying above again but with some different code: https://gist.github.com/tomhopper/faa24797bb44addeba79
#' 
#' -------------------------------------------------------------

plot1 <- df.stet %>%
  select(linear.ad4, d15n, d15n.5pt) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad4, y = d15n), size = 0.5, alpha = 0.75) +
  geom_point(shape = 16, color = "#f46d43", alpha = 0.4, size = 1.85) +
  geom_line(aes(y=d15n.5pt), size = 1.0, color = 'black') +
  ylab(n) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_blank(),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
        # axis.line = element_line(colour = "black"),
        # panel.border = element_rect(colour = "black", fill=NA, size=0.75))
plot1

plot2 <- df.stet %>%
  select(linear.ad4, d13c, d13c.5pt) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad4, y = d13c), size = 0.5, alpha = 0.75) +
  geom_point(shape = 16, color = "#4575b4", alpha = 0.4, size = 1.85) +
  geom_line(aes(y=d13c.5pt), size = 1.0, color = 'black') +
  ylab(c) +
  theme_classic() +
  xlab(x) +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_text(size=10, color = "black"),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
        # axis.line = element_line(colour = "black"),
        # panel.border = element_rect(colour = "black", fill=NA, size=0.75))
plot2

grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))


################
# Jacksonville #
################
# Disk 3
par(pty = "s", mfrow = c(1,2))
par(pty = "s")
plot(d15n ~ linear.ad, df.sav,
     xlab = x,
     ylab = n,
     type = "l",
     xlim = c(0, 1400),
     col = alpha("black", 0.0))
lines(rollmean(d15n, 3, na.pad = TRUE) ~ linear.ad, df.jack,
      col = "red", lwd = 2, type = "l")
lines(rollmean(d15n, 3, na.pad = TRUE) ~ linear.ad, df.sav,
      col = "blue", lwd = 2, type = "l")
lines(rollmean(d15n, 15, na.pad = TRUE) ~ linear.ad, df.stet, col = "purple", lwd = 2)

plot(d13c ~ linear.ad, df.jack,
     xlab = x,
     ylab = c,
     type = "l",
     col = alpha("black", 0.4))
lines(rollmean(d13c, 3, na.pad = TRUE) ~ linear.ad, df.jack, # rollmean(d13c, 3, na.pad = TRUE))
      col = "red", lwd = 2)
lines(rollmean(d13c, 3, na.pad = TRUE) ~ linear.ad, df.stet,
      col = "purple", lwd = 2)
lines(rollmean(d13c, 3, na.pad = TRUE) ~ linear.ad, df.sav,
      col = "blue", lwd = 2)

j <- ggplot(df.jack, aes(linear.ad, d13c))
j + geom_line() +
  # geom_point(shape = 21, size = 3, color = "black", fill = "red") +
  # theme_classic() +
  # xlim(1400, 2010) +
  # scale_x_continuous(breaks = seq(800,2010,by=100), limits = c(800, 2010)) +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_text(size=10, color = "black"),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75)) +
  xlab(x) +
  ylab(c)

# par(pty = "s")
# plot(d15n ~ linear.ad2, df.jack, # Using overall growth rate
#      xlab = "Calendar Years (C.E.)",
#      ylab = n,
#      type = "l",
#      col = alpha("black", 0.4))
# lines(d15n.5pt ~ linear.ad2, df.jack,
#       col = "red", lwd = 2)
# lines(d15n.3pt ~ linear.ad, df.sav,
#       col = "blue", lwd = 2)

# par(pty = "s")
# plot(d15n ~ rev(test), df.jack, # Using overall growth rate but starting from oldest peel instead of youngest
#      xlab = "Calendar Years (C.E.)",
#      ylab = n,
#      type = "l", xlim = c(-1000,1500),
#      col = alpha("black", 0.4))
# lines(d15n.5pt ~ rev(test), df.jack,
#       col = "red", lwd = 2)
# lines(d15n.3pt ~ linear.ad, df.sav,
#       col = "blue", lwd = 2)
# lines(d15n.3pt ~ linear.ad, df.stet,
#       col = "purple", lwd = )2

par(pty = "s")
plot(d13c ~ linear.ad, df.jack,
     xlab = "Calendar Years (C.E.)",
     ylab = c,
     type = "l",
     xlim = c(-1000,2000),
     col = alpha("black", 0.4))
abline(h=mean(df.jack$d13c))
lines(d13c.3pt ~ linear.ad, df.jack,
      col = "red", lwd = 2)
lines(d13c ~ linear.ad, df.sav, col = "blue")
points(d13c ~ linear.ad, df.sav, col = "blue")

par(mfrow = c(1,2), pty = "s")
plot(d15n ~ linear.ad, df.jack,
     xlab = x,
     ylab = n,
     type = "l",
     col = alpha("black", 0.4),
     lwd = 2)
lines(d15n.3pt ~ linear.ad, df.jack2, col = "blue", lwd = 2)
points(d15n.3pt ~ binded2$yrs, df.jack2, col = "blue", type = "o")

plot(d13c ~ linear.ad, df.jack,
     xlab = x,
     ylab = c,
     type = "l",
     col = alpha("black", 0.4),
     lwd = 2)

lines(d13c.3pt ~ linear.ad, df.jack2, col = "red", lwd = 2)
points(d13c.3pt ~ linear.ad, df.jack2, col = "red")

# Disk 1
par(mfrow = c(1,2), pty = "s")
plot(d15n ~ linear.ad, df.jack2,
     xlab = "Calendar Years (C.E.)",
     ylab = c,
     type = "o",
     col = alpha("black", 0.4))
# lines(d15n.3pt ~ linear.ad, df.jack2,
#       col = "red", lwd = 2)
# points(d15n ~ linear.ad, df.jack, pch = 23, col = alpha("black", 0.4))
lines(d15n.3pt ~ linear.ad, df.jack, col = alpha("blue", 0.4))

par(pty = "s")
plot(d13c ~ linear.ad, df.jack2,
     xlab = "Calendar Years (C.E.)",
     ylab = c,
     type = "o",
     col = alpha("black", 0.4))
lines(d13c.3pt ~ linear.ad, df.jack2,
      col = "red", lwd = 2)
points()
lines()


plotj1 <- df.jack %>%
  select(linear.ad, d15n, d15n.3pt) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d15n), size = 0.5, alpha = 0.75) +
  geom_point(shape = 21, color = "black", fill = "blue", alpha = 0.3) +
  geom_line(aes(y=d15n.3pt), size = 1.0) +
  ylab(n) +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(color = "black"))
plotj1

plotj2 <- df.jack %>%
  select(linear.ad, d13c, d13c.3pt) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d13c), size = 0.5, alpha = 0.75) +
  geom_point(shape = 22, color = "black", fill = "red", alpha = 0.3) +
  geom_line(aes(y=d13c.3pt), size = 1.0) +
  ylab(c) +
  theme_classic() +
  xlab(x) +
  theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"))
plotj2

grid.newpage()
grid.draw(rbind(ggplotGrob(plotj1), ggplotGrob(plotj2), size = "last"))


##################
# Savannah Banks #
##################
par(pty = "s")
plot(d15n ~ linear.ad, df.sav,
     xlab = "Calendar Years (C.E.)",
     ylab = n,
     type = "l",
     col = alpha("black", 0.4),
     xlim = c(0,1500))
lines(rollmean(d15n, 3, na.pad = TRUE) ~ linear.ad, df.sav,
      col = "red", lwd = 2)
lines(rollmean(d15n, 3, na.pad = TRUE) ~ linear.ad, df.jack,
      col = "blue", lwd = 2)
lines(rollmean(d15n, 3, na.pad = TRUE) ~ linear.ad, df.stet,
      col = "purple", lwd = 2)

par(pty = "s")
plot(d13c ~ linear.ad, df.sav,
     xlab = "Calendar Years (C.E.)",
     ylab = c,
     type = "l",
     col = alpha("black", 0.0))
lines(rollmean(d13c, 3, na.pad = TRUE) ~ linear.ad, df.sav,
      col = "red", lwd = 2)
lines(rollmean(d13c, 3, na.pad = TRUE) ~ linear.ad, df.jack,
      col = "blue", lwd = 2)
lines(rollmean(d13c, 5, na.pad = TRUE) ~ linear.ad, df.stet,
      col = "purple", lwd = 2)



plots1 <- df.sav %>%
  select(linear.ad, d15n) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d15n), size = 0.5, alpha = 0.75) +
  geom_point(shape = 21, color = "black", fill = "blue", alpha = 0.4) +
  geom_line(aes(y=d15n), size = 0.25) +
  ylab(n) +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(color = "black"))
plots1

plots2 <- df.sav %>%
  select(linear.ad, d13c) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d13c), size = 0.5, alpha = 0.75) +
  geom_point(shape = 22, color = "black", fill = "red", alpha = 0.4) +
  geom_line(aes(y=d13c), size = 0.25) +
  ylab(c) +
  theme_classic() +
  xlab(x) +
  theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"))
plots2

grid.newpage()
grid.draw(rbind(ggplotGrob(plots1), ggplotGrob(plots2), size = "last"))



#####################
# Jacksonville-4684 #
#####################
par(pty = "s")
plot(d15n ~ linear.ad, df.jack4684,
     xlab = "Calendar Years (C.E.)",
     ylab = n,
     # xlim=c(1300, 2005),
     type = "l",
     col = alpha("red", 0.8))
lines(rollmean(d15n, 3, na.pad = TRUE) ~ linear.ad, df.jack4684,
      col = "black", lwd = 1.25)
points(rollmean(d15n, 3, na.pad = TRUE) ~ linear.ad, df.stet,
      col = "blue", lwd = 1.5)

# with xyplot (lattice)
xyplot(d15n ~ linear.ad, 
       data = df.jack4684,
       ylab = n,
       xlab = x,
       type = "o",
       pch = 21, cex = 1.5, col = "black", fill = "maroon")
xyplot(d13c ~ linear.ad, 
       data = df.jack4684,
       ylab = n,
       xlab = x,
       ylim = c(-17, -14),
       type = "o",
       pch = 21, cex = 1.5, col = "black", fill = "maroon")

par(pty = "s")
plot(d15n ~ linear.ad, df.jack4684,
     xlab = x,
     ylab = n,
     type = "l",
     xlim = c(-250,2000),
     ylim = c(7, 10),
     col = alpha("black", 0.01))
lines(d15n.3pt ~ linear.ad, df.jack4684, col = '#1a9850', lwd = 2)
lines(d15n.3pt ~ linear.ad, df.jack, col = '#d73027', lwd = 2)

par(pty = "s")
plot(d13c ~ linear.ad, df.jack4684,
     xlab = "Calendar Years (C.E.)",
     ylab = c,
     type = "o",
     col = alpha("black", 0.4))
lines(d13c.3pt ~ linear.ad, df.sav,
      col = "red", lwd = 2)


ggplot(data=df.jack4684, aes(linear.ad, d15n)) +
  theme_bw() +
  xlab(x) +
  ylab(n) +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_text(size=10, color = "black")) +
  geom_line(color = alpha("black", 0.75)) +
  geom_point(shape=21, fill = "gray", size = 2) +
  # geom_text(data=df.jack4684, aes(label=df.jack4684$sample.no.), hjust=0, vjust=-1) +
  # xlim(1200, 2000) +
  geom_point(data=df.jack4684[df.jack4684$sample.no. == "36" | df.jack4684$sample.no. == "39"
                          | df.jack4684$sample.no. == "41"
                          | df.jack4684$sample.no. == "85" | df.jack4684$sample.no. == "155",],
             pch=21, fill="red", size=2.75)

ndata %>%
  filter(Sample.ID2 == "Jacksonville-4684") -> t.ndata
obj2 <- xyplot(d15n ~ linear.ad, data = df.jack4684, type = "l", lwd = 2, xlim = c(1500, 2005),
               ylab = expression(delta^{15}*"N (\u2030)"), xlab = x, col = '#f03b20')
obj1 <- xyplot(Phe ~ Year.AD, data = t.ndata, pch = 1, cex = 2, col = alpha("#d95f0e", 1), lwd = 15, xlim = c(1500, 2005),
               xlab = x, ylab = expression(delta^{15}*"N"[" Phe"]*" (\u2030)"))
doubleYScale(obj2, obj1, add.ylab2 = TRUE, style1= NULL, style2 = NULL)

#' Same idea as above but more scientifically sound since it does not
#' do a double Y axis, which is controversial in some circles (I tend to agree, mostly)
#' 


pl1 <- df.jack4684 %>%
  select(linear.ad, d15n) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d15n), size = 0.5, alpha = 0.75) +
  geom_point(shape = 21, color = "black", fill = "blue", alpha = 0.3) +
  geom_line(aes(y=rollmean(d15n, 3, na.pad = TRUE)), size = 0.5) +
  ylab(n) +
  xlim(1600,2005) +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(color = "black"))
pl1

pl2 <- t.ndata %>%
  select(Year.CE, Phe) %>%
  na.omit() %>%
  ggplot(aes(x = Year.CE, y = Phe)) +
  geom_point(shape = 23, color = "black", fill = "red", alpha = 0.75, size = 4) +
  # geom_line(aes(y=Phe), size = 1.0) +
  ylab(phe) +
  theme_classic() +
  xlab(x) +
  xlim(1600,2005) +
  ylim(1,8) +
  theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"))
pl2

grid.newpage()
grid.arrange(rbind(ggplotGrob(pl1), ggplotGrob(pl2), size = "last"))



##################### 
# Jacksonville-4686 #
#####################

par(pty = "s")
plot(d15n ~ distance, df.jack4686,
     xlab = x,
     ylab = n,
     type = "o",
     col = alpha("black", 0.8))

par(pty = "s")
plot(d13c ~ distance, df.jack4686,
     xlab = "Calendar Years (C.E.)",
     ylab = c,
     type = "o",
     col = alpha("black", 0.4))
lines(d13c.3pt ~ distance, df.jack4686,
      col = "red", lwd = 2)

zz <- ggplot(df.jack4686, aes(x=distance, y=d15n))
zz + geom_line() +
  geom_point() +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_text(size=10, color = "black"),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75)) +
  xlab("Distance from edge (mm)") +
  ylab(n)


##################################
## Load Jacksonville-4907 BAM1  ## 
## bamboo coral data            ##
##################################

bam <- read.csv('C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff_jack4907 bam1 by treatment.csv')
bam <- read.csv('~/Google Drive/projects/rproj/seus/data/schiff_jack4907 bam1 by treatment.csv')








############################
## Lattice plots showing  ## 
## bulk d13C vs bulk d15N ##
############################
df.bulk %>%
  filter(coral.id == "jack-4907-bc1-d1"  | coral.id == "sav-4902-bc1-unk" | coral.id == "stet-4904-bc1-d2") %>%
  select(coral.id, d13c, d15n, linear.ad) %>%
  na.omit(.) -> t.df

df.stet %>%
  filter(linear.ad < 1860) -> df.stet
t.lm <- lm(d15n ~ d13c, df.stet)
summary(t.lm)

plot(d15n ~ d13c, df.stet, col = df.stet$Col, cex = 1.5, pch = 20)
qplot(d13c, d15n, data=df.stet, colour=df.stet$linear.ad) + 
  scale_colour_gradient(low="red", high="blue") +
  labs(colour = "Calendar Years (C.E.)",
       x = c,
       y = n) +
  # theme(aspect.ratio = 1) +
  theme_bw()

mycol <- brewer.pal(6, "Blues")
xyplot(d15n ~ d13c | coral.id,
       group = linear.ad,
       data = t.df,
       ylab = n,
       xlab = c,
       strip = strip.custom(factor.levels = c("Jacksonville Lithoherms", "Savannah Banks","Stetson Banks"),
                            par.strip.text=list(col="white", font=1.75)),
       par.settings = list(superpose.symbol = list(pch = 21,
                                                   cex = 1,
                                                   col = "black",
                                                   fill = "lightblue"),
                           strip.background = list(col = mycol[6])))



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

ndata %>%
  filter(Region == "SEUS") %>%
  melt(., "Sample.ID") %>%
  droplevels(.) %>%
  filter(variable %in% tr.srcaa) %>% # Success! 03-04-2019
  xyplot(value ~ variable,
         .,
         panel = function(x,  y, ...){
           panel.xyplot(x, y, ...)
           # panel.abline(h = 9, lty = 1)
           # panel.abline(h = 9.25, lty = 2)
           # panel.abline(h = 8.75, lty = 2)
           panel.abline(v = 7.5, lty = 2)
         },
         cex=1.5,
         ylim=c(-20, 40), # ylim and xlim are "first class" parameters and don't need to be in scales=list()
         group = Sample.ID,
         xlab = NULL,
         ylab = expression({delta}^15*"N (\u2030)"),
         auto.key=list(columns=2, cex = 0.75),
         par.settings = my.settings2)


ggplot(value ~ variable, .)


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



# same plot in baseR -- wokring on this
plot.data <- ndata.melt[c(136:720),]
plot.data$value <- as.numeric(plot.data$value)
plot(value ~ variable1, plot.data)

######################
## Trophic Dynamics ## 
######################
xyplot(TP ~ Sum.V | Region, # SumV vs. Trophic Position
       ndata[1:29,],
       group = Sample.ID2,
       pch = 21,
       cex = 1.5,
       layout = c(1,2),
       # auto.key=list(columns=c(1,2), cex = 0.75),
       xlab = expression(paste(Sigma,"V")),
       ylab = "Trophic Position (Glu - Phe)")

# For only the SEUS and GOM regions
ndata %>%
  filter(Type == "Black Coral - GOM" | Type == "Black Coral - SEUS") -> bc.ndata
  
p <- ggplot(bc.ndata, aes(x=Sum.V, y=TP, shape = Region, fill = Region), color = "black")
p + geom_point(size = 3) +
  facet_wrap(~Region, ncol = 1) +
  scale_shape_manual(values = c(23, 23)) +
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

ndata %>%
  filter(Region == "SEUS") %>%
  droplevels(.) -> t.ndata

ggplot(aes(x = Year.CE, y = Sum.V), data = t.ndata) +
  geom_point(aes(fill = factor(Sample.ID2), shape = factor(Sample.ID2)), size = 3.5) +
  # scale_color_grey() +
  scale_fill_manual(values=c("#636363","#bdbdbd","#f0f0f0")) +
  scale_shape_manual(values=c(22,23,24))+
  theme_classic() +
  theme(legend.title=element_blank()) +
  xlab(x) +
  ylab(expression(paste(Sigma,"V"))) +
  geom_vline(xintercept=950, linetype = 'dotted') +
  geom_vline(xintercept=1250, linetype = 'dotted') +
  geom_text(x=1100, y=3.5, label = "MCA")
  

########################
## TP through time    ##
########################
ndata %>%
  filter(Region == "SEUS") %>%
  droplevels(.) %>%
  xyplot(TP ~ Year.CE,
         data = .,
         group = Sample.ID2,
         cex = 1.5,
         # xlim = c(750, 1500),
         par.settings = my.settings,
         auto.key = list(columns=c(2), cex = 0.95),
         xlab = x,
         ylab = "Trophic Position (Glu - Phe)")

ggplot(aes(x = Year.CE, y = TP), data = t.ndata) +
  geom_point(aes(fill = factor(Sample.ID2), shape = factor(Sample.ID2)), size = 3.5) +
  scale_fill_manual(values=c("#636363","#bdbdbd","#f0f0f0")) +
  scale_shape_manual(values=c(22,23,24))+
  theme_classic() +
  theme(legend.title=element_blank()) +
  # xlab(x) +
  xlab("Year CE") +
  ylab("Trophic Position (Glu - Phe)") +
  ylim(1,3) +
  # geom_vline(xintercept=950, linetype = 'dotted') +
  # geom_vline(xintercept=1250, linetype = 'dotted') +
  geom_text(x=1100, y=3.5, label = "MCA")


######################
## Phe through time ##
######################
ndata %>% 
  # filter(Phe > 1) %>%
  # filter(Phe < 13) %>%
  # filter(Sample.ID2 == "Jacksonville-4907" | Sample.ID2 == "Savannah Banks-4902") %>%
  filter(Region == "SEUS") %>%
  droplevels(.) %>%
  xyplot(Phe ~ Year.CE,
       data = .,
       group = Sample.ID2,
       cex = 1.5,
       par.settings = my.settings,
       auto.key = list(columns=c(2), cex = 0.75),
       xlab = x,
       ylab = expression({delta}^15*"N"[" Phe"]*" (\u2030)"))
ndata %>%
  filter(Region == "SEUS") -> t.ndata
t.ndata$bp <- 1950 - t.ndata$Year.CE

ggplot(aes(x = Year.CE, y = SrcAA), data = t.ndata) +
  geom_point(aes(fill = factor(Sample.ID2), shape = factor(Sample.ID2)), size = 3.5) +
  scale_fill_manual(values=c("#636363","#bdbdbd","#f0f0f0")) +
  scale_shape_manual(values=c(22,23,24))+
  theme_classic() +
  theme(legend.title=element_blank()) +
  # xlab(x) +
  xlab("Year BP") +
  ylab(phe) +
  geom_vline(xintercept=950, linetype = 'dotted') +
  geom_vline(xintercept=1250, linetype = 'dotted') +
  geom_text(x=1100, y=3.5, label = "MCA")

ndata %>% # Showing Phe through time, but only for select specimens (modify code as needed)
  filter(Phe > 1) %>%
  filter(Phe < 13) %>%
  filter(Sample.ID2 == "Jacksonville-4907" | Sample.ID2 == "Savannah Banks-4902" | Sample.ID2 == "Jacksonville-4684") -> t.ndata

par(mfrow=c(1,2))
xyplot(Phe ~ Year.AD,
     data = t.ndata,
     groups = Sample.ID2,
     # xlim =c(1500,1900),
     xlab = x,
     xlim = c(-250,1500),
     ylab = phe,
     cex = 3,
     ylim = c(4, 14))


plot(SrcAA ~ Year.CE, 
     t.ndata,
     col = alpha("black", 1.0),
     xlab = x,
     xlim = c(1500, 2005),
     ylim = c(7, 10),
     ylab = n,
     lwd = 1.25)
lines(forecast::ma(d15n, order = 1, centre = TRUE) ~ linear.ad, df.stet, col = "red")
lines(forecast::ma(d15n, order = 1, centre = TRUE) ~ linear.ad, df.jack4684, col = "blue")
  # xyplot(Phe ~ Year.AD,
  #        data = .,
  #        group = Sample.ID2,
  #        cex = 1.5,
  #        par.settings = my.settings,
  #        auto.key = list(columns=c(2), cex = 0.75),
  #        xlab = x,
  #        ylab = expression({delta}^15*"N"[" Phe"]*" (\u2030)"))

#################################
## Average Src-AA through time ##
#################################
ndata %>%
  dplyr::select(Gly, Ser, Phe) -> t.data # Select the Src-AA columns
ndata$SrcAA <- rowMeans(t.data, na.rm = TRUE)


################################
## Glu vs. Phe                ##
################################
xyplot(Glu ~ Phe,
       data = ndata,
       group = Region,
       par.settings = my.settings,
       # pch = c(16, 17, 18, 19, 8),
       # pch = c(21,22,23,24,25),
       cex = 1.5,
       auto.key = list(columns = c(2),
                       text = c("POM (Central NA)", "POM (Eastern NA)",
                                   "Black Coral - GOM", "Black Coral - SEUS", "POM (Western NA)")),
       ylim = c(0, 25),
       xlim = c(-10, 15),
       ylab = expression({delta}^15*"N"[" Glu"]*" (\u2030)"),
       xlab = expression({delta}^15*"N"[" Phe"]*" (\u2030)"),
       panel = function(...) {
         panel.abline(coef = c(3.4,1), lty = "dashed") # TP = 1
         panel.abline(coef = c(11,1), lty = "dashed") # TP = 2
         panel.abline(coef = c(18.6,1), lty = "dashed") # TP = 3
         panel.abline(coef = c(26.2,1), lty = "dashed") # TP = 4
         panel.abline(coef = c(33.8,1), lty = "dashed")
         panel.xyplot(...)
       })

d <- ggplot(ndata, aes(x = Phe, y = Glu))
d + geom_point() + 
  xlim(-10, 15) +
  ylim(0,25) +
  geom_abline(intercept = 3.4, slope = 1) +
  geom_abline(intercept = 11, slope = 1) +
  geom_abline(intercept = 18.6, slope = 1) +
  geom_abline(intercept = 26.2, slope = 1) +
  geom_abline(intercept = 33.8, slope = 1)

#' Remember, Trophic Position is calculated thus:
#' TP = (d15Nglu - d15Nphe - 3.4)/7.6 + 1


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

seus_carbon %>% 
  # select(-c(Bulk)) %>%
  melt(., "Group.ID") -> melt # Nifty trick using magrittr pipe operator, %>%

melt$variable1 <- factor(melt$variable, 
                               levels=c(eaa.neaa)) 

xyplot(value ~ variable1,
       melt,
       panel = function(x,  y, ...){
         panel.xyplot(x, y, ...)
         # panel.abline(h = -15.9, lty = 1) # Average bulk value
         # panel.abline(h = (-15.9+0.54), lty = 2) # + standard deviation
         # panel.abline(h = (-15.9-0.54), lty = 2) # - standard deviation
         # panel.abline(h = -11.4, lty = 1)
         # panel.abline(h = (-11.4+0.8), lty = 2)
         # panel.abline(h = (-11.4-0.8), lty = 2)
         panel.abline(v = 5.5, lty = 2)
       },
       cex=1.5,
       ylim=c(-30, 10), # ylim and xlim are "first class" parameters and don't need to be in scales=list()
       group = Group.ID,
       xlab = NULL,
       ylab = expression({delta}^13*"C (\u2030)"),
       auto.key=list(columns=3, cex = 0.55),
       par.settings = my.settings3)

# Same plot using ggplot
melt %>% 
  drop_na() -> melt
aa <- ggplot(melt, aes(variable1, value))
aa + geom_point(aes(color = Group.ID, shape = Group.ID), size = 3, stroke = 1.0) +
  theme_classic() +
  scale_shape_manual(NULL, values=c(0:10)) +
  scale_color_manual(NULL, values = c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')) +
  xlab(NULL) +
  # scale_x_discrete(limits=c('Gly', 'Ser', 'Asx', 'Pro', 'Ala', 'Thr', 'Ile', 'Val', 'Phe', 'Leu')) + # ordered like in McMahon et al (2018)
  ylab(expression({delta}^13*"C (\u2030)")) +
  theme(axis.text=element_text(size=12, color = "black"),
        legend.position = "top") +
  geom_vline(xintercept = 5.5, linetype="dashed", 
             color = "black", size=0.75)

#' 
#' We want to plot just the normalized means with standard errors for 
#' each group (group = Sample.ID2). I do this below
#' 

melt(set, "Group") %>% # This is a good plot to publish, 3-2-2019
  xyplot(value ~ variable,
         .,
         panel = function(x,  y, ...){
           panel.xyplot(x, y, ...)
           panel.abline(h = -15.9, lty = 1)
           panel.abline(h = (-15.9+0.54), lty = 2)
           panel.abline(h = (-15.9-0.54), lty = 2)
           },
         cex=1.5,
         ylim=c(-15, 30), # ylim and xlim are "first class" parameters and don't need to be in scales=list()
         group = Group,
         xlab = NULL,
         ylab = expression({delta}^13*"C (\u2030)"),
         auto.key=list(columns=2, cex = 0.95),
         par.settings = my.settings)

melt(set, "Group") %>%
  ggplot(., aes(y=value, x=variable, fill = Group), color = "black") +
  # scale_color_manual(values = c('#2c7fb8', '#7fcdbb')) +
  scale_fill_manual(values = c('#2c7fb8', '#7fcdbb')) +
  scale_shape_manual(values = c(21,23)) +
  labs(y = expression("Normalized "*{delta}^13*"C (\u2030)"), x = NULL) +
  geom_point(aes(shape = Group), size = 4.5) -> p1
melt(set, "Group") %>%
  ggplot(., aes(y=value, x=variable, color = Group)) +
  scale_color_manual(values = c('#2c7fb8', '#7fcdbb')) +
  # scale_fill_manual(values = c('#2c7fb8', '#7fcdbb')) +
  scale_shape_manual(values = c(21,23)) +
  labs(y = expression("Normalized "*{delta}^13*"C (\u2030)"), x = NULL) +
  geom_point(aes(shape = Group), size = 4.5, stroke = 2) -> p2
p2 + theme_bw()

#' 
#' Below are exploratory charts
#' Not final ones
#' 

xyplot(EAA ~ Year.CE, # Average EAA through time compared to bulk d13C
       data = seus_carbon,
       type = "o",
       pch = 21,
       cex = 1.5)

xyplot(Gly ~ Year.CE, # Average EAA through time compared to bulk d13C
       data = seus_carbon,
       type = "o",
       # ylim = c(-25, -15),
       pch = 21,
       cex = 1.5)

plot(Bulk ~ Year.CE,
     data = seus_carbon,
     type = "o", ylim = c(-18, -14),
     ylab = expression({delta}^13*"C (\u2030)"),
     xlab = "Calendar Year (C.E).")
points(Glx ~ Year.CE,
       data = seus_carbon,
       col = "blue")

obj2 <- xyplot(EAA ~ Year.CE, data = seus_carbon, pch = 21, cex = 1.5, type = "o", xlim = c(1300, 2005),
               ylab = expression("Reconstructed RPP "*delta^{13}*"C (\u2030)"), xlab = x)
obj1 <- xyplot(d13c ~ linear.ad, data = df.stet, type = "o", col = alpha("black", 0.4), lwd = 1.5, xlim = c(1300, 2005),
               xlab = x, ylab = expression(delta^{13}*"C (\u2030)"))
doubleYScale(obj1, obj2, add.ylab2 = TRUE)

obj2 <- xyplot(EAA ~ Year.CE, data = seus_carbon, pch = 21, cex = 1.5, type = "o", xlim = c(1300, 2005),
               ylab = expression("EAA "*delta^{13}*"C (\u2030)"), xlab = x)
obj1 <- xyplot(d13c.5pt ~ linear.ad, data = df.stet, type = "l", col = alpha("black", 0.4), lwd = 1.5, xlim = c(1300, 2005),
               xlab = x, ylab = expression(delta^{13}*"C (\u2030)"))
doubleYScale(obj1, obj2, add.ylab2 = TRUE)

cplot1 <- seus_carbon %>%
  dplyr::select(EAA, Year.CE) %>%
  na.omit() %>%
  ggplot(aes(x = Year.CE, y = EAA), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.65) +
  # geom_line(aes(y=rollmean(EAA, 3, na.pad = TRUE)), size = 0.65, color = '#225ea8') +
  
  ylab(n) +
  xlab(NULL) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(1300,2005)) +
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

cplot2 <- df.stet %>%
  dplyr::select(linear.ad, d13c) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d13c), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.65) +
  geom_line(aes(y=rollmean(d13c, 3, na.pad = TRUE)), size = 0.65, color = '#081d58') +
  
  ylab(c) +
  theme_classic() +
  xlab(x) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(1300,2005)) +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_text(size=10, color = "black"),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

grid.newpage()
grid.draw(rbind(ggplotGrob(cplot1), ggplotGrob(cplot2), size = "last"))

#' Fig


#' Note 3-1-2019
#' Make sure to make some charts with the C-CSIAA data and send it off to Nancy
#' 
#' 

#########################
## Plotting Mol % data ##
#########################
mol %>%
  melt(., "Sample") %>%
  bwplot(value ~ variable,
         data = .,
         horiz = FALSE,
         ylab = 'Amino acids (mol, %)',
         pch = "|",
         par.settings = bw_theme)
mol %>%
  melt(., "Sample") %>%
  barchart(value ~ variable,
         data = .,
         horiz = FALSE,
         ylab = 'Amino acids (mol, %)',
         col = "darkgrey")


