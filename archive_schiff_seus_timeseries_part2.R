#' -------------------------------------------------
#' Topic: Time series models using bulk data, updated
#' John Schiff
#' 02/08/2019
#' -------------------------------------------------
#' 
library(lattice)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(data.table)
library(Rssa)
#' Here, I am redoing some of the code for CN age models,
#' based on the code in schiff_seus_timeseries.R
#' I am importing the bulk data from a new file that does
#' not already have the estimated ages for each sample
#' attached. I will do that here as an extension of the
#' linear/clam age model scripts.

# Load in the data
# Create data tables for each coral of interest
pathmac1 <- '~/Google Drive/projects/rproj/seus/data/schiff cn age models 02-05-2019.csv'
pathwin1 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff cn age models 02-05-2019.csv'
bulk.import <- read.csv(pathmac1) # If on MacOS, don't alter this dataframe
bulk.import <- read.csv(pathwin1) # If on Windows, don't alter this dataframe
bulk <- bulk.import
bulk <- as.data.table(bulk)
# bulk.stet <- bulk %>% filter(coral.id == 'stet-4904-bc1-d2') # Using dplyr package
bulk.stet <- bulk[coral.id == 'stet-4904-bc1-d2'] # Using data.table package
bulk.jackusgs <- bulk[coral.id == 'jack-4907-bc1-d1']
bulk.jack <- bulk[coral.id == 'jack-4907-bc1-d3']
bulk.jack4684 <- bulk[coral.id == 'jack-4684-bc-unk']
bulk.jack4686 <- bulk[coral.id == 'jack-4686-bc-d1-t1']
bulk.sav <- bulk[coral.id == 'sav-4902-bc1-unk']

###########################################################################
# Now we need to add ages to these coral data tables.                     #
# This is where the distance column will be useful.                       #
# Essentially, we are creating a new column with estimated ages           #
# Let's do this with Stetson first by using depths from seus_clam_ages.R  #
###########################################################################

bulk.stet$linear.bp.min95 <- stet.depths$min95.
bulk.stet$linear.bp.max95 <- stet.depths$max95.
bulk.stet$linear.bp.fit <- stet.depths$best
bulk.stet$year.ad <- 1950 - bulk.stet$linear.bp.fit
bulk.stet$min.ad <- 1950 - bulk.stet$linear.bp.min95
bulk.stet$max.ad <- 1950 - bulk.stet$linear.bp.max95

ts.stet <- arrange(desc(bulk.stet$d15n.vs.air))
ts.stet <- ts(bulk.stet$d15n.vs.air, end = c(2005, 6))

bulk.jack4684$linear.bp.min95 <- jack4684.depths$min95.
bulk.jack4684$linear.bp.max95 <- jack4684.depths$max95.
bulk.jack4684$linear.bp.fit <- jack4684.depths$best
bulk.jack4684$year.ad <- 1950 - bulk.jack4684$linear.bp.fit
bulk.jack4684$min.ad <- 1950 - bulk.jack4684$linear.bp.min95
bulk.jack4684$max.ad <- 1950 - bulk.jack4684$linear.bp.max95

# bulk.jack$linear.bp.min95   These are all based on models in clam
# bulk.jack$linear.bp.max95
# bulk.jack$linear.bp.fit
# bulk.jack$spline.bp.min95
# bulk.jack$spline.bp.max95
# bulk.jack$spline.bp.fit
# bulk.jack$year.ad
# bulk.jack$min.ad
# bulk.jack$max.ad
bulk.jack$yearad.linear <- bind$yrs # This is the combined chronology from two differnt growth rates

############################################################################
## Combine two growth rate vectors into one to add to bulk.jack dataframe ##
############################################################################
distance <- bulk.jack$distance..mm.*1000
# dis1 <- bulk.jack %>% filter(distance > 2500)
# dis2 <- bulk.jack %>% filter(distance < 2500)
yrs1 <- (1950-588) - (distance/g2) # 588 is median year AD for the edge of the coral, g2 from linear_age_models.R script file
yrs2 <- (1950-2062) - (distance/g3) # 2062 is median year for the coral at distance 2558 microns from edge
y <- cbind(distance,yrs1,yrs2)
y <- as.data.frame(y)

first <- y[1:85,]
first <- first[,-c(3)]
second <- y[86:nrow(y),]
second <- second[,-c(2)]
names(first)[2] <- paste("yrs")
names(second)[2] <- paste("yrs")
bind <- rbind(first,second)
length(bind$distance)
############################################
## Creating linear growth rate chronology ##
## for Savannah Banks, BC1                ##
############################################

distance <- bulk.sav$distance..mm.*1000
yrs1 <- (1950-)




####################################################################
## After adding age estimation to a data sheet ...                ##
## ... save the data sheet as a csv to avoid having to redo this. ##
## Then import that data sheet.                                   ##
####################################################################

#' Saving Stetson Chronology data
write.csv(stet.table, file = "C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/stetson bulk read_only.csv")



################################################################
#' First step is to create moving averages of the isotope data 
#' Various packages do this (zoo, TTR, smooth, forecast)
#' For now I'll use SMA (simple moving average) in the
#' TTR package since it seems to have broad support
#' from statistical community.
################################################################
stet.table <- as.data.table(bulk.stet)
stet.table <- stet.table[, c("X.n", "X.c", "c.n", "machine") := NULL] # These columns are currently empty
stet.table <- na.omit(stet.table) # Remove NAs from d15N and d13C columns
nitrogen <- stet.table$d15n.vs.air
carbon <- stet.table$d13c.vs.vpdb
sma.nitrogen <- forecast::ma(nitrogen, n = 3) # Doing this the longer way
sma.carbon <- forecast::ma(carbon, n = 3) # ... for some reason TTR::SMA was being finicky
stet.table$d15n.ra <- sma.nitrogen
stet.table$d13c.ra <- sma.carbon

stet.ts <- ggplot(stet.table, aes(x=year.ad,y=d15n.ra)) +
  # theme_bw() +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(x="Calendar Year (C.E.)",
       y=expression(delta^{15}*"N (\u2030)")) +
  ylim(6,10) +
  scale_x_continuous(limits=c(350, 2005),
                     breaks=seq(-750, 2005, by=250)) +
  annotate("rect",xmin=1600,xmax=1850,
           ymin=-Inf,ymax=Inf,fill="#9ecae1",color=NA,size=0.25,alpha=0.25) +
  annotate("rect",xmin=950,xmax=1250,
           ymin=-Inf,ymax=Inf,fill="#e9a3c9",color=NA,size=0.25,alpha=0.25) +
  # annotate("rect",xmin=-900,xmax=-300,
  #          ymin=-Inf,ymax=Inf,fill="darkseagreen1",color="black",size=0.25,alpha=0.25) +
  # annotate("rect",xmin=-805,xmax=-795,
           # ymin=-Inf,ymax=Inf,fill="darkorange2",color="black",size=0.25,alpha=0.25) +
  # annotate("text", x=1075, y=7, label="MCA") +
  # annotate("text", x=1700, y=7, label="LIA") +
  # annotate("text", x=525, y=9.5, label="Bond 1") +
  scale_color_manual(values=c('#225ea8'),
                     name="Black Coral",
                     breaks=c('jack-4907-bc1-d3', 'sav-4902-bc1-unk', 'stet-4904-bc1-d2'),
                     labels=c('Jacksonville', 'Savannah', 'Stetson Banks')) +
  geom_line(aes(color=coral.id),size=0.75) +
  geom_vline(xintercept=600, color="black", linetype="dashed")
  # geom_vline(xintercept=1850, color="black")

stet.ts
# ggsave("stetson_d15n_5pt.pdf", width=7, height=4)
ggsave("stetson_d15n_5pt.png", width=6, height=4)

#' ----------------------------------------
#' Figure. Savannah BC1 through time
#' 
#' 
#' ----------------------------------------

bulk.sav$linear.bp.min95 <- sav.depths$min95.
bulk.sav$linear.bp.max95 <- sav.depths$max95.
bulk.sav$linear.bp.fit <- sav.depths$best
bulk.sav$year.ad <- 1950 - bulk.sav$linear.bp.fit

sav.ts <- ggplot(bulk.sav, aes(x=year.ad,y=d13c.ra)) +
  geom_line()
sav.ts

#' ----------------------------------------
#' Figure. Jacksonville-4686 through time
#' ----------------------------------------
#' 
jack4684.table <- as.data.table(bulk.jack4684)
jack4684.table <- jack4684.table[, c("X.n", "X.c", "c.n", "machine") := NULL] # These columns are currently empty
jack4684.table <- na.omit(jack4684.table) # Remove NAs from d15N and d13C columns
nitrogen <- jack4684.table$d15n.vs.air
carbon <- jack4684.table$d13c.vs.vpdb
sma.nitrogen <- TTR::SMA(nitrogen, n = 1) # Doing this the longer way
sma.carbon <- TTR::SMA(carbon, n = 1) # ... for some reason TTR::SMA was being finicky
jack4684.table$d15n.ra <- sma.nitrogen
jack4684.table$d13c.ra <- sma.carbon

##################################################
##                                              ##
## Analysis of Isotope Records after converting ##
## to R time series objects                     ##
##                                              ##
##################################################

#' NOTE (02-25/2019)
#' I think this is more trouble than it is worth.
#' Instead I am going to do an SSA on an object and then bind
#' the reconstructed vector to the overall dataframe
#' e.g., bind recon$F1 to bulk.stet
#' 

##################################################
##                                              ##
## Creating dummy data frames from bulk records ##
## with years and then attaching reconstructed  ##
## SSA objects to those.                        ##
##                                              ##
##################################################

# Stetson Banks coral

View(bulk.stet)
dt <- bulk.stet %>% select(d15n.vs.air, d13c.vs.vpdb, linear.ad)
dt <- na.omit(dt)

t.ssa <- ssa(dt$d15n.vs.air, L = 20)
plot(t.ssa)
plot(t.ssa, "series", groups = 1:20)
plot(t.ssa, "vectors", groups = 1:20)
plot(t.ssa, "paired", groups = 1:20)
plot(t.ssa, "wcor", groups = 1:20)
t.recon <- reconstruct(t.ssa, groups = list(1, 2, 3, 4, c(5:10)))
plot(t.recon, plot.method = "xyplot", type = "cumsum")

plot(d15n.vs.air ~ linear.ad, dt, col = alpha("black", 0.4), type = "l", ylim = c(0, 10))
lines((t.recon$F3 +t.recon$F2 + t.recon$F1) ~ dt$linear.ad, col = "blue", lwd = 2)
lines((t.recon$F3 +t.recon$F2 + t.recon$F1)+0.2 ~ dt$linear.ad, col = "black", lwd = 1.5, lty="dashed")
lines((t.recon$F3 +t.recon$F2 + t.recon$F1)-0.2 ~ dt$linear.ad, col = "black", lwd = 1.5, lty = "dashed")
lines(forecast::ma(dt$d15n.vs.air, order = 3, centre = TRUE) ~ dt$linear.ad, col = "purple", lwd = 2)

View(bulk.jack)
dt.j <- bulk.jack %>% select(d15n.vs.air, d13c.vs.vpdb, yearad.linear)

t.ssa <- ssa(dt.j$d15n.vs.air, L = 20)
plot(t.ssa)
plot(t.ssa, "series", groups = 1:20)
plot(t.ssa, "vectors", groups = 1:20)
plot(t.ssa, "paired", groups = 1:20)
plot(t.ssa, "wcor", groups = 1:20)
t.recon <- reconstruct(t.ssa, groups = list(1, c(2,3), 4, 5, 6, c(7,8)))
plot(t.recon, plot.method = "xyplot", type = "cumsum")

plot(d15n.vs.air ~ yearad.linear, data = dt.j, type = "l")
lines(t.recon$F1 +t.recon$F2 + t.recon$F3 ~ dt.j$yearad.linear, col = "blue", lwd = 2)
##################
## UNUSED CODE  ##
##################

# Below is old code to calculate the moaverage average (using the zoo package)
# First thing we do is apply a moving average
# bulk.stet$d15n.ra <- rollapply(bulk.stet$d15n.vs.air,
#                              width=5, # control the degree of smoothing
#                              align='center',
#                              FUN=mean,
#                              partial=TRUE)
# bulk.stet$d13c.ra <- rollapply(bulk.stet$d13c.vs.vpdb,
#                                width=5, # control the degree of smoothing
#                                align='center',
#                                FUN=mean,
#                                partial=TRUE)

# bulk.sav$d15n.ra <- rollapply(bulk.sav$d15n.vs.air,
#                               width = 1,
#                               align = 'center',
#                               FUN = mean,
#                               partial = TRUE)
# 
# bulk.sav$d13c.ra <- rollapply(bulk.sav$d13c.vs.vpdb,
#                               width = 1,
#                               align = 'center',
#                               FUN = mean,
#                               partial = TRUE)

# bulk.youngjack$d15n.ra <- rollapply(bulk.youngjack$d15n.vs.air,
#                               width = 1,
#                               align = 'center',
#                               FUN = mean,
#                               partial = TRUE)
# 
# bulk.youngjack$d13c.ra <- rollapply(bulk.youngjack$d13c.vs.vpdb, # d13C looks contaminated somehow with Jacksonville-4684
#                               width = 1,
#                               align = 'center',
#                               FUN = mean,
#                               partial = TRUE)
