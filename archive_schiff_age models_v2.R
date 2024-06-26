#' ------------
#' John Schiff
#' 02/26/2019
#' ------------

#' Importing radiocarbon data and 
#' doing linear growth rate 
#' calculations.
#' 
#' You will attach these chronologies
#' to the bulk records in a separate 
#' script file.
#' 

library(clam)
library(dplyr)

# Create pathnames
path1 <- '~/Documents/GitHub/data/schiff radiocarbon 02-05-2019.csv'
path2 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff radiocarbon 02-05-2019.csv'
path3 <- '/home/john/Desktop/data/schiff radiocarbon 02-05-2019.csv'

df <- read.csv(path1, header = TRUE)

radiocarbon <- df

##################################
## Create separate dataframes.  ##
##################################
r.jack <- radiocarbon %>% filter(Coral == 'jack-4907-bc1-d1')
r.jack2 <- radiocarbon %>% filter(Coral == 'jack-4684-bc1')
r.stet <- radiocarbon %>% filter(Coral == 'stet-4904-bc1-d5')
r.sav <- radiocarbon %>% filter(Coral == 'sav-4902-bc1')
r.stet2 <- radiocarbon %>% filter(Coral == 'stet-4904-bc1')

####################################
## Attach yrmin, yrmax estimates  ## 
## from clam calibration output   ##
####################################

# Jacksonville
ages1 <- read.delim('C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/clam/Cores/4907usgs/4907usgs_cleaned.txt', 
                     header = TRUE, sep = "\t", dec = ".") # 'Cleaned' file is created yourself
ages1 <- read.delim('~/Google Drive/projects/rproj/seus/data/clam/Cores/4907usgs/4907usgs_cleaned.txt', 
                     header = TRUE, sep = "\t", dec = ".")
r.jack <- cbind(r.jack, ages1)
r.jack$mean <- rowMeans(r.jack[c('yrmin', 'yrmax')], na.rm = TRUE)

# Stetson
ages2 <- read.delim('C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/clam/Cores/stetson/stetson_cleaned.txt', 
                     header = TRUE, sep = "\t", dec = ".") # 'Cleaned' file is created yourself
ages2 <- read.delim('~/Google Drive/projects/rproj/seus/data/clam/Cores/stetson/stetson_cleaned.txt', 
                     header = TRUE, sep = "\t", dec = ".")
ages_stet2 <- read.delim('C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/clam/Cores/stetson2/stetson2_cleaned.txt')

r.stet <- cbind(r.stet, ages2)
r.stet$mean <- rowMeans(r.stet[c('yrmin', 'yrmax')], na.rm = TRUE)

r.stet2 <- cbind(r.stet2, ages_stet2) 
r.stet2$mean <- rowMeans(r.stet2[c('yrmin', 'yrmax')], na.rm = TRUE)

# Savannah 
ages3 <- read.delim('C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/clam/Cores/4902usgs/4902usgs_cleaned.txt', 
                     header = TRUE, sep = "\t", dec = ".") # 'Cleaned' file is created yourself
ages3 <- read.delim('~/Google Drive/projects/rproj/seus/data/clam/Cores/4902usgs/4902usgs_cleaned.txt', 
                     header = TRUE, sep = "\t", dec = ".")
r.sav <- cbind(r.sav, ages3)
r.sav$mean <- rowMeans(r.sav[c('yrmin', 'yrmax')], na.rm = TRUE)

################################################
## Error ranges for each radiocarbon sample   ##
################################################
r.sav$range <- r.sav$yrmax - r.sav$yrmin
r.stet$range <- r.stet$yrmax - r.stet$yrmin
r.jack$range <- r.jack$yrmax - r.jack$yrmin


##################################
##
##
## CALCULATE LINEAR MODELS for  ##
## radiocarbon data sets        ##
##
##
##################################

# Jacksonville-4907
lm.jack <- lm(r.jack$mean ~ r.jack$Distance..um.) # The is the OVERALL linear model, not split into two
lm.jack <- lm(r.jack$X14C.Age ~ r.jack$Distance..um.)

#' --------------------------------------------
#' NOTE:
#' Jacksonville-4907 needs to be split into
#' TWO different dataframes and then use those
#' to come up with two different lienar models
#' and growth rates
#' --------------------------------------------
#' 

# Split into two dataframes. 
# jsplit1 <- r.jack %>% filter(Distance..um. < 2500)
# jsplit2 <- r.jack %>% filter(Distance..um. > 2500)
jsplit1 <- r.jack %>% filter(Distance..um. <= 1890) # test
jsplit2 <- r.jack %>% filter(Distance..um. >= 1890) # test

jsplit1 <- r.jack %>% filter(Distance.microns <= 1229) # official, 08/13/2019
jsplit2 <- r.jack %>% filter(Distance.microns >= 1229) # official, 08/13/2019

# Linear models for split data
jlm1 <- lm(mean ~ Distance.microns, jsplit1)
jlm2 <- lm(mean ~ Distance.microns, jsplit2)
jlm1.mm <- lm(mean ~ mm, jsplit1)
jlm2.mm <- lm(mean ~ mm, jsplit2)
summary(jlm1)
summary(jlm2)
summary(lm.jack)

resid <- resid(lm.jack)
plot(mean ~ Distance..um., r.jack, type = "o",
     xlab = expression(paste("Distance", " (",mu,"m)")),
     ylab = "Median Age (Yrs BP)",
     pch = 21,
     cex = 1.5,
     bg = "#bdbdbd")
axis(side = 3, labels = FALSE, tck = -0.015)
axis(side = 4, labels = FALSE, tck = -0.015)
# clip(2000,8200,1750,3000)
abline(jlm2, col = "#41b6c4", lty = "dashed", lwd = 2)
# clip(0,2500,500,2200)
abline(jlm1, col = "#253494", lty = "dashed", lwd = 2)
# abline(lm.jack, col = "#253494", lty = "dashed", lwd = 2)
abline(lm.jack, col = "black", lty = "dashed", lwd = 1.5)




# Stetson-4904
# Note: leaving some points out because of variability. Argue that it agrees with other Stetson radiocarbon set
# Update: 06/18/2019, deciding to leave those points in.

r.stet.edit <- r.stet[-c(25:28),]
lm.stet <- lm(mean ~ Distance..um., r.stet) # linear model for Stetson, edit: changed 08/09/2019
# lm.stet <- lm(r.stet$mean ~ r.stet$Distance..um.) # linear model for Stetson

# To NOT leave out the points, and do a regression model on all of them:
lm.stet.all <- lm(r.stet$mean ~ r.stet$Distance..um.) # linear model for Stetson
lm.stet.all <- lm(r.stet$X14C.Age ~ r.stet$Distance..um.) # linear model for Stetson

summary(lm.stet.all)

# Stetson-4904 
# Note: 07/12/2019, splitting the Stetson record into two grwoth rates like with Jacksonville to see how it varies
r.stet %>% 
  filter(Distance.microns <= 4998) -> stet.split1
r.stet %>%
  filter(Distance.microns >= 4998) -> stet.split2

View(stet.split1)
View(stet.split2)

plot(mean ~ Distance..um., stet.split1)
plot(mean ~ Distance..um., stet.split2)

s1 <- lm(mean ~ Distance.microns, stet.split1)
s2 <- lm(mean ~ Distance.microns, stet.split2)
s1.mm <- lm(mean ~ mm, stet.split1)
s2.mm <- lm(mean ~ mm, stet.split2)
summary(s1)
summary(s2)




# Savannah-4902
lm.sav <- lm(mean ~ Distance..um., r.sav) # linear model for Savannah
lm.sav <- lm(r.sav$X14C.Age ~ r.sav$Distance..um.) # linear model for Savannah
lm.sav.mm <- lm(mean ~ mm, r.sav)
summary(lm.sav)

################################
##                            ##
## PLOTTING linear age models ##
##                            ##
################################

################
# Jacksonville #
################

plot(mean ~ Distance.microns, r.jack, type = "o",
     xlab = expression(paste("Distance", " (",mu,"m)")),
     ylab = "Median Age (Yrs BP)",
     pch = 21,
     cex = 1.5,
     bg = "#bdbdbd")
axis(side = 3, labels = FALSE, tck = -0.015)
axis(side = 4, labels = FALSE, tck = -0.015)
# clip(2000,8200,1750,3000)
abline(jlm2, col = "#41b6c4", lty = "dashed", lwd = 2)
# clip(0,2500,500,2200)
abline(jlm1, col = "#253494", lty = "dashed", lwd = 2)
# abline(lm.jack, col = "#253494", lty = "dashed", lwd = 2)

#################
# Stetson Banks #
#################

plot(mean ~ Distance..um., r.stet, type = "o",
     xlab = expression(paste("Distance", " (",mu,"m)")),
     ylab = "Median Age (Yrs BP)",
     pch = 21,
     cex = 1.25,
     bg = "#bdbdbd")
# abline(lm.stet, col = "#000000", lty = "dashed", lwd = 1.25)
r2 <- summary(lm.stet.all)$adj.r.squared
# abline(lm.stet.all, col = "black", lwd = 1.25)
abline(lm.stet, col = "black", lty = "dashed")
mylabel <- bquote(italic(R)^2 == .(format(r2, digits = 3)))
text(x = 6000, y = 400, labels = mylabel)
legend("bottomright", bty = "n", legend = paste("R2 is", format(summary(lm.stet.all)$adj.r.squared, digits = 3)))
legend("bottomright", bty = "n", legend = bquote(italic(R)^2 == .(format(r2, digits = 3))))



plot(D14C ~ Distance..um., r.stet, type = "o",
     xlab = expression(paste("Distance", " (",mu,"m)")),
     ylab = "Median Age (Yrs BP)",
     pch = 21,
     cex = 1.5,
     bg = "#bdbdbd")

plot(D14C ~ Distance..um., r.stet2, type = "o",
     xlab = expression(paste("Distance", " (",mu,"m)")),
     ylab = "Median Age (Yrs BP)",
     pch = 21,
     cex = 1.5,
     bg = "#bdbdbd")
abline(lm(r.stet2$mean ~ r.stet2$Distance..um.), col = "#41b6c4", lty = "dashed", lwd = 2)

plot(D14C ~ Distance..um., r.jack, type = "o",
     xlab = expression(paste("Distance", " (",mu,"m)")),
     ylab = "D14C",
     pch = 21,
     cex = 1.5,
     bg = "#bdbdbd")

plot(D14C ~ Distance..um., r.jack2, type = "o",
     xlab = expression(paste("Distance", " (",mu,"m)")),
     ylab = "D14C",
     pch = 21,
     cex = 1.5,
     bg = "#bdbdbd")

plot(mean ~ Distance.microns, r.stet, type = "o",
     xlab = expression(paste("Distance", " (",mu,"m)")),
     ylab = "Median Age (Yrs BP)",
     pch = 21,
     cex = 1.25,
     bg = "#bdbdbd")
abline(lm.stet)
abline(s1, col = "black", lwd = 1.25)
abline(s2, col = "black", lwd = 1.25)

######################
#                 
#
# Savannah Banks    #
#
#
######################

plot(mean ~ Distance..um., r.sav, type = "o",
     xlab = expression(paste("Distance", " (",mu,"m)")),
     ylab = "Median Age (Yrs BP)",
     pch = 21,
     cex = 1.5,
     bg = "#bdbdbd")
abline(lm.sav, col = "#41b6c4", lty = "dashed", lwd = 2)

######################################
##                                  ##
##                                  ##
##                                  ##
## Calculating linear growth rates  ##
## Part 1                           ##
##                                  ##
##                                  ##
##                                  ##
######################################

#####################
#                   #
#                   #
# Jacksonville-4907 #
#                   #
#                   #
#####################

growth1 <- as.numeric(1/jlm1$coefficients[2]) # Slower growth rate
growth2 <- as.numeric(1/jlm2$coefficients[2]) # Faster growth rate
growth3 <- as.numeric(1/lm.jack$coefficients[2]) # Overall growth rate
print(growth1)
print(growth2)
print(growth3)

# Create a column of ages to attach to bulk record dataframe in separate script file
distance <- df.jack$distance*1000
yrs1 <- (1950-min(r.jack$mean)) - (distance/growth1) # 588 (minimum) is median year AD for the edge of the coral
yrs2 <- (1950-max(r.jack$mean)) - (distance/growth2) # 1380 is median year for the coral at distance 512 microns from edge
t.dat <- cbind(distance,yrs1,yrs2)
t.dat <- as.data.frame(t.dat)
t.dat$diff <- t.dat$yrs1 - t.dat$yrs2

first <- t.dat[1:67,]
first <- first[,-c(3)]
second <- t.dat[68:nrow(t.dat),]
second <- second[,-c(2)]
names(first)[2] <- paste("yrs")
names(second)[2] <- paste("yrs")
binded <- rbind(first,second)
length(binded$distance)

# Vector of linear ages from overall growth rate
yrs <- (1950-min(r.jack$mean)) - (distance/5)
test <- (1950-max(r.jack$mean)) + (distance/6.5)

# Below I am using the predict function to predict the age of the coral
t.data <- df.jack
t.data$Distance.microns <- t.data$distance*1000
t.data$predict1 <- predict(jlm1, t.data)
t.data$predict2 <- predict(jlm2, t.data)

t.data$predict1 <- 1950 - t.data$predict1 # to convert to years AD
t.data$predict2 <- 1950 - t.data$predict2 # to convert to years AD
t.data$diff <- t.data$predict1 - t.data$predict2 # find the inflection point, in thise case between ...

t.data %>%
  select(Distance.microns, predict1, predict2) -> t.data

first <- t.data[1:41,]
first <- first[,-c(3)]
second <- t.data[42:nrow(t.data),]
second <- second[,-c(2)]
names(first)[2] <- paste("yrs")
names(second)[2] <- paste("yrs")
binded <- rbind(first,second)
length(binded$Distance.microns) # check length that it is what you expect

#####################
#                   #
#                   #
# Jacksonville-4907 # Now do the same thing again with the other Jackonsville-4907 disk (Disk 1, at USGS in California)
#                   #
#                   #
#####################

distance2 <- jack2$distance..mm.*1000
yrs1 <- (1950-min(r.jack$mean)) - (distance2/growth1) # 588 is median year AD for the edge of the coral, g2 from linear_age_models.R script file
yrs2 <- (1950-1941) - (distance2/growth2) # 2062 is median year for the coral at distance 2558 microns from edge
t.dat <- cbind(distance2,yrs1,yrs2)
t.dat <- as.data.frame(t.dat)
t.dat$diff <- t.dat$yrs1 - t.dat$yrs2

first <- t.dat[1:21,]
first <- first[,-c(3)]
second <- t.dat[22:nrow(t.dat),]
second <- second[,-c(2)]
names(first)[2] <- paste("yrs")
names(second)[2] <- paste("yrs")

# k <- data.frame(distance2 = 540, yrs = ((480 + 303)/2)) # Slightly alter the method; this is the inflection point
# binded2 <- rbind(first, k, second)
binded2 <- rbind(first, second)
length(binded2$distance)

t.data <- jack2
t.data$Distance..um. <- t.data$distance..mm.*1000
t.data$predict1 <- predict(jlm1, t.data)
t.data$predict2 <- predict(jlm2, t.data)

t.data$predict1 <- 1950 - t.data$predict1
t.data$predict2 <- 1950 - t.data$predict2
t.data$diff <- t.data$predict1 - t.data$predict2 # find the inflection point, in thise case between 23 and 24

t.data %>%
  select(Distance..um., predict1, predict2) -> t.data

first <- t.data[1:23,]
first <- first[,-c(3)]
second <- t.data[24:nrow(t.data),]
second <- second[,-c(2)]
names(first)[2] <- paste("yrs")
names(second)[2] <- paste("yrs")
binded2 <- rbind(first,second)
length(binded2$Distance..um.)

#####################
#                   #
#                   #
# Stetson-4904      #
#                   #
#                   #
#####################

growth.stet1 <- as.numeric(1/lm.stet$coefficients[2]) # Rsquared = 0.94, residual standard error 94 years
growth.stet2 <- as.numeric(1/lm.stet.all$coefficients[2]) # Rsquared = 0.86, residual standard error 156 years

stet.linear.ad1 <- 2005 - ((stetson$distance..mm.*1000)/growth.stet1) # 'stetson' from timeseries script file
stet.linear.ad2 <- 2005 - ((stetson$distance..mm.*1000)/growth.stet2)

iodine.rate <- 11 # growth rate from the iodine SSA-MTM output, about 11.5 but has error.
stet.linear.ad4 <- 2005 - ((stetson$distance..mm.*1000)/iodine.rate)

print(growth.stet1)
print(growth.stet2)

test <- 2005 - ((stetson$distance..mm.*1000)/growth.stet)

# Update: 07/12/2019, growth rates for the split regression mdoels
s1.growth <- as.numeric(1/s1$coefficients[2])
s2.growth <- as.numeric(1/s2$coefficients[2])

print(s1.growth)
print(s2.growth)

lin1 <- -55 + ((df.stet$distance*1000)/s1.growth)
# lin2 <- 2005 - ((stetson$distance..mm.*1000)/s2.growth)
lin2 <- (-55 + 752.5) + (df.stet$distance*1000/s2.growth)


f <- cbind(stetson$distance..mm.*1000,lin1, lin2)
f <- as.data.frame(f)
f$diff <- f$lin1 - f$lin2 # this identifies the inflection point more easily, at distance = 8500 microns


first <- f[1:252,]
second <- f[253:nrow(f),]

first <- first[, -c(3)]
second <- second[,-c(2)]

names(first)[2] <- paste("yrs")
names(second)[2] <- paste("yrs")

stet.binded <- rbind(first, second)
length(stet.binded[,1])

stet.linear.ad3 <- 1950 - stet.binded$yrs


# predict the ages in Stetson based on linear regresion 8-14-2019 --> abandoned because it does not go to present (2005 AD)
t.stetson <- stetson
t.stetson$Distance.microns <- t.stetson$distance..mm.*1000
t.stetson$predict1 <- predict(s1, t.stetson)
t.stetson$predict2 <- predict(s2, t.stetson)

t.data$predict1 <- 1950 - t.data$predict1
t.data$predict2 <- 1950 - t.data$predict2
t.data$diff <- t.data$predict1 - t.data$predict2 # find the inflection point, in thise case between 23 and 24

t.data %>%
  select(Distance..um., predict1, predict2) -> t.data

first <- t.data[1:23,]
first <- first[,-c(3)]
second <- t.data[24:nrow(t.data),]
second <- second[,-c(2)]
names(first)[2] <- paste("yrs")
names(second)[2] <- paste("yrs")
binded2 <- rbind(first,second)
length(binded2$Distance..um.)

#################
# Savannah-4902 #
#################

growth.sav <- as.numeric(1/lm.sav$coefficients[2])
linear.ad2 <- (1950 - min(r.sav$mean)) - ((sav$distance..mm.*1000)/growth.sav)

t.data <- sav
t.data$Distance..um. <- t.data$distance..mm.*1000
vector <- predict(lm.sav, t.data)

sav$Distance..um. <- sav$distance..mm.*1000
sav$predict1 <- predict(lm.sav, sav)
sav.predict.ad <- 1950 - sav$predict1

######################################
##                                  ##
##                                  ##
##                                  ##
## Calculating linear growth rates  ##
## Part 2                           ##
##                                  ##
##                                  ##
##                                  ##
######################################

####################################
##                                ##
##                                ##
##                                ##
## BOMB SPIKE linear growth rates ##
## for Jacksonville-4684 coral    ##
##                                ##
##                                ##
##                                ##
####################################

rising <- r.jack2 %>% filter(Distance.microns < 700 & Distance.microns > 269)
descending <- r.jack2 %>% filter(Distance.microns < 300)
overall <- r.jack2 %>% filter(Distance.microns < 700)
whole <- r.jack2 %>% filter(Distance.microns > 545)

# Three linear models: rising, descending and overall
lm.rising <- lm(D14C ~ Distance.microns, data = rising)
lm.desc <- lm(D14C ~ Distance.microns, data = descending)
lm.overall <- lm(D14C ~ Distance.microns, data = overall)

lm.whole <- lm(X14C.Age ~ Distance.microns, data = whole)
# lm.whole.mm <- lm(X14C.Age ~ mm, data = whole)
# Use predict() function to create a line, NOTE: You don't need to do it this way, so I am commenting it out
# rise <- predict(lm.rising, newdata = rising)
# desc <- predict(lm.desc, newdata = descending)

t <- as.numeric(1/lm.whole$coefficients[2])


# Plot the distance vs. radiocarbon figure
par(pty = "s")
plot(D14C ~ Distance.microns, data = r.jack2,
     type="o",
     cex = 1.5,
     pch = 21,
     bg = '#bdbdbd',
     xlim = c(0, 2500),
     ylim = c(-75, 150),
     ylab = expression(paste(Delta^{14}, "C (\u2030)")),
     xlab = expression(paste("Distance", " (",mu,"m)")))
# lines(rise ~ Distance..um., data = rising, col = "#41b6c4", lwd = 2) # predicted based on lm.rising
# lines(desc ~ Distance..um., data = descending, col = "#253494", lwd = 2) # predicted based on lm.descending
clip(250, 750, -55, 140)
abline(lm.rising, col = "black", lwd = 1, lty = "dashed")
clip(0, 500, -50, 140)
abline(lm.desc, col = 'black', lwd = 1, lty = "dashed")


# Make a useful data table of bomb spike growth rates
tbl <- rbind(descending, rising)
bomb.table <- data.frame("Distance" = (tbl$Distance..um. / 1000),
                         "D14C" = tbl$D14C,
                         "Assigned.Age" = c(2004, NA, NA, 1975, 1975, NA, NA, 1957))
bomb.table <- bomb.table[-c(4),]
bomb.table




#' Figure below is for Results section. It is
#' a linear regression of each radiocarbon age model
#' with Rsquared value.
#' 

par(pty = "s", mfrow = c(1,3)) # edit to adjust for the number of radiocarbon age models

# Stetson-4904 BC1
plot(mean ~ Distance.microns, r.stet, type = "o",
     xlab = expression(paste("Distance", " (",mu,"m)")),
     ylab = "Years BP",
     pch = 21,
     cex = 1.25,
     bg = "#bdbdbd")
# abline(lm.stet, col = "#000000", lty = "dashed", lwd = 1.25)
r2 <- summary(lm.stet.all)$adj.r.squared
abline(lm.stet.all, col = "black", lwd = 1.25)
mylabel <- bquote(italic(R)^2 == .(format(r2, digits = 3)))
text(x = 6000, y = 400, labels = mylabel)
# legend("bottomright", bty = "n", legend = paste("R2 is", format(summary(lm.stet.all)$adj.r.squared, digits = 3)))
legend("bottomright", bty = "n", legend = bquote(italic(R)^2 == .(format(r2, digits = 3))))

# Jacksonville-4907 BC1
plot(mean ~ Distance..um., r.jack, type = "o",
     xlab = expression(paste("Distance", " (",mu,"m)")),
     ylab = "Years BP",
     pch = 21,
     cex = 1.25,
     bg = "#bdbdbd")
r2j <- summary(lm.jack)$adj.r.squared
abline(lm.jack, col = "black", lwd = 1.25)
mylabel <- bquote(italic(R)^2 == .(format(r2, digits = 3)))
# text(x = 6000, y = 400, labels = mylabel)
# legend("bottomright", bty = "n", legend = paste("R2 is", format(summary(lm.stet.all)$adj.r.squared, digits = 3)))
legend("bottomright", bty = "n", legend = bquote(italic(R)^2 == .(format(r2j, digits = 3))))

# Savannah
plot(mean ~ Distance..um., r.sav, type = "o",
     xlab = expression(paste("Distance", " (",mu,"m)")),
     ylab = "Years BP",
     pch = 21,
     cex = 1.25,
     bg = "#bdbdbd")
r2s <- summary(lm.sav)$adj.r.squared
abline(lm.sav, col = "black", lwd = 1.25)
mylabel <- bquote(italic(R)^2 == .(format(r2, digits = 3)))
# text(x = 6000, y = 400, labels = mylabel)
# legend("bottomright", bty = "n", legend = paste("R2 is", format(summary(lm.sav)$adj.r.squared, digits = 3)))
legend("bottomright", bty = "n", legend = bquote(italic(R)^2 == .(format(r2s, digits = 3))))



################
# Growth rates #
################

gr1 <- 270/(2004-1975)
gr2 <- (682-270)/(1975-1957) # 23 microns per year
gr3 <- 682/(2004-1957) # 15 microns per year
gr4 <- as.numeric(1/lm.whole$coefficients[2])

print(gr1)
print(gr2)
print(gr3)
print(gr4)

m <- mean(c(gr1, gr2, gr3))
sd(c(gr1, gr2, gr3))

ad1 <- 2004 - ((jack4684$distance..mm.*1000)/gr3) # 15 microns per year, from rounding up average overall growth rate
ad2 <- 2004 - ((jack4684$distance..mm.*1000)/gr2) # 23 microns per year, this also fits with overall linear trend
ad3 <- 2004 - ((jack4684$distance..mm.*1000)/gr1) # 9 microns per year
ad4 <- 2004 - ((jack4684$distance..mm.*1000)/20) # Test, based on sheet from Nancy
ad5 <- 2004 - ((jack4684$distance..mm.*1000)/16) # from 'm'

t.df <- cbind(ad1, ad2, ad3, ad4, ad5)
t.df <- as.data.frame(t.df)

########################################
##                                    ##
## SECTION 2: calibrations using the  ##
## clam package.                      ##
##                                    ##
########################################

# Assign a folder where you will store your age model output

core.dir1 <- "~/Google Drive/projects/rproj/seus/data/clam/Cores"
core.dir2 <- "C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/clam/Cores"

########################
## Useful data tables ##
########################

####################
## Miscellaneous  ##
####################

##############################
## Root mean squared error  ##
## (RMSE)                   ##
##############################

rmse <- sqrt(sum(residuals(lm)^2) / df.residual(lm))
print(rmse)


