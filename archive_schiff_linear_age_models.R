## ---
## Topic: Linear age models for black corals based on radiocarbon data
## John Schiff
## Date: October 22, 2018
# ---

###########################################################################
## This code is for developing robust age models based on deep-sea coral ##
## radiocarbon data. The age models are linear; non-linear/smooth spline ##
## age models can be done using the clam R package                       ##
###########################################################################

library(tidr)
library(dplyr)

path1 <- '~/Google Drive/projects/rproj/seus/data/schiff radiocarbon 02-05-2019.csv'
path2 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff radiocarbon 02-05-2019.csv'
radiocarbon <- read.csv(path2)

################################################
## 1. Linear age model for Jacksonville-4907  ##
##                                            ##
## This one is annoying because of the        ##
## variable growth rate in the black coral    ##
################################################

filter1 <- radiocarbon %>% filter(Coral == 'jack-4907-bc1-d1')
x1 <- filter1$Distance..um.
y1 <- filter1$mean

plot(y1 ~ x1)
cor(y1, x1)
linearMod1 <- lm(y1 ~ x1)
print(linearMod1)
summary(linearMod1)

plot(y1 ~ x1)
abline(linearMod1, col="blue")

t.ages <- read.delim('C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/clam/Cores/4907usgs/4907usgs_cleaned.txt', 
                     header = TRUE, sep = "\t", dec = ".") # 'Cleaned' file is created yourself
t.ages <- read.delim('~/Google Drive/projects/rproj/seus/data/clam/Cores/4907usgs/4907usgs_cleaned.txt', 
                     header = TRUE, sep = "\t", dec = ".")
filter1 <- cbind(filter1, t.ages)
filter1$mean <- rowMeans(filter1[c('yrmin', 'yrmax')], na.rm = TRUE)
filter1$mean <- round(filter1$mean)

# What I am going to do here is split it up into two linear growth rates
jack4907.split <- filter1 %>% select(Coral, ID.2, Distance..um., D14C, X14C.Age, error.1, yrmin, yrmax, mean, probability)
jack4907.split$distance <- jack4907.split$Distance..um.
jsplit1 <- jack4907.split %>% filter(distance > 2500)
jsplit2 <- jack4907.split %>% filter(distance < 2500)

jlm1 <- lm(mean ~ distance, jsplit1)
jlm2 <- lm(mean ~ distance, jsplit2)
summary(jlm1)
summary(jlm2)
print(jlm1$coefficients)
print(jlm2$coefficients)

plot(mean ~ distance, data = jack4907.split, tck = +0.015)
axis(side = 3, labels = FALSE, tck = +0.015)
axis(side = 4, labels = FALSE, tck = +0.015)
clip(0,2500,500,2200)
abline(jlm2)
clip(2000,8200,1750,3000)
abline(jlm1)

growthrate1 <- as.numeric(1/linearMod1$coefficients[2]) # Overall growth rate, in microns per year
g2 <- as.numeric(1/jlm2$coefficients[2]) # Growth rate at end of its life
g3 <- as.numeric(1/jlm1$coefficients[2]) # growth rate for majority of coral lifespan

print(g2)
print(g3)

###################################################
## 2. Linear age model for Savannah-4902         ##
##                                               ##
##                                               ##
##                                               ##
###################################################

filter2 <- radiocarbon %>% filter(Coral == 'sav-4902-bc1')
x2 <- filter2$Distance..um.
y2 <- filter2$mean
plot(y2 ~ x2)
cor(y2, x2)
linearMod2 <- lm(y2 ~ x2)
print(linearMod2)
summary(linearMod2)

plot(y2 ~ x2)
abline(linearMod2, col="blue")

t.ages <- read.delim('C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/clam/Cores/4902usgs/4902usgs_cleaned.txt', 
                     header = TRUE, sep = "\t", dec = ".") # 'Cleaned' file is created yourself
t.ages <- read.delim('~/Google Drive/projects/rproj/seus/data/clam/Cores/4902usgs/4902usgs_cleaned.txt', 
                     header = TRUE, sep = "\t", dec = ".")
filter2 <- cbind(filter2, t.ages)
filter2$mean <- rowMeans(filter2[c('yrmin', 'yrmax')], na.rm = TRUE)
filter2$mean <- round(filter2$mean)

savannah <- filter2 %>% select(Coral, ID.2, Distance..um., D14C, X14C.Age, error.1, yrmin, yrmax, mean, probability)
savannah$distance <- savannah$Distance..um. /1000

slm <- lm(mean ~ distance, data = savannah)
summary(slm)

new <- data.frame(distance = bulk.sav$distance..mm.)
m <- predict(slm, newdata = new, interval = "predict", level = 0.95, se.fit = TRUE) # se.fit = TRUE
m <- as.data.frame(m)
plot(mean ~ distance, data=savannah)
lines((m$fit.fit + m$se.fit) ~ bulk.sav$distance..mm.)
lines((m$fit.fit - m$se.fit) ~ bulk.sav$distance..mm.)

###################################################
## 3. Linear age model for Stetston-4904 (Beta)  ##
## from Beta Analytics                           ##
##                                               ##
##                                               ##
###################################################

filter3 <- radiocarbon %>% filter(Coral == 'stet-4904-bc1')
x3 <- filter3$Distance..um.
y3 <- filter3$mean
plot(y3 ~ x3)
cor(y3, x3)
linearMod3 <- lm(y3 ~ x3)
print(linearMod3)
summary(linearMod3)

plot(y3 ~ x3)
abline(linearMod3, col="blue")

t.ages <- read.delim('~/Google Drive/projects/rproj/seus/data/clam/Cores/stetson2/stetson2_cleaned.txt', 
                     header = TRUE, sep = "\t", dec = ".")
filter3 <- cbind(filter3, t.ages)
filter3$mean <- rowMeans(filter3[c('yrmin', 'yrmax')], na.rm = TRUE)
filter3$mean <- round(filter3$mean)

gr.stet2 <- as.numeric(1/linearMod3$coefficients[2])

######################################
## Ages based on linear growth rate ##
######################################
micron.depths <- stet.depths$depth * 1000
linear.ad <- 2005 - (micron.depths/gr.stet2)
linear.ad <- round(linear.ad)
bulk.stet$linear.ad2 <- linear.ad

# plot_residuals(linearMod1)

###################################################
## 4. Linear age model for Stetston-4904 (LLNL)  ##
##                                               ##
##                                               ##
##                                               ##
###################################################

filter4 <- radiocarbon %>% filter(Coral == 'stet-4904-bc1-d5')
# x4 <- filter4[-c(25:28),]$Distance..um.
# y4 <- filter4[-c(25:28),]$mean

x4 <- filter4$Distance..um.
y4 <- filter4$mean

plot(y4 ~ x4)
cor(y4, x4)
linearMod4 <- lm(y4 ~ x4)
print(linearMod4)
summary(linearMod4)

t.ages <- read.delim('C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/clam/Cores/stetson/stetson_cleaned.txt', 
                     header = TRUE, sep = "\t", dec = ".") # 'Cleaned' file is created yourself
t.ages <- read.delim('~/Google Drive/projects/rproj/seus/data/clam/Cores/stetson/stetson_cleaned.txt', 
                     header = TRUE, sep = "\t", dec = ".")

filter4 <- cbind(filter4, t.ages)
filter4$mean <- rowMeans(filter4[c('yrmin', 'yrmax')], na.rm = TRUE)
filter4$mean <- round(filter4$mean)

stetbanks <- filter4 %>% select(Coral, ID.2, Distance..um., D14C, X14C.Age, error.1, yrmin, yrmax, mean, probability)
stetbanks$distance <- stetbanks$Distance..um. /1000

stetlm <- lm(mean ~ distance, data = stetbanks)
summary(stetlm)

new <- data.frame(distance = c.stet1$distance..mm.)
n <- predict(stetlm, newdata = new, interval = "predict", level = 0.95, se.fit = TRUE) # se.fit = TRUE
n <- as.data.frame(n)
plot(mean ~ distance, data=stetbanks)
# lines(mean ~ distance, data=stetbanks)
# lines((n$fit.fit + n$se.fit) ~ c.stet1$distance..mm.)
# lines((n$fit.fit - n$se.fit) ~ c.stet1$distance..mm.)
lines((n$fit.fit + 156) ~ c.stet1$distance..mm.)
lines((n$fit.fit - 156) ~ c.stet1$distance..mm.)

plot_model(stetlm)

##############################
## Root mean squared error  ##
## (RMSE)                   ##
##############################

rmse <- sqrt(sum(residuals(lm)^2) / df.residual(lm))
print(rmse)

######################################
## Ages based on linear growth rate ##
######################################

micron.depths <- stet.depths$depth * 1000
linear.ad <- 2005 - (micron.depths/as.numeric(1/linearMod4$coefficients[2]))
linear.ad <- round(linear.ad)
bulk.stet$linear.ad <- linear.ad

#' 
#' Plotting a chart that shows the age ranges of all corals
#' as a box and whisker plot. Makes it easy to see
#' how they overlap.
#' 


###################################################
## Linear age model for Jacksonvile-4684 via a   ##
## bomb-spike model                              ##
##                                               ##
##                                               ##
###################################################

jc <- radiocarbon %>% filter(Coral == "jack-4684-bc1")
par(pty = "s")
plot(D14C ~ Distance..um., data = jc,
     type="o",
     xlim = c(0, 2500),
     ylim = c(-75, 150),
     ylab = expression(paste(Delta^{14}, "C")),
     xlab = expression(paste("Distance", " (",mu,"m)")))
lines(rise ~ Distance..um., data = rising, col = "#41b6c4", lwd = 2) # predicted based on lm.rising
lines(desc ~ Distance..um., data = descending, col = "#253494", lwd = 2) # predicted based on lm.descending
# abline(lm.rising, col = "#41b6c4", lwd = 2) # from lm() models below
# abline(lm.desc, col = "#253494", lwd = 2)
rising <- jc %>% filter(Distance..um. < 700 & Distance..um. > 269)
descending <- jc %>% filter(Distance..um. < 300)
lm.rising <- lm(D14C ~ Distance..um., data = rising)
lm.desc <- lm(D14C ~ Distance..um., data = descending)
rise <- predict(lm.rising, newdata = rising)
desc <- predict(lm.desc, newdata = descending)
tbl <- rbind(descending, rising)
bomb.table <- data.frame("Distance" = (tbl$Distance..um. / 1000),
                            "D14C" = tbl$D14C,
                            "Assigned.Age" = c(2004, NA, NA, 1975, 1975, NA, NA, 1957))
bomb.table <- bomb.table[-c(4),]
bomb.table
gr1 <- 270/(2004-1975)
gr2 <- (682-270)/(1975-1957)
gr3 <- 682/(2004-1957)
print(gr1)
print(gr2)
print(gr3)
mean(c(gr1, gr2, gr3))
sd(c(gr1, gr2, gr3))

##############################
## Lifespan vs Growth rate  ##
##############################

lifespan1 <- max(filter1$mean) - min(filter1$mean) # Jacksonville 4907
lifespan2 <- max(filter2$mean) - min(filter2$mean) # Savannah 4902
lifespan3 <- 2005 - (1950 - max(filter4[-c(2),]$mean)) # Stetson 4904, collected alive in 2005
lifespan4 

growthrate1 <- growthrate1
growthrate2 <- as.numeric(1/linearMod4$coefficients[2])
growthrate3 <- as.numeric(1/linearMod2$coefficients[2])

data <- data.frame(coral = c("Jack", "Stet", "Savannah"),
                   lifespan = c(lifespan1, lifespan3, lifespan2),
                   growth.rate = c(growthrate1, growthrate2, growthrate3))

plot(lifespan ~ growth.rate, data)
