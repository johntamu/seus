#'
#' John Schiff
#' 02/28/2019
#'
#'

library(dplyr)

######################################
## This is for all data analysis    ##
## regarding N-CSIAA data. I split  ##
## it into multiple parts.          ##
######################################

# Create pathnames -- Always to check if they conflict with other script files that are loaded

###############
# SEUS Region #
###############
path.win1 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff ncsiaa 09-01-2018.csv'
path.mac1 <- '~/Google Drive/projects/rproj/seus/data/schiff ncsiaa 02-06-2019.csv'
string <- '~/Documents/GitHub/data/schiff ncsiaa 05-28-2020.csv'

##################
# Gulf of Mexico #
##################
# This data is N-CSIAA data from Disks 1 and 7 from ONE black coral specimen in the Gulf of Mexico
# Disks 1 and 7 encompass different periods of time

path1 <- '~/Google Drive/projects/rproj/seus/data/prouty gom n-csiaa.csv'
path2 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/prouty gom n-csiaa.csv' # The eventual path for the Google Drive folder on the Windows computer
path3 <- '~/Documents/GitHub/data/prouty gom n-csiaa-means-v2.csv'

path4 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/prouty gom n-csiaa-means-v2.csv'
# Load the data
gom <- read.csv(path4)
gom <- read.csv(path3)

seus <- read.csv(string)

#####################################################
# POM from elsewhere (e.g., Western North Atlantic) #
#####################################################
pathwin <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff ncsiaa data from mompean 2016.csv'
ndata <- read.csv(pathwin)

############################################
# N-CSIA from other deep-sea coral studies #
############################################

pathdsc <- 'path goes here'

################################
## Calculate Trophic Position ##
################################

gom$trophic.position <- (((gom$Glu - gom$Phe) - 3.4)/7.6) + 1
seus$trophic.position <- (((seus$Glx - seus$Phe)-3.4)/7.6) + 1
ndata$trophic.position <- (((ndata$Glx - ndata$Phe) - 3.4)/7.6) + 1


####################
## Calculate SumV ## *************************************
####################
#########################################################################################################
## SumV = (1/7)*Sum(Abs(ChiAA))                                                                        ##
## (1/7) <-- assuming calculating from seven Trophic AAs (if 6, use 1/6 etc.; some pubs remove Valine) ##
## ChiAA <- deviation of each Trophic AA --> d15N-AA - Avg d15N-Trophic AA (from all seven)            ##
## Take the absolute value of each ChiAA                                                               ##
## Sum them all and multiply by (1/7)                                                                  ##
#########################################################################################################

#' Note: You could combine the datasets now and get that out of the way before calculations.
#' However I find that this method allows me more control.
#' 

##################
# Gulf of Mexico #
##################

trophic <- gom %>% select(X.1, Ala, Asp, Glu, Ile, Leu, Pro, Val)
avg.trophic <- rowMeans(trophic[, c(2:8)])
trophic$Dev.Ala <- abs(trophic$Ala - avg.trophic) # Remember absolute values!
trophic$Dev.Asp <- abs(trophic$Asp - avg.trophic)
trophic$Dev.Glu <- abs(trophic$Glu - avg.trophic)
trophic$Dev.Ile <- abs(trophic$Ile - avg.trophic)
trophic$Dev.Leu <- abs(trophic$Leu - avg.trophic)
trophic$Dev.Pro <- abs(trophic$Pro - avg.trophic)
trophic$Dev.Val <- abs(trophic$Val - avg.trophic)
trophic$Sum.Chi <- rowSums(trophic[, c(9:15)]) # Make sure you are selecting the correct columns, can probably use select() or filter() to do this

gom$Sum.V <- (1/7)*trophic$Sum.Chi

########
# SEUS #
########
trophic.seus <- seus %>% select(sample.id, Ala, Asx, Glx, Ile, Leu, Pro, Val)
avg.trophic2 <- rowMeans(trophic.seus[, c(2:8)])

trophic.seus$Dev.Ala <- abs(trophic.seus$Ala - avg.trophic2) # Remember absolute values!
trophic.seus$Dev.Asp <- abs(trophic.seus$Asp - avg.trophic2)
trophic.seus$Dev.Glu <- abs(trophic.seus$Glu - avg.trophic2)
trophic.seus$Dev.Ile <- abs(trophic.seus$Ile - avg.trophic2)
trophic.seus$Dev.Leu <- abs(trophic.seus$Leu - avg.trophic2)
trophic.seus$Dev.Pro <- abs(trophic.seus$Pro - avg.trophic2)
trophic.seus$Dev.Val <- abs(trophic.seus$Val - avg.trophic2)
trophic.seus$Sum.Chi <- rowSums(trophic.seus[, c(9:15)]) # Make sure you are selecting the correct columns
test <- rowSums(ndata.tr[, c(11:16)])
seus$Sum.V <- (1/7)*trophic.seus$Sum.Chi

#################
# Other Regions #
#################

ndata.tr <- ndata %>% select(Zone, Size.fraction, Ala, Asx, Glx, Ile, Leu, Val, Pro)
ndata.tr$Avg.Tr <- rowMeans(ndata.tr[, c(3:9)])
ndata.tr$Dev.Ala <- abs(ndata.tr$Ala - ndata.tr$Avg.Tr) # Remember absolute values!
ndata.tr$Dev.Asp <- abs(ndata.tr$Asx - ndata.tr$Avg.Tr)
ndata.tr$Dev.Glu <- abs(ndata.tr$Glx - ndata.tr$Avg.Tr)
ndata.tr$Dev.Ile <- abs(ndata.tr$Ile - ndata.tr$Avg.Tr)
ndata.tr$Dev.Leu <- abs(ndata.tr$Leu - ndata.tr$Avg.Tr)
ndata.tr$Dev.Pro <- abs(ndata.tr$Pro - ndata.tr$Avg.Tr)
ndata.tr$Dev.Val <- abs(ndata.tr$Val - ndata.tr$Avg.Tr)
ndata.tr$Sum.Chi <- rowSums(ndata.tr[, c(11:17)]) # Make sure you are selecting the correct columns, can probably use select() or filter() to do this
ndata$Sum.V <- (1/7)*ndata.tr$Sum.Chi

########
# SumV #
########

mean.sumv <- mean(gom$Sum.V) # Should be about 3
mean.sumv2 <- mean(seus$Sum.V) # Should be about 3.5 -- Use to check if you change the code
mean.sumv3 <- mean(ndata$Sum.V)

print(mean.sumv) 
print(mean.sumv2)
print(mean.sumv3)

# Now let's combine this into a new dataframe with GOM, SEUS, and now this POM data
# You can do this in R with some code (as I did earlier), but I'm lazy and did this in Excel
write.csv(seus, file = 'seus.csv')
write.csv(gom, file = 'gom.csv')
write.csv(ndata, file = 'ndata_all.csv')

# Read in the new, cleaned up file that has been combined with others
ndata <- read.csv("cleaned_ndata.csv")



##################################
## Principal Component Analysis ##
##################################
#' I haven't seen this used with
#' d15N-AA data yet so I have no references.
#' However, could be worth trying and seeing what
#' we get.
#' 
#' I am following the same idea with carbon:
#' Use it on Source-AA data. 
#' 

##########################
## Statistical Analyses ##
##########################

#' Identifying if differences 
#' between coral specimens are significant

ndata %>%
  filter(Sample.ID2 == "Jacksonville-4907" | Sample.ID2 == "Jacksonville-4684" | Sample.ID2 == "Savannah Banks-4902") -> corals

# seus %>%
#   filter(Sample.ID2 == "Jacksonville-4907" | Sample.ID2 == "Jacksonville-4684" | Sample.ID2 == "Savannah Banks-4902") -> corals

phe.aov <- aov(Phe ~ Sample.ID2, corals) # One-way ANOVA with Phe
sraa.aov <- aov(SrcAA ~ Sample.ID2, corals) # One-way ANOVA with avg Sr-AA
tp.aov <- aov(TP ~ Coral, corals)
sumv.aov <- aov(Sum.V ~ Sample.ID2, corals)

summary(phe.aov)
summary(sraa.aov)
summary(tp.aov)
summary(sumv.aov)


#' The results above tell us that there are significant 
#' differences between the groups, but not what 
#' they are

TukeyHSD(phe.aov)
TukeyHSD(sraa.aov)
TukeyHSD(tp.aov)
TukeyHSD(sumv.aov)

anc <- aov(Sum.V ~ Bulk*TP, corals)
summary(anc)

#####################
# Linear regression #
#####################

# d15N-Bulk vs AA parameters

reg <- lm(Sum.V ~ TP, corals)
summary(reg)

# Average Trophic AA for corals
corals$TrAA <- rowMeans(corals[,c(7:12)])
reg <- lm(TrAA ~ SrcAA, corals)
reg <- lm(Sum.V ~ SrcAA, corals)
reg <- lm(SrcAA ~ Phe, corals)
summary(reg)

################################
## Predicted bulk d15N        ##
## based on Trophic Position  ##
## in CSIAA samples           ##
################################

# Calculation: d15N (predicted) = (TP * d15N-source) + d15N-source or d15N-phe (try both)
nitrate <- 5 # d15N of nitrate based on depth profiles
src <- 4.5 # d15N of source nitrate or Phe depending on test you want to do
enrichment <- 3 # How much d15N enriches with each trophic transfer, in permil

test <- (1.9 * enrichment) + nitrate
print(test)

corals$d15N.pred <- (corals$TP * enrichment) + nitrate

ggplot(corals, aes(Year.CE, d15N.pred)) +
  geom_point() +
  theme_classic() +
  ylab(n) +
  xlab(x)

ggplot(corals, aes(Year.CE, Bulk)) +
  geom_point() +
  theme_classic() +
  ylab(n) +
  xlab(x)

ggplot(corals, aes(Year.CE, Nfix)) +
  geom_point() +
  theme_classic() +
  geom_vline(xintercept = 950, lty = 'dashed') +
  geom_vline(xintercept = 1250, lty = 'dashed') +
  geom_vline(xintercept = 1550, lty = 'dotted') +
  geom_vline(xintercept = 1850, lty = 'dotted') +
  ylab(n) +
  xlab(x)

ggplot(corals, aes(Year.CE, TP)) +
  geom_point() +
  theme_classic() +
  ylab("Trophic Position (Glu - Phe)") +
  xlab(x)

################################
## Calculated fraction Nfix   ##
## based on Bulk d15N         ##
## in CSIAA samples           ##
################################

d15N.Psink <- corals$Bulk
d15N.Nfix <- -1
d15N.nitrate <- 4.5

Nfix <- 1 - ((d15N.Psink - d15N.Nfix)/(d15N.nitrate - d15N.Nfix))
print(Nfix)
corals$Nfix <- Nfix

#' ------------------------------------------
#' 
#' Below is code to add the most recent ages
#' to the N-CSIAA data
#' 
#' ------------------------------------------
#' 

ndata <- read.csv("~/Documents/GitHub/data/cleaned_ndata.csv") # Contains SEUS black coral data, GOM black corals, and some POM data from elsewhere

# Add chronologies to N-CSIAA data

# Filter bulk dataframes that were loaded into the figures script file
df.sav %>%
  filter(sample.no. == 2 | sample.no. == 7 | sample.no. == 15 | sample.no. == 33 | sample.no. == 61 |
           sample.no. == 66 | sample.no. == 115) -> sav.n
df.jack %>%
  filter(sample.no. == 23 | sample.no. == 55 | sample.no. == 56 | 
           sample.no. == 134) -> jack.n
df.jack4684 %>%
  filter(sample.no. == 36 | sample.no. == 39 | sample.no. == 41 | sample.no. == 85 | 
           sample.no. == 155) -> jack4684.n

z <- rbind(sav.n, jack.n, jack4684.n)

ndata %>%
  filter(Region != "SEUS") -> t.df1

t.vector <- c(z$linear.ad, t.df1$Year.CE)

ndata$Year.CE <- t.vector # Now we have an updated ndata dataframe with current years
ndata$yrBP <- 1950 - ndata$Year.CE

write.csv(ndata, '~/Documents/GitHub/data/cleaned_ndata_12-22-2019.csv')
