## ---
## Topic: Gulf of Mexico black coral N-CSIAA data, used in Prouty et al 2015
## John Schiff
## Date: January 21, 2019
## ---

##########################################################
## This code will be used to analyze the Gulf of Mexico ##
## black coral N-CSIAA data. I want to see how it       ##
## compares to my Gulf Stream black coral data.         ##
##########################################################

#' Note: In the future, work to create functions that do some of the
#' calculations for us. For example, Trophic Position, SumV, etc.

# First load approrpriate libraries and load the data
library(ggplot2)
library(lattice)
library(reshape2)
library(zoo)
library(dplyr)
# library(tidyr)

# Load the GOM black coral N-CSIAA data
path1 <- '~/Google Drive/projects/rproj/seus/data/prouty gom n-csiaa.csv'
path2 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/prouty gom n-csiaa.csv' # The eventual path for the Google Drive folder on the Windows computer
path3 <- '~/Google Drive/projects/rproj/seus/data/prouty gom n-csiaa-means-v2.csv'
path4 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/prouty gom n-csiaa-means-v2.csv'
# gom <- read.csv(path2)
gom <- read.csv(path4)
gom <- read.csv(path3)
# gom$year.ad <- 


# This data is N-CSIAA data from Disks 1 and 7 from ONE black coral specimen in the Gulf of Mexico
# Disks 1 and 7 encompass different periods of time

######################### 
## Trophic Position    ##
#########################

gom$trophic.position <- (((gom$Glu - gom$Phe) - 3.4)/7.6) + 1
# gom$Disk <- NA


#########################################################################################################
## Coming up with a way to quickly calculate SumV across a lot of different samples                    ##
## Each sample will have d15N data from Trophic AA: Ala, Asp, Glu, Ile, Leu, Pro, Val)                 ##
## SumV = (1/7)*Sum(Abs(ChiAA))                                                                        ##
## (1/7) <-- assuming calculating from seven Trophic AAs (if 6, use 1/6 etc.; some pubs remove Valine) ##
## ChiAA <- deviation of each Trophic AA --> d15N-AA - Avg d15N-Trophic AA (from all seven)            ##
## Take the absolute value of each ChiAA                                                               ##
## Sum them all and multiply by (1/7)                                                                  ##
#########################################################################################################

# First, remove all non-Trophic AAs
# gom.tr <- gom %>% select((!!as.name(column)) == 'Ala' | 'Asp' | 'Glu' | 'Ile' | 'Leu' |
gom.tr <- gom %>% select(X.1, Ala, Asp, Glu, Ile, Leu, Pro, Val)
gom.tr$Avg.Tr <- rowMeans(gom.tr[, c(2:8)])
# Make a column of Ala - Avg Tr and go from there!
gom.tr$Dev.Ala <- abs(gom.tr$Ala - gom.tr$Avg.Tr) # Remember absolute values!
gom.tr$Dev.Asp <- abs(gom.tr$Asp - gom.tr$Avg.Tr)
gom.tr$Dev.Glu <- abs(gom.tr$Glu - gom.tr$Avg.Tr)
gom.tr$Dev.Ile <- abs(gom.tr$Ile - gom.tr$Avg.Tr)
gom.tr$Dev.Leu <- abs(gom.tr$Leu - gom.tr$Avg.Tr)
gom.tr$Dev.Pro <- abs(gom.tr$Pro - gom.tr$Avg.Tr)
gom.tr$Dev.Val <- abs(gom.tr$Val - gom.tr$Avg.Tr)
gom.tr$Sum.Chi <- rowSums(gom.tr[, c(10:16)]) # Make sure you are selecting the correct columns, can probably use select() or filter() to do this
gom.tr$Sum.V <- (1/7)*gom.tr$Sum.Chi


#############################################################################
## Below is code for combining the GOM and Gulf Stream N-CSIAA             ##
## datasets to do appropriate comparisons. For now, I am importing         ##
## data from the CSIAA csv file referenced in the schiff_seus_timeseries.R ##
## script file. Always make sure to check the original CSV file.           ##
#############################################################################

path.win1 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff ncsiaa 09-01-2018.csv'
path.mac1 <- '~/Google Drive/projects/rproj/seus/data/schiff ncsiaa 02-06-2019.csv'
csiaa <- read.csv(path.win1)
csiaa <- read.csv(path.mac1)
csiaa$year.ad <- (1950-csiaa$year.bp.clam)
csiaa.tr <- csiaa %>% select(sample.id, Ala, Asx, Glx, Ile, Leu, Pro, Val)
csiaa.tr$Avg.Tr <- rowMeans(csiaa.tr[, c(2:8)])
csiaa$year.ad <- (1950-csiaa$calc.yr.bp)

csiaa.tr$Dev.Ala <- abs(csiaa.tr$Ala - csiaa.tr$Avg.Tr) # Remember absolute values!
csiaa.tr$Dev.Asp <- abs(csiaa.tr$Asx - csiaa.tr$Avg.Tr)
csiaa.tr$Dev.Glu <- abs(csiaa.tr$Glx - csiaa.tr$Avg.Tr)
csiaa.tr$Dev.Ile <- abs(csiaa.tr$Ile - csiaa.tr$Avg.Tr)
csiaa.tr$Dev.Leu <- abs(csiaa.tr$Leu - csiaa.tr$Avg.Tr)
csiaa.tr$Dev.Pro <- abs(csiaa.tr$Pro - csiaa.tr$Avg.Tr)
csiaa.tr$Dev.Val <- abs(csiaa.tr$Val - csiaa.tr$Avg.Tr)
csiaa.tr$Sum.Chi <- rowSums(csiaa.tr[, c(10:16)]) # Make sure you are selecting the correct columns
csiaa.tr$Sum.V <- (1/7)*csiaa.tr$Sum.Chi

# Mean SumV
sumv.gom <- mean(gom.tr$Sum.V) # Should be about 3
sumv.gs <- mean(csiaa.tr$Sum.V) # Should be about 3.5 -- Use to check if you change the code

print(sumv.gom)
print(sumv.gs)

# Calculate Trophic Position

gom.tr$Trophic.Position <- (((gom$Glu - gom$Phe)-3.4)/7.6) + 1
csiaa.tr$Trophic.Position <-  (((csiaa$Glx - csiaa$Phe)-3.4)/7.6) + 1

# We want to use Mol% in some calculations. We have this information from the Stetson coral


##########################
## Trophic Dynamics     ##
##########################

# disk1.tr <- mean(gom.tr) 
# gom.tr$disk7 

# xyplot(gom$Glu ~ gom$Avg.Tr,
#        group = gom$Disk,
#        auto.key = TRUE)
# 
# xyplot(Glu ~ Avg.Tr,
#        csiaa.tr,
#        group = coral.id,
#        auto.key = TRUE,
#        pch = 19, cex = 1.5)
# 
# # Merge Gulf of Mexico and Gulf Stream dataframes together
# # rename sample ID columns so you can do this
# 
# csiaa.tr <- csiaa.tr %>%
#   rename(
#     Sample.ID = coral.id
#   )
#   
# gom.tr <- gom.tr %>%
#   rename(
#     Sample.ID = X.1
#   )
# 
# total <- rbind(csiaa.tr, gom.tr) # Easy!
# list1 <- 1:16
# list2 <- 17:30
# 
# gs <- rep("Gulf Stream", length(list1))
# gmex <- rep("Gulf of Mexico", length(list2))
# # list3 <- data.frame(gs, gmex)
# 
# total$Region <- c(gs, gmex)
# 
# xyplot(Glu ~ Avg.Tr | Region,
#        total,
#        group = Sample.ID,
#        pch = 19,
#        cex = 1.5)
# 
# # Make new ID column
# savannah <- rep("Savannah Banks-4902", length(1:7))
# jacksonville <- rep("Jacksonville-4907", length(8:11))
# jacksonville2 <- rep("Jacksonville-4684", length(12:16))
# gulfmexico <- rep("Gulf of Mexico-Disk 7", length(17:20))
# gulfmexico2 <- rep("Gulf of Mexico-Disk 1", length(21:29))
# 
# total$Sample.ID2 <- c(savannah, jacksonville, jacksonville2, gulfmexico, gulfmexico2, "POM")
# 
# xyplot(Glu ~ Avg.Tr | Region,
#        total,
#        group = Sample.ID2,
#        pch = 19,
#        cex = 1.5,
#        auto.key = TRUE,
#        layout = c(2,1))


##############################################################
## We are going to do the same process as above             ##  
## (combining data frames, etc.) but with the main          ##
## imported data sets (i.e., before doing the "csiaa.tr"    ##
## or "gom.tr" step in the code above. I think one way      ##
## to do this would be to just add the Source AAs back into ##
## the "csiaa.tr" and "gom.tr" files and rename them as new ##
## dataframes.                                              ##
##############################################################

csiaa.all <- csiaa.tr
gom.all <- gom.tr

# Add in Source-AAs
csiaa.all$Phe <- csiaa$Phe
csiaa.all$Gly <- csiaa$Gly
csiaa.all$Ser <- csiaa$Ser
csiaa.all$Tyr <- csiaa$Tyr
csiaa.all$Lys <- csiaa$Lys
csiaa.all$Avg.Src <- ((csiaa$Gly + csiaa$Ser + csiaa$Phe)/3)

gom.all$Phe <- gom$Phe
gom.all$Gly <- gom$Gly
gom.all$Ser <- gom$Ser
gom.all$Tyr <- gom$Tyr
gom.all$Lys <- gom$Lys
gom.all$Avg.Src <- ((gom$Gly + gom$Ser + gom$Phe)/3)

csiaa.all <- csiaa.all %>%
  rename(
    Sample.ID = sample.id
  )

gom.all <- gom.all %>%
  rename(
    Sample.ID = X.1
  )

seus <- rbind(csiaa.all, gom.all)

gs <- rep("Gulf Stream", length(1:16))
gmex <- rep("Gulf of Mexico", length(17:30))

savannah <- rep("Savannah Banks-4902", length(1:7))
jacksonville <- rep("Jacksonville-4907", length(8:11))
jacksonville2 <- rep("Jacksonville-4684", length(12:16))
gulfmexico <- rep("Gulf of Mexico-Disk 7", length(17:20))
gulfmexico2 <- rep("Gulf of Mexico-Disk 1", length(21:29))

seus$Sample.ID2 <- c(savannah, jacksonville, jacksonville2, gulfmexico, gulfmexico2, "POM")

seus$Region <- c(gs, gmex)

xyplot(Trophic.Position ~ Sum.V | Region,
       seus,
       group = Sample.ID2,
       pch = 21,
       cex = 1.5,
       auto.key=list(columns=2, cex = 0.75),
       xlab = expression(paste(Sigma,"V")),
       ylab = "Trophic Position (Glu - Phe)")

# xyplot(Avg.Src ~ Avg.Tr | Region,
#        t.seus,
#        group = Sample.ID2,
#        pch = 21,
#        cex = 1.5,
#        auto.key=TRUE)

xyplot(Glu ~ Avg.Tr | Region,
       seus,
       group = Sample.ID2,
       pch = 21,
       cex = 1.5,
       auto.key=TRUE)

# Put a timeframe to each sample if we want to continue doing temporal analyses
f <- rep(NA, length(17:30))
seus$Year.AD <- c(csiaa$year.ad, f)
ff <- c(356, 312, 269, 225, 2003, 2000, 1997, 1995, 1990, 1985, 1980, 1975, 1969, NA)
seus$Year.AD <- c(csiaa$year.ad, ff)

# Add bulk d15N to the seus dataset
seus$Bulk.N <- c(csiaa$bulk.N, f)

####################################
##      Temporal analysis         ##
##  of NCSIA from Gulf of Mexico  ##
##  and Gulf Stream               ##
####################################

xyplot(Phe ~ Year.AD | Region,
       seus,
       group = Sample.ID2,
       pch = 21,
       cex = 1.5,
       auto.key=TRUE)

xyplot(Sum.V ~ Year.AD | Region,
       seus,
       group = Sample.ID2,
       pch = 21,
       cex = 1.5,
       auto.key=TRUE)

xyplot(Trophic.Position ~ Year.AD | Region,
       seus,
       group = Sample.ID2,
       pch = 21,
       cex = 1.5,
       auto.key=TRUE)

xyplot(Avg.Src ~ Year.AD | Region,
       seus,
       group = Sample.ID2,
       pch = 21,
       cex = 1.5,
       auto.key=TRUE)

xyplot(Glu ~ Year.AD | Region,
       seus,
       group = Sample.ID2,
       pch = 21,
       cex = 1.5,
       auto.key=TRUE)

xyplot(Sum.V ~ Year.AD,
       seus,
       group = Sample.ID2,
       auto.key = list(columns=3, cex = 0.5),
       ylab = expression(paste(Sigma,"V")),
       xlab = "Calendar Year (C.E.)",
       pch = 21,
       cex = 1.5)
       # panel.xyplot(Sum.V ~ Year.AD,
       #              abline=c(0,500)))


# panel = function(x, y...) {
#   panel.abline(v=500)
#   panel.abline(v=600)
# },

xyplot(Bulk.N ~ Sum.V,
       seus)
xyplot(Sum.V ~ Avg.Src,
       seus)

linearMod <- lm(Bulk.N ~ Sum.V, data=seus)
linearMod <- lm(Trophic.Position ~ Sum.V, data=seus)
linearMod <- lm(Phe ~ Avg.Src, data=seus)
linearMod <- lm(Gly ~ Avg.Tr, data=seus)
linearMod <- lm(Phe ~ Trophic.Position, data=seus)
linearMod <- lm(Sum.V ~ Bulk.N + Phe + Gly, data=seus)

summary(linearMod)
plot(linearMod)

linearMod <- lm(Bulk.N ~ Trophic.Position, data=seus)
linearMod <- lm(Bulk.N ~ Avg.Src, data=seus)
linearMod <- lm(Bulk.N ~ Avg.Tr, data=seus)
linearMod <- lm(Bulk.N ~ Gly + Ser + Phe + Avg.Src + Sum.V,
          data = seus)

xyplot(Phe ~ Trophic.Position,
       seus,
       ylab = expression({delta}^15*"N Phe"),
       xlab = "Trophic Position",
       pch = 21,
       cex = 1.5,
       col = "orange")


##################################################################
## This code is so we can categorize our amino acids            ##
## by Source, Trophic, Essential, Non-Essential, since          ##
## this is not in the original CSV datasheet (though you could  ##
## also just edit that).                                        ##
##################################################################
library(data.table) # Need to install this package 

seus.melt <- melt(seus, id.vars=c('Sample.ID', 'Sample.ID2'))
seus.melt <- data.table(seus.melt)
# v[, Group.1 := ifelse(variable == c('Phe','Gly','Ser','Lys','Tyr'), "Source", NA)]
# v[, Group.1 := ifelse(variable %in% c('Phe','Gly','Ser','Lys','Tyr'), "Source", NA)]
# v[, Group.1 := ifelse(variable %in% c('Glu','Asp','Ala','Ile','Leu','Pro','Val'), "Trophic", NA)]
seus.melt[, Group.1 := ifelse(variable %in% c('Phe','Gly','Ser','Lys','Tyr'), "Source",
                       ifelse(variable %in% c('Glu','Asp','Ala','Ile','Leu','Pro','Val'), "Trophic", NA))]
seus.melt[, Group.2 := ifelse(variable %in% c('Phe', 'Thr', 'Ile', 'Leu','Val','Lys'), "Essential",
                       ifelse(variable %in% c('Asp', 'Glu', 'Pro', 'Ala', 'Ser', 'Gly'), "Non-Essential", NA))]


# seus.melt <- data.frame(lapply(seus.melt, as.character), stringsAsFactors=FALSE)

seus.melt <- seus.melt %>% filter(variable %in% c('Phe','Gly','Ser','Lys','Tyr','Glu','Asp','Ala','Ile','Leu','Pro','Val'))

seus.melt$variable1 <- factor(seus.melt$variable, levels=c('Glu', 'Asp', 'Ala', 'Ile', 'Leu', 'Pro', 'Val', 'Gly', 'Ser', 'Lys', 'Tyr', 'Phe')) # For Trophic/Source AA ordering

xyplot(value ~ variable1,
       seus.melt,
       panel = function(x,  y, ...){
         panel.xyplot(x, y, ...)
         panel.abline(h = 9, lty = 1)
         panel.abline(h = 9.25, lty = 2)
         panel.abline(h = 8.75, lty = 2)
       },
       pch=21,
       cex=1.25,
       ylim=c(-20, 40), # ylim and xlim are "first class" parameters and don't need to be in scales=list()
       group = Sample.ID2,
       xlab = NULL,
       ylab = expression({delta}^15*"N (\u2030)"))

plot(value ~ variable, data=seus.melt)


##############################################################                     
## Calculating and Plotting amino acid means                ##
## For each variable (coral specimen, POM, etc.).           ##
## This is good if you want to show what the average d15N   ##
## for each AA is for a given coral specimen with multiple  ##
## peels used for CSIA.                                     ##
##############################################################  

# Remember that 'seus' has Sample.ID2 which is by specimen (or disk, etc.) instead ofby peel
# Sample.ID2 is essentially some sort of secondary characteristic of the sample, such as region, Disk, location, etc.

aa <- c('Sample.ID2','Phe','Gly','Ser','Lys','Tyr','Glu','Asp','Ala','Ile','Leu','Pro','Val') # Don't need Sample.ID1

seus2 <- seus %>% select(aa) # New dataframe with just the Sample.ID1 and AAs

# One example online had a manual way of doing by using colMeans() function, but the way below might be simpler

seus.means <- aggregate(seus2[,3:14], list(Sample.ID2=seus2$Sample.ID2), mean, na.rm = TRUE) # Using aggregate() function for this, but make sure to ignore NAs
# In the above code, list(Sample.ID2 = ...) is just so that the first column doesn't have the name "Group.1"
# Otherwise it can just be list(seus2$Sample.ID2)

seus.means.melt <- melt(seus.means, "Sample.ID2")
seus.means.melt$variable1 <- factor(seus.means.melt$variable, 
                              levels=c('Glu', 'Asp', 'Ala', 'Ile', 'Leu', 'Pro', 'Val', 'Gly', 'Ser', 'Lys', 'Tyr', 'Phe')) # For Trophic/Source AA ordering

##############
## Plotting ##
##############

xyplot(value ~ variable1,
       seus.means.melt,
       panel = function(x,  y, ...){
         panel.xyplot(x, y, ...)
         panel.abline(h = 9, lty = 1)
         panel.abline(h = 9.25, lty = 2)
         panel.abline(h = 8.75, lty = 2)
       },
       pch=21,
       cex=1.5,
       ylim=c(-20, 40), # ylim and xlim are "first class" parameters and don't need to be in scales=list()
       group = Sample.ID2,
       xlab = NULL,
       ylab = expression({delta}^15*"N (\u2030)"),
       auto.key=list(columns=2, cex = 0.75))

bwtheme <- standard.theme("pdf", color = FALSE)
xyplot(value ~ variable1,
       seus.means.melt,
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
       par.settings = bwtheme)



xyplot(Phe ~ Avg.Src | Region,
       seus[1:29,],
       panel = function(x,  y, ...){
         panel.xyplot(x, y, ...)
         panel.lmline(x, y,
                      lty=2)
       },
       layout = c(2,1),
       par.settings = bwtheme,
       cex = 1.5,
       ylab = expression({delta}^15*"N Phe (\u2030)"),
       xlab = expression({delta}^15*"N Source"['AA']*'(\u2030)'))


################################################
## Now we are going to import code from other ## 
## published studies for a data comparison    ##
################################################

pathwin <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff ncsiaa data from mompean 2016.csv'
# pathmac <- 
ndata <- read.csv(pathwin)
ndata.tr <- ndata %>% select(Zone, Size.fraction, Ala, Asx, Glx, Ile, Leu, Val, Pro)
ndata.tr$Avg.Tr <- rowMeans(ndata.tr[, c(3:9)])
# Make a column of Ala - Avg Tr and go from there!
ndata.tr$Dev.Ala <- abs(ndata.tr$Ala - ndata.tr$Avg.Tr) # Remember absolute values!
ndata.tr$Dev.Asp <- abs(ndata.tr$Asx - ndata.tr$Avg.Tr)
ndata.tr$Dev.Glu <- abs(ndata.tr$Glx - ndata.tr$Avg.Tr)
ndata.tr$Dev.Ile <- abs(ndata.tr$Ile - ndata.tr$Avg.Tr)
ndata.tr$Dev.Leu <- abs(ndata.tr$Leu - ndata.tr$Avg.Tr)
ndata.tr$Dev.Pro <- abs(ndata.tr$Pro - ndata.tr$Avg.Tr)
ndata.tr$Dev.Val <- abs(ndata.tr$Val - ndata.tr$Avg.Tr)
ndata.tr$Sum.Chi <- rowSums(ndata.tr[, c(11:17)]) # Make sure you are selecting the correct columns, can probably use select() or filter() to do this
ndata.tr$Sum.V <- (1/7)*ndata.tr$Sum.Chi
ndata.all <- ndata.tr
ndata.all$Phe <- ndata$Phe
ndata.all$Gly <- ndata$Gly
ndata.all$Ser <- ndata$Ser
ndata.all$Tyr <- ndata$Tyr
ndata.all$Lys <- ndata$Lys
ndata.all$Lys <- ndata$Thr # Not actually a Source AA, but in the "metabolic" group
ndata.all$Bulk <- ndata$Bulk
ndata.all$Avg.Src <- ((ndata$Gly + ndata$Ser + ndata$Phe)/3)

# Now let's combine this into a new dataframe with GOM, SEUS, and now this POM data
# You can do this in R with some code (as I did earlier), but I'm lazy and did this in Excel
write.csv(ndata.all, file = 'ndata_all.csv')

# Now, after cleaning up in Excel and combining with the other two
globe <- read.csv("global_ncsiaa.csv") # Keeping the new filename different so I don't accidentally overwrite it when running old code
globe$Trophic.Position <- (((globe$Glx - globe$Phe)-3.4)/7.6) + 1
globe$Avg.Tr <- rowMeans(globe[, c(2:8)])

globe$Dev.Ala <- abs(globe$Ala - globe$Avg.Tr)
globe$Dev.Asx <- abs(globe$Asx - globe$Avg.Tr)
globe$Dev.Glx <- abs(globe$Glx - globe$Avg.Tr)
globe$Dev.Ile <- abs(globe$Ile - globe$Avg.Tr)
globe$Dev.Leu <- abs(globe$Leu - globe$Avg.Tr)
globe$Dev.Pro <- abs(globe$Pro - globe$Avg.Tr)
globe$Dev.Val <- abs(globe$Val - globe$Avg.Tr)
globe$Sum.Chi <- rowSums(globe[, c(10:16)]) # Make sure you are selecting the correct columns
globe$Sum.V <- (1/7)*globe$Sum.Chi

t.globe <- globe %>% dplyr::filter(Region != "Eastern North Atlantic" & Region != "Central North Atlantic")
# t.globe <- 

globe.melt <- melt(globe[1:45,], "Sample.ID2")
globe.melt$variable1 <- factor(globe.melt$variable, 
                                    levels=c('Glu', 'Asp', 'Ala', 'Ile', 'Leu', 'Pro', 'Val', 'Gly', 'Ser', 'Lys', 'Tyr', 'Phe', 'Thr')) # For Trophic/Source AA ordering


##############
## Plotting ##
##############

myColors <- brewer.pal(9, "YlGnBu")
my.settings <- list(
  superpose.symbol=list(pch = 21,
                        col = "black",
                        fill=myColors[1:8]),
  strip.background=list(col=myColors[8]),
  strip.border=list(col='black'))

bwtheme <- standard.theme("pdf", color = FALSE)
xyplot(value ~ variable1,
       globe.melt,
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
       par.settings = bwtheme)

xyplot(Trophic.Position ~ Sum.V | Region,
       t.globe, # Only include Western North Atlantic POM data
       group = Sample.ID2,
       cex = 1.5,
       # pch = 21,
       # col = "black",
       # bg = c(myColors),
       # pch=21,
       auto.key=list(columns=2, 
                     cex = 0.75,
                     space = "top",
                     points = TRUE),
       par.settings = my.settings,
       par.strip.text = list(col = 'white', font = 1),
       xlab = expression({Sigma}*"V"),
       ylab = "Trophic Position (Glu - Phe)")
# layout=c(1,3))


xyplot(Trophic.Position ~ Sum.V | Region,
       t.globe, # Only include Western North Atlantic POM data
       group = Sample.ID2,
       cex = 1.5,
       # pch = c(21:24)[as.numeric(globe$Region)],
       pch=21,
       # auto.key=list(columns=4, cex = 0.75),
       xlab = expression({Sigma}*"V"),
       ylab = "Trophic Position (Glu - Phe)")
       # layout=c(1,3))

xyplot(Trophic.Position ~ Sum.V | Type,
       t.globe, # Only include Western North Atlantic POM data
       group = Sample.ID2,
       cex = 1.5,
       # pch = c(21:24)[as.numeric(globe$Region)],
       pch=21,
       # auto.key=list(columns=4, cex = 0.75),
       xlab = expression({Sigma}*"V"),
       ylab = "Trophic Position (Glu - Phe)")
# layout=c(1,3))


xyplot(Dev.Val ~ Dev.Glx | Region,
       t.globe, # Only include Western North Atlantic POM data
       group = Sample.ID2,
       cex = 1.5,
       # pch = c(21:24)[as.numeric(globe$Region)],
       pch=21,
       # auto.key=list(columns=4, cex = 0.75),
       xlab = expression({Sigma}*"V"),
       ylab = "Valine Deviation (Val - Avg. Trophic)")

# layout=c(1,3))

