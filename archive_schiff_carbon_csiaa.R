## ---
## Topic: Carbon CSIAA data analysis from Gulf Stream black corals
## John Schiff
## Date: January 23, 2019
## ---

#############################################################
## This code will be used to analyze the Gulf Stream/SEUS  ##
## black coral C-CSIAA data. Any temporal work will likely ##
## be scripted between here and the seus_timeseries.R      ##
## script file.                                            ##
#############################################################

library(ggplot2)
library(lattice)
library(reshape2)
library(zoo)
library(dplyr)

path1 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff c-csiaa stetson.csv'
path2 <- '~/Google Drive/projects/rproj/seus/data/schiff c-csiaa stetson.csv'

carbon <- read.csv(path1) # Check original excel sheet for standard deviations, etc. 
carbon <- read.csv(path2)

carbon$Year.AD <- c(1984, 1936, 1877, 1587) # Ages are subject to change, be careful
carbon$Avg.EAA <- rowMeans(carbon[,c(5,7:9,14)]) # Lys not included; not great in chromatogram
carbon$Avg.NEAA <- rowMeans(carbon[, c(3,4,6,11,12,13)])

# Reconstructed Primary Production
# Using equation from Vokhshoori et al, 2014
carbon$RPP <- (0.9274*carbon$Avg.EAA) + 0.0769

##############
## Plotting ##
##############

xyplot(Avg.EAA ~ Year.AD,
       carbon,
       pch=19,
       cex=1.5,
       ylab = expression(delta^{13}*"C Average EAA"),
       xlab = 'Calendar Year (C.E.)')      

xyplot(Phe ~ Year.AD,
       carbon,
       pch=19,
       cex=1.5,
       ylab = expression(delta^{13}*"C Phe"),
       xlab = 'Calendar Year (C.E.)')

xyplot(Avg.NEAA ~ Year.AD,
       carbon,
       pch=19,
       cex=1.5,
       ylab = expression(delta^{13}*"C Average NEAA"),
       xlab = 'Calendar Year (C.E.)')

# Throw in the bulk d13C for the analyzed peels
carbon$bulk.d13C <- c(-16.99, -15.79, -15.31, -15.39)

##############
## Plotting ##
##############
xyplot(Avg.EAA ~ bulk.d13C,
       carbon,
       pch=19,
       cex=1.5,
       ylab = expression(delta^{13}*"C Average EAA"),
       xlab = expression(delta^{13}*"C Bulk"))

#############################################
## Code to create a barplot with Mol% data ##
#############################################

aa <- c('Coral.ID','Phe.Mol.','Gly.Mol.','Ser.Mol.','Lys.Mol.','Thr.Mol.','Glu.Mol.','Asp.Mol.','Ala.Mol.','Ile.Mol.','Leu.Mol.','Pro.Mol.','Val.Mol.')
carbon2 <- carbon %>% select(aa) 
carbon.means <- aggregate(carbon2[,2:13], list(Coral.ID=carbon2$Coral.ID), mean, na.rm = TRUE)
# Right now the column names are still 'AA.Mol.' which is kind of annoying, so let's fix that
# Note that there are multiple ways to do this, as usual
# colnames(carbon.means)[which(names(carbon.means) == "Phe.Mol.")] <- "Phe"
library(data.table)
setnames(carbon.means, 
         old = c('Glu.Mol.','Asp.Mol.','Ala.Mol.','Ile.Mol.','Leu.Mol.','Pro.Mol.','Val.Mol.','Gly.Mol.','Ser.Mol.','Lys.Mol.','Thr.Mol.','Phe.Mol.'),
         new = c('Glu', 
                 'Asp', 
                 'Ala', 
                 'Ile', 
                 'Leu', 
                 'Pro', 
                 'Val', 
                 'Gly', 
                 'Ser', 
                 'Lys', 
                 'Thr',
                 'Phe'))

carbon.melt <- melt(carbon.means, 'Coral.ID')

# Calculate standard deviation for every column using plyr
# plyr is an older package and there is likely a successor in the newer package 'dplyr'
library(plyr)
carbon.stdev <- ddply(carbon2, .(Coral.ID), colwise(sd))
setnames(carbon.stdev, 
         old = c('Glu.Mol.', 'Asp.Mol.', 'Ala.Mol.', 'Ile.Mol.', 'Leu.Mol.', 'Pro.Mol.', 'Val.Mol.', 'Gly.Mol.', 'Ser.Mol.', 'Lys.Mol.', 'Thr.Mol.', 'Phe.Mol.'),
         new = c('Glu', 'Asp', 'Ala',  'Ile', 'Leu', 'Pro', 'Val', 'Gly', 'Ser', 'Lys', 'Thr', 'Phe'))

carbon.melt$variable1 <- factor(carbon.melt$variable,
                                levels=c('Glu','Asp','Ala','Ile','Leu','Pro','Val','Gly','Ser','Lys','Thr','Phe')) # For Trophic/Source AA ordering

##############
## Plotting ##
##############
# Bar chart
barchart(value ~ variable1, # Make sure lattice is loaded
         carbon.melt,
         horiz = FALSE,
         ylab = 'Amino acids (mol, %)',
         col = "red",
         ylim=c(0,40))

carbon.melt2 <- melt(carbon2, "Coral.ID")
bwplot(value ~ variable,
         carbon.melt2,
         horiz = FALSE,
         ylab = 'Amino acids (mol, %)',
         col = "red")



####################################################
## Import the CSIAA data from Larsen et al (2013) ##
####################################################

newpath1 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff c-csiaa data from schiff_larsen_in prog.csv'
newpath2 <- '~/Google Drive/projects/rproj/seus/data/schiff c-csiaa data from schiff_larsen_in prog.csv'

training <- read.csv(newpath1)

#' Remember to test for normality across all data columns
#' before doing these analyses. PCA, LDA, et.c require
#' certain assumptions about the data.

# Example: lshap <-lapply(myData[,-c(1,2)], shapiro.test)

##################################################
## Based on Larsen et al (2013), this procedure ##
## requires the following packages: MASS,       ##
## cluster, vegan (for PCA)                     ##
##################################################

##################################################
## Based on McMahon et al (2015), we also need  ## 
## the package SIAR (or simmr) for the Bayesian ##
## stable isotope mixing model.                 ##
##################################################

# Helpful: https://rpubs.com/brouwern/veganpca

# The "training" data has standard deviation data, which is good, but we don't need right now
# Focus on Essential amino acids (for below, see Larsen et al (2013))
# EAA: Leu, Ile, Lys, Phe, Thr, Val, Met, His (but we don't have Hist and Met)
# NEAA: Asx, Glx, Ala, Gly, Tyr

training <- training[order(training$Group.ID2),] # Sort by name in Group 2 just to avoid any mixup later
training <- training %>% select(Group.ID, Group.ID2, Ala, Asx, Glx, Gly, Ile, Leu, Lys, Phe, Thr, Tyr, Val)
# training <- training[c(1:109),]

# In literature, PCA is done on normalized AA data
# eaa.avg <- rowMeans(training[,c(3:8)], na.rm = TRUE)
aa.avg <- rowMeans(training[,c(3:13)], na.rm = TRUE)

##############################################################
## NOTE: Did calculations, and AA-n in Larsen et al (2013)  ##
## is different. AA-n is every AA normalized to the overall ## 
## sample mean. NOT the average of the EAA only like in     ##             
## Schiff et al (2014).                                     ##
##############################################################

# training$Leu.norm <- training$Leu - eaa.avg
# training$Ile.norm <- training$Ile - eaa.avg

train.norm <- training[,3:13] - aa.avg
train.norm <- cbind(training[,1:2], train.norm)
train.norm <- train.norm %>% filter(Group.ID2 != "Collembola" & 
                                      Group.ID2 != "Daphnia" & 
                                      Group.ID2 != "Fish" & 
                                      Group.ID2 != "Seston" &
                                      Group.ID2 != "Soil" &
                                      Group.ID2 != "Mussels" &
                                      Group.ID2 != "Leiopathes")
train.norm <- na.omit(train.norm)

# larsen.pca2 <- princomp(train.norm2[-c(1:2)], cor = TRUE) # This recreates the PCA in Larsen (2013) Fig 1
train.dropped <- droplevels(train.norm) # Remove 'levels' not being used, since we subsetted the data
larsen.pca <- rda(train.dropped[-c(1:2)], scale = TRUE) # This recreates the PCA in Larsen as well, using vegan package

group.names <- levels(train.dropped$Group.ID2)

##############
## Plotting ##
##############

biplot(larsen.pca,
       display = c("species"),
       type = c("text"),
       col = c("#252525", "black"),
       xlab = "PC1 (32%)",
       ylab = "PC2 (22%)")
points(larsen.pca,
       pch = as.numeric(train.dropped$Group.ID2),
       col = "black",
       bg = as.numeric(train.dropped$Group.ID2))
legend("bottomleft",
       lty = NULL,
       legend = group.names,
       cex = 0.5,
       pch = unique(train.dropped$Group.ID2))

# points(train.pca2,
#        pch = c(21,22,23,24,25,9)[as.numeric(train.dropped$Sample.ID.2)],
#        col = "black",
#        bg = c('#d73027','#fc85d9','#fee090','#e0f3f8','#91bfdb','#4575b4')[as.numeric(train.dropped$Sample.ID.2)])

# ordihull(train.pca2,
#          draw = "polygon",
#          alpha = 25,
#          group = train.dropped$Sample.ID.2,
#          col = c('#d73027',
#                  '#fc8d59',
#                  '#fee090',
#                  '#e0f3f8',
#                  '#91bfdb',
#                  '#4575b4'))

# legend("bottomleft",
#        col = NULL,
#        fill = c('#d73027',
#          '#fc8d59',
#          '#fee090',
#          '#e0f3f8',
#          '#91bfdb',
#          '#4575b4'),
#        lty = NULL,
#        legend = group.names,
#        cex = 0.5,
#        pch = 21:24)

# 
# palette(c('#d73027',
#           '#fc85d9',
#           '#fee090',
#           '#e0f3f8',
#           '#91bfdb',
#           '#4575b4'))


################################################
## Combine Larsen (2013) training data        ##
## with our own sample data (deep-sea coral,  ##
## sediment traps, etc.)                      ##
################################################

# Helpful link: https://blog.bioturing.com/2018/06/18/how-to-read-pca-biplots-and-scree-plots/

trainData <- read.csv(newpath1)
trainData <- trainData[order(trainData$Group.ID2),] # Sort by name in Group 2 just to avoid any mixup later
trainData <- trainData %>% select(Group.ID, Group.ID2, Ala, Asx, Glx, Gly, Ile, Leu, Lys, Phe, Thr, Val) # My samples do not have Tyr
aa.avg <- rowMeans(trainData[,c(3:12)], na.rm = TRUE)


t.train <- trainData[,3:12] - aa.avg
t.train <- cbind(trainData[,1:2], t.train)
t.train <- t.train %>% filter(Group.ID2 != "Collembola" & 
                                Group.ID2 != "Daphnia" & 
                                Group.ID2 != "Fish" & 
                                Group.ID2 != "Seston" &
                                Group.ID2 != "Soil" &
                                Group.ID2 != "Mussels" &
                                Group.ID2 != "Leiopathes")
t.train <- na.omit(t.train)
# 
# t.train <- t.train %>% filter(Group.ID.2 == "Collembola" & 
#                                 Group.ID2 == "Daphnia" & 
#                                 Group.ID2 == "Fish" & 
#                                 Group.ID2 == "Seston" &
#                                 Group.ID2 == "Soil" &
#                                 Group.ID2 == "Mussels")

#' Finding: When including Isidella, the values are far away 
#' from the rest of them. This could be due to that the normalized
#' AA for this are normalized by the average AA overall and not
#' by the average EAA for each sample.
#' When done with Leiopathes, they are close to terrestrial plants (?).

t.train <- droplevels(t.train)
t.train <- t.train[order(t.train$Group.ID2),]
group.names <- levels(t.train$Group.ID2)

train.pca <- rda(t.train[-c(1:2)], scale = TRUE)
train.pca2 <- princomp(t.train[-c(1:2)], cor = TRUE)
# train.pca3 <- prcomp(train.dropped[-c(1:2)])


##############
## Plotting ##
##############

biplot(train.pca,
       display = c("species"),
       type = c("text"),
       col = c("#252525", "black"),
       # xlab = "PC1 (32%)",
       # ylab = "PC2 (22%)",
       main = "PCA")

# biplot(train.pca2,
#        display = c("sites", "species"),
#        type = c("text", "points"),
#        col = c("#252525", "black"),
#        main = "PCA")
points(train.pca,
       pch = as.numeric(t.train$Group.ID2),
       col = c('#66c2a5',
               '#fc8d62',
               '#8da0cb',
               '#e78ac3',
               '#a6d854',
               '#ffd92f',
               '#e5c494',
               '#b3b3b3')[as.numeric(t.train$Group.ID2)],
       bg = as.numeric(t.train$Group.ID2))
legend("topright",
       lty = NULL,
       legend = group.names,
       cex = 0.5,
       pch = unique(t.train$Group.ID2),
       col = c('#66c2a5',
                     '#fc8d62',
                     '#8da0cb',
                     '#e78ac3',
                     '#a6d854',
                     '#ffd92f',
                     '#e5c494',
                     '#b3b3b3')) # Consider unique() instead of as.numeric
# legend("bottomleft",
#        lty = 1,
#        legend = group.names,
#        cex = 0.5,
#        col = c('#66c2a5',
#                '#fc8d62',
#                '#8da0cb',
#                '#e78ac3',
#                '#a6d854',
#                '#ffd92f',
#                '#e5c494',
#                '#b3b3b3'))
ordihull(train.pca2,
         draw = "polygon",
         alpha = 50,
         group = train.dropped$Group.ID2,
         col = c('#66c2a5',
                 '#fc8d62',
                 '#8da0cb',
                 '#e78ac3',
                 '#a6d854',
                 '#ffd92f',
                 '#e5c494',
                 '#b3b3b3'))

# NOTE: Check to make sure the order of 'levels' after using the levels() function is the same as in the actual data sheet

#' ------------------------------------------ 
#' The code below will be used for 
#' LDA with the 'training' data set
#' and with my data. 
#' 
#' Per Larsen (2013), the LDA is
#' performed on the d13C-EAA.
#' 
#' The first thing we will do is re-create
#' the LDA figure (Fig. 2) from Larsen (2013)
#' ------------------------------------------

trainData <- read.csv(newpath1)
trainData <- trainData[order(trainData$Group.ID2),] # Sort by name in Group 2 just to avoid any mixup later
trainData <- trainData %>% select(Group.ID, Group.ID2, Ala, Asx, Glx, Gly, Ile, Leu, Lys, Phe, Thr, Tyr, Val)
# trainData <- trainData %>% select(Group.ID, Group.ID2, Ile, Leu, Lys, Phe, Thr, Val) # Isolate EAAs
trainData.larsen <- trainData %>% filter(Group.ID2 == "Bacteria" |
                                           Group.ID2 == "Fungi" |
                                           # Group.ID2 == "Isidella" |
                                           Group.ID2 == "Macroalgae" | 
                                           Group.ID2 == "Microalgae" | 
                                           Group.ID2 == "Plants" |
                                           Group.ID2 == "Seagrasses")
                                           # Group.ID2 == "Leiopathes")
trainData.larsen <- droplevels(trainData.larsen)
#' First thing I will do here though is re-do
#' the PCA from Larsen, but only using the EAA
#' rather than all of them.

trainData.larsen <- na.omit(trainData.larsen)

eaa.avg <- rowMeans(trainData.larsen[,c(3:8)], na.rm = TRUE)
larsen.norm <- trainData.larsen[,-c(1,2)] - eaa.avg
larsen.norm <- cbind(trainData.larsen[,c(1,2)], larsen.norm)

# lshap <-lapply(trainData.larsen[,-c(1,2)], shapiro.test)
# lres <- sapply(lshap, '[', c("statistic", "p.value"))

larsen.norm <- droplevels(larsen.norm)
group.names <- levels(larsen.norm$Group.ID2)

larsen.pca <- rda(larsen.norm[-c(1:2)], scale = TRUE)

##############
## Plotting ##
##############

biplot(larsen.pca,
       display = c("species"),
       type = c("text"),
       col = c("#252525", "black"),
       main = "PCA")

points(larsen.pca,
       pch = as.numeric(larsen.norm$Group.ID2),
       col = c('#66c2a5',
               '#fc8d62',
               '#8da0cb',
               '#e78ac3',
               '#a6d854',
               '#ffd92f',
               '#e5c494',
               '#b3b3b3')[as.numeric(larsen.norm$Group.ID2)],
       bg = as.numeric(larsen.norm$Group.ID2))
legend("topright",
       lty = NULL,
       legend = group.names,
       cex = 0.5,
       pch = unique(larsen.norm$Group.ID2),
       col = c('#66c2a5',
               '#fc8d62',
               '#8da0cb',
               '#e78ac3',
               '#a6d854',
               '#ffd92f',
               '#e5c494',
               '#b3b3b3'))

#' Now for the actual LDA recreation

larsen.lda <- lda(Group.ID2 ~ Ile + Leu + Lys + Phe + Thr + Val,
                  data=trainData.larsen,
                  CV = FALSE)

newData <- trainData %>% filter(Group.ID2 == "Isidella")
newData <- droplevels(newData)
predicted <- predict(larsen.lda, newData)

#' Now that we have successfully recreated the PCA and LDA with the Larsen (2013)
#' dataset, we can use it on new data. Here is where you think about the
#' following things: What are the carbon sources to my corals? We can
#' make a custome dataset with d13C-AA data from the sources we are
#' interested in. 
#' 
#' In my case, we think about the sources of carbon to black corals and
#' bamboo corals: Bacteria, microalgae, macroalgae. Likely no terrestrial 
#' plants. From there, what about N2-fixers and non-N2 fixers?
#' 
#' d13C of POM from GOM (Demopolous et al 2010) is around -20 permil.
#' 
#' 
#' What I will do first is create my own 'modified' training set for PCA and LDA
#' based on what are the mostly likely contributors to carbon source for black
#' corals in SEUS and Gulf of Mexico. 
