#' -------------------------------------------------------------
#' Topic: Continuing C-CSIAA work from the part 1 script file 
#' John Schiff
#' Date: January 31, 2019
#' -------------------------------------------------------------
#' 
#' Now that we have successfully recreated the PCA and LDA with the Larsen (2013)
#' dataset, we can use it on new data. Here is where you think about the
#' following things: What are the carbon sources to my corals? We can
#' make a custom dataset with d13C-AA data from the sources we are
#' interested in. 
#' 
#' In my case, we think about the sources of carbon to black corals and
#' bamboo corals: Bacteria, microalgae, macroalgae. Likely no terrestrial 
#' plants. From there, what about N2-fixers and non-N2 fixers?
#' 
#' d13C of POM from GOM (Demopolous et al 2010) is around -20 permil.
#' 
#' What I will do first is create my own 'modified' training set for PCA and LDA
#' based on what are the mostly likely contributors to carbon source for black
#' corals in SEUS and Gulf of Mexico. 
#' 
#' This script will be redoing some of the code from Part 1, but
#' the code should be a bit more refined.

library(dplyr)
library(MASS)
library(vegan)

# Helpful link: https://rpubs.com/Nolan/298913

#' Load the data
newpath1 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff c-csiaa data from schiff_larsen_in prog.csv'
newpath2 <- '~/Google Drive/projects/rproj/seus/data/schiff c-csiaa data from schiff_larsen_in prog.csv'
trainData <- read.csv(newpath1)
trainData <- read.csv(newpath2)
trainData <- trainData[1:141,]
trainData <- trainData[order(trainData$Group.ID2),] # Sort by name in Group 2 to avoid any mixups
larsenData <- trainData %>% dplyr::filter(Group.ID2 == "Bacteria" |
                                     Group.ID2 == "Fungi" |
                                     Group.ID2 == "Macroalgae" | 
                                     Group.ID2 == "Microalgae" | 
                                     Group.ID2 == "Plants" |
                                     Group.ID2 == "Seagrasses")
seusDataC <- trainData %>% dplyr::filter(Group.ID2 == "Leiopathes")
seusDataC$Peel <- c(7, 22, 40, 109)
seusDataC$Coral <- c(rep("Stetson Banks-4904"))

seusDataC <- seusDataC %>% dplyr::select(Coral, Peel, Ala, Asx, Glx, Gly, Ile, Leu, Lys, Phe, Pro, Ser, Thr, Tyr, Val)

################################################
## Principal Components Analysis with Avg AA  ##
################################################

AA <- trainData %>% dplyr::select(Group.ID2, Phe, Thr, Ile, Leu, Val, Lys, Asx, Glx, Ala, Gly, Tyr)
AA <- na.omit(AA)
avg.aa <- rowMeans(AA[,c(2:ncol(AA))])


# Isolate EAAs that we care about in this case
larsenData <- larsenData %>% dplyr::select(Group.ID, Group.ID2, Ile, Leu, Phe, Thr, Val)
larsenData <- larsenData %>% dplyr::filter(Group.ID2 != "Fungi" & Group.ID2 != "Plants")
larsenData <- droplevels(larsenData)
larsenData <- na.omit(larsenData) # Remove rows with NAs



larsen.lda <- lda(Group.ID2 ~ Ile + Leu + Phe + Thr + Val,
                  data = larsenData)

  
  
  
  



