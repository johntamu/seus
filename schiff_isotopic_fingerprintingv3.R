#' Created 08/08/2019
#'
#'
#'

library(dplyr)
library(MASS)
library(vegan)

newpath1 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff c-csiaa data from schiff_larsen_in prog.csv'
newpath2 <- '~/Documents/GitHub/data/schiff c-csiaa data from schiff_larsen_in prog.csv'
newpath3 <- '/home/john/Desktop/rproj/seus/data/schiff c-csiaa data from schiff_larsen_in prog.csv'

load.data <- newpath2

carbon <- read.csv(load.data)
head(carbon)

carbon <- carbon[order(carbon$Group.ID2),]

