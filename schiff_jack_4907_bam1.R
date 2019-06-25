# ---
# Topic: Looking at bulk data from a bamboo coral, Jacksonville 4907 BAM1
# John Schiff
# Date: January 22, 2019
# ---

#########################################################
## This bamboo coral is from Jacksonville lithoherms   ##
## and also from dive 4907. I have gotten bulk isotope ##
## data from the coral so far, but am not sure on the  ##
## age. The coral disks look relatively young and not  ##
## diagenetically altered.                             ##
#########################################################

library(ggplot2)
library(lattice)
library(reshape2)
library(zoo)
library(dplyr)
library(tidyr)

# Import the data
bam.path1 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff 4907bam1 bulk data.csv'
bam.path2 <- '~/Google Drive/projects/rproj/seus/data/schiff 4907bam1 bulk data.csv'

bam1 <- read.csv(bam.path1)
bam1 <- read.csv(bam.path2)

bam1$avg.d13c <- rowMeans(bam1[,c(20:34)], na.rm =TRUE)

d15N <- bam1$avg.d15n
d13C <- bam1$avg.d13c

xyplot(d15N ~ sample,
       bam1,
       type = "o",
       pch = 20,
       cex = 1.5,
       xlab = "Sample ID",
       ylab = expression(delta^{15}*"N"))

xyplot(d13C ~ sample,
       bam1,
       type = "o",
       pch = 20,
       cex = 1.5,
       xlab = "Sample ID",
       ylab = expression(delta^{13}*"C"))

######################################
## Tidying the data because the     ##
## original file is a bit of a mess ##
######################################

library(tidyr)
library(reshape2)
n15 <- read.csv('C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff 4907bam1 d15n.csv')
c13 <- read.csv('C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff 4907bam1 d13c.csv')
n15 <- read.csv('~/Google Drive/projects/rproj/seus/data/schiff 4907bam1 d15n.csv')
c13 <- read.csv('~/Google Drive/projects/rproj/seus/data/schiff 4907bam1 d13c.csv')
n15 <- melt(n15, id.vars = "distance.mm")
c13 <- melt(c13, id.vars = "distance.mm")
table1 <- n15[c(63:nrow(n15)),]
table2 <- c13[c(63:nrow(n15)),]
bam <- cbind(table1,table2$value)
names(bam) <- c("Distance", "Treatment", "d15N", "d13C")
bam <- na.omit(bam)
bam$d13C <- as.character(bam$d13C) # For some reason, in the process these became characters or factors instead of numbers
bam$d13C <- as.numeric(bam$d13C)
bam$d15N <- as.numeric(bam$d15N)
bam$Distance <- as.numeric(bam$Distance)

#'
#'
#' Write the newly created dataframe, 
#' with treatmeants, to a .csv file
#' 
write.csv(bam, 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff_jack4907 bam1 by treatment.csv')
#' 

n_lattice <- lattice::bwplot(d15N ~ Distance, 
                data = bam,
                horiz = FALSE,
                xlab = "Distance (mm).",
                ylab = expression(paste(delta^{15},"N (\u2030)")))
c_lattice <- lattice::bwplot(d13C ~ Distance, 
                data = bam,
                horiz = FALSE,
                xlab = "Distance (mm)",
                ylab = expression(paste(delta^{13},"C (\u2030)")))
n_lattice
c_lattice

bw_theme <- trellis.par.get()
bw_theme$box.dot$pch <- "|"
bw_theme$box.rectangle$col <- "black"
bw_theme$box.rectangle$lwd <- 2
bw_theme$box.rectangle$fill <- "grey90"
bw_theme$box.umbrella$lty <- 1
bw_theme$box.umbrella$col <- "black"
bw_theme$plot.symbol$col <- "grey40"
bw_theme$plot.symbol$pch <- "*"
bw_theme$plot.symbol$cex <- 2
bw_theme$strip.background$col <- "grey80"

n_bw <- update(n_lattice, par.settings = bw_theme)
c_w <- update(c_lattice, par.settings = bw_theme)
par(mfrow=c(2,1))
n_bw
c_w
par(pty = "s")
grid.arrange(n_bw,c_w, nrow = 2, ncol=1)
# bam <- melt(bam1[,-c(3:7)], id.vars = "sample")
# var <- bam$variable
# bam$treatment <- NA 
# Code below is obsolete, but could be useful in other situations
# bam$treatment <- with(bam, ifelse(variable %in% "d15n.a", "a",
#                                   ifelse(variable %in% "d15n.b", "b",
#                                          ifelse(variable %in% "d15n.c", "c",
#                                                 ifelse(variable %in% "d15n.d", "d",
#                                                        ifelse(variable %in% "d15n.e", "e",
#                                                               ifelse(variable %in% "d15n.e", "e", "continue")))))))

lattice::bwplot(cyl.f ~ mpg, data = mtcars)
