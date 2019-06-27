## ---
## Topic: Carbon CSIAA with only Stetson corals (no Larsen data or data from elsewhere)
## John Schiff
## Date: 02/15/2019
## ---

path1 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff c-csiaa stetson.csv'
path2 <- '~/Google Drive/projects/rproj/seus/data/schiff c-csiaa stetson.csv'

# "C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/stetson bulk read_only.csv" <-- stetson bulk data output

carbon <- read.csv(path1) # Check original excel sheet for standard deviations, etc. 
carbon <- read.csv(path2)

# Samples: 7, 22, 40, 109, 8, 25, 46, 103, 135, 170
yr <- stet.table$year.ad
carbon$Year.AD <- c(yr[7], yr[22], yr[40], yr[109], yr[8], yr[25], yr[46], yr[103], yr[135], yr[170])
EAA <- carbon %>% dplyr::select(Phe, Thr, Ile, Leu, Val) # Not including Lys because its peaks weren't great (large SD)
NEAA <- carbon %>% dplyr::select(Asp, Glu, Pro, Ala, Ser, Gly)
carbon$Avg.EAA <- rowMeans(EAA) 
carbon$Avg.NEAA <- rowMeans(NEAA)
carbon$RPP <- (0.9274*carbon$Avg.EAA) + 0.0769
carbon <- dplyr::arrange(carbon, Sample.No.)
# norm <- EAA - carbon$Avg.EAA
# norm$Avg.EAA <- rowMeans(norm)

xyplot(Phe ~ Year.AD,
       carbon,
       type = "o",
       pch=19,
       cex=1.5,
       ylab = expression(delta^{13}*"C Average EAA"),
       xlab = 'Calendar Year (C.E.)')   
