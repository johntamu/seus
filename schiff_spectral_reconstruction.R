#' 
#' John Schiff
#' 3/1/2019
#' 
#' 
#' This is a script that writes .csv files of isotope 
#' vectors to use in the kSpectra (aka SSA-MTM toolkit)
#' and do spectral analysis without R code.
#' 
#' I do this mainly to check what I get there vs. what I get
#' in the Rssa package. 
#' 
#' This assumes that the data has already been loaded
#' from another script file.
#' 
#' NOTE: 3/10/2019
#' The Rssa package references the techniques and math used
#' by the creators of the SSA-MTM Toolkit
#' 
#' For SSA, I recommend using the Rssa package (or kSpectra, Brendan has the full license)
#' because it is a bit of a pain using the Toolkit
#' 
#' However, the toolkit is still useful for using the Multitaper Method (MTM).
#' 
#' 
# Stetson-4904
write.table(df.stet$d15n, 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/toolkit_vectors/stetson_d15n.txt', 
            row.names = FALSE, col.names = FALSE, sep = "\t", na = "NaN")
write.table(df.stet$d13c, 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/toolkit_vectors/stetson_d13c.txt', 
            row.names = FALSE, col.names = FALSE, sep = "\t", na = "NaN")
# Jacksonville-4907
write.table(df.jack$d15n, 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/toolkit_vectors/jack_d15n.txt', 
            row.names = FALSE, col.names = FALSE, sep = "\t", na = "NaN")
write.table(df.jack$d13c, 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/toolkit_vectors/jack_d13c.txt', 
            row.names = FALSE, col.names = FALSE, sep = "\t", na = "NaN")
# Savannah-4907
write.table(df.sav$d15n, '~/Google Drive/projects/rproj/seus/data/toolkit_vectors/sav_d15n.txt', 
            row.names = FALSE, col.names = FALSE, sep = "\t", na = "NaN")
write.table(df.sav$d13c, '~/Google Drive/projects/rproj/seus/data/toolkit_vectors/sav_d13c.txt', 
            row.names = FALSE, col.names = FALSE, sep = "\t", na = "NaN")

# Stetson Iodine
# write.table()


#' After analyzing in SSA-MTM toolkit,
#' we import the vectors so we can plot
#' them in here and make them look good.

# Stetson-4904
# Window Length = 40
stet_ssa_n <- read.delim('~/Google Drive/projects/rproj/seus/data/toolkit_output/stet_n.txt', header = FALSE) # SSA on bulk d13C
stet_ssa_c <- read.delim('~/Google Drive/projects/rproj/seus/data/toolkit_output/stet_c.txt', header = FALSE) # SSA on bulk d15N
df.stet$ssa.recon.d15n <- stet_ssa_n$V1 # Since it reads in the text file as a data frame
df.stet$ssa.recon.d13c <- stet_ssa_c$V1
  
# Jacksonville-4907
jack_ssa_n <- read.delim('~/Google Drive/projects/rproj/seus/data/toolkit_output/jack_n.txt', header = FALSE) # SSA on bulk d13C
jack_ssa_c <- read.delim('~/Google Drive/projects/rproj/seus/data/toolkit_output/jack_c.txt', header = FALSE) # SSA on bulk d15N
df.jack$ssa.recon.d15n <- jack_ssa_n$V1
df.jack$ssa.recon.d13c <- jack_ssa_c
# Savannah-4902
sav_ssa_n <- read.delim('~/Google Drive/projects/rproj/seus/data/toolkit_output/sav_n.txt', header = FALSE) # SSA on bulk d13C
sav_ssa_c <- read.delim('~/Google Drive/projects/rproj/seus/data/toolkit_output/sav_c.txt', header = FALSE) # SSA on bulk d15N
df.sav$ssa.recon.d15n <- sav_ssa_n$V1
df.sav$ssa.recon.d13c <- sav_ssa_c

############
# Plotting # Copy-paste this to figures script file
############

plot(d15n ~ linear.ad, data = df.stet,
     type = "l",
     col = alpha("black", 0.4))
lines(ssa.recon.d15n ~ linear.ad, df.stet, col = "red", lwd = 2)

plot(ssa.recon.d15n ~ linear.ad, df.jack, col = 'red', type = "l")
lines(ssa.recon.d15n ~ linear.ad, df.sav, col = 'blue')


############################ Imported from another script 3/10/2019
##                        ##
##                        ##
## SSA with Rssa package  ##
##                        ##
##                        ##
############################
library(Rssa)

head(df.stet)
head(df.jack)
head(df.sav)
head(df.jack4686)



