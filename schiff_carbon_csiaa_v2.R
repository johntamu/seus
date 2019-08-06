#' 
#' John Schiff
#' 02-28-2019
#' 
#' Data analyses with d13C-CSIAA data
#' 

library(ggplot2)
library(lattice)
library(reshape2)
library(dplyr)
library(MASS)
library(vegan)

# Create pathnames
# Stetson CSIAA data by itself (although next CSV sheet includes Stetson data as well)
path1 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff c-csiaa stetson.csv'
path2 <- '~/Google Drive/projects/rproj/seus/data/schiff c-csiaa stetson.csv'
path3 <- '/home/john/Desktop/rproj/seus/data/schiff c-csiaa stetson.csv'

# Stetson data plus Larsen (2013) and other data
newpath1 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff c-csiaa data from schiff_larsen_in prog.csv'
newpath2 <- '~/Documents/GitHub/data/schiff c-csiaa data from schiff_larsen_in prog.csv'
newpath3 <- '/home/john/Desktop/rproj/seus/data/schiff c-csiaa data from schiff_larsen_in prog.csv'

load.data <- newpath2

# Load the data
carbon <- read.csv(load.data) # This dataset has Stetson and others in it, not just SEUS data

##############################
## SEUS C-CSIAA time series ## 
## and other analyses       ##
##############################

# Create a dataframe for Stetson only
df <- carbon %>% filter(Group.ID2 == "Leiopathes")
# Filter out standard deviations. Don't need this for now
seus_carbon <- df %>% 
  dplyr::select(Group.ID, Bulk, Ala, Asx, Glx, Gly, Ser, Pro, Ile, Leu, Lys, Phe, Thr, Tyr, Val, Year.AD)

# Determine Average EAA and average NEAA
EAA <- seus_carbon %>% dplyr::select(Phe, Thr, Ile, Leu, Val) # Not including Lys because its peaks weren't great (large SD)
NEAA <- seus_carbon %>% dplyr::select(Asx, Glx, Pro, Ala, Ser, Gly)

# Attach 
seus_carbon$EAA <- rowMeans(EAA) 
seus_carbon$NEAA <- rowMeans(NEAA)
# seus_carbon$RPP <- (0.9274*seus_carbon$Avg.EAA) + 0.0769
seus_carbon$RPP <- predict(vokh.lm, newdata = seus_carbon) # Using the vokh.lm model calculated below to predict RPP (instead of just doing equation of the line)

#' We now have a data frame with only SEUS C-CSIAA data
#' and we can use that to do time series analysis 
#' or other methods without utilizing the other C-CSIAA
#' datasets.
#' 
#' I will save this as a .csv file that can be used in the 
#' script file for plotting.
#' 
seus_carbon <- subset(seus_carbon, select = -c(Year.AD)) # We don't need the imported ages anymore
seus_carbon$Sample <- c(7,8,22,25,40,46,103,109,135,170)
df.stet %>%
  filter(sample.no. %in% seus_carbon$Sample) -> t.df
seus_carbon$Year.AD <- t.df$linear.ad2
write.csv(seus_carbon, "C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff c-csiaa stetson_cleaned.csv")

#'
#' Some statistical analysis
#' 

eaa.reg <- lm(EAA ~ Year.AD, seus_carbon[c(1:5),]) # Only on samples impacted by Seuss effect
neaa.reg <- lm(NEAA ~ Year.AD, seus_carbon[c(1:5),])
summary(eaa.reg) # Significant
summary(neaa.reg) # Not significant

t.test(seus_carbon$EAA, seus_carbon$NEAA) # Welch's t-test


##############################################
## Dataset used in Vokhshoori et al (2014)  ##
## which combines Larsen and Lehman data    ##
##############################################

vokh <- carbon %>% filter(Source == "Larsen" | Source == "Lehman")
vokh <- vokh %>% filter(Group.ID2 == "Eukaryote" | Group.ID2 == "Prokaryote")
vokh <- vokh %>% select(Group.ID, Group.ID2, Source, Bulk, Phe, Thr, Ile, Leu, Val, Lys, Asx, Glx, Pro, Ala, Ser, Gly)
vokh$EAA <- rowMeans(vokh %>%
                       select(Phe, Thr, Ile, Leu, Val))
# Recreate Vokhshoori (2014) linear regression curve
vokh.lm <- lm(Bulk ~ EAA, data = vokh)
summary(vokh.lm)
plot(vokh.lm)
plot(Bulk ~ EAA, data = vokh)

##############################################
## Dataset used in Schiff et al (2014)      ##
## which combines Larsen and Lehman data    ##
##############################################
#' Note: This data is slightly different from Vokhshoori (2014)
#' but it is still valid. Natasha and I just used different
#' phytoplankton samples but they had similar sample
#' names (e.g., "Santa Cruz dinoflagellate"). Nothing
#' to be alarmed over.
#' 


##########################
## EAA Normalized Data  ##
##########################

# Using Vokhshoori (2014) dataset, modified as in Schiff (2014) --- actually, only Lehman samples which includes Euk and Prok
vokh %>%
  dplyr::filter(Group.ID2 == "Eukaryote") %>%
  dplyr::select(-c(Group.ID, Group.ID2, Source, Bulk, EAA)) -> euk
vokh %>%
  dplyr::filter(Group.ID2 == "Prokaryote") %>%
  dplyr::select(-c(Group.ID, Group.ID2, Source, Bulk, EAA)) -> prok

vokh %>%
  dplyr::filter(Group.ID2 == "Eukaryote") %>%
  dplyr::select(EAA) -> euk_EAA # Isolate the EAA for this data
vokh %>%
  dplyr::filter(Group.ID2 == "Prokaryote") %>%
  dplyr::select(EAA) -> prok_EAA # Isolate the EAA for this data

vokh %>%
  dplyr::filter(Group.ID2 == "Eukaryote") %>%
  dplyr::select(Group.ID, Group.ID2) -> euk_IDs # Isolate sample IDs. There is probably a more "efficient" way to do this but I'm lazy
vokh %>%
  dplyr::filter(Group.ID2 == "Prokaryote") %>%
  dplyr::select(Group.ID, Group.ID2) -> prok_IDs # Isolate sample IDs

# Eukaryote normaliztion
euk.norm <- euk - euk_EAA$EAA # Since euk_EAA is technically
euk.means <- colMeans(euk.norm)
euk.sd <- apply(euk.norm, 2, sd)
euk.norm <- cbind(euk_IDs, euk.norm) # Put it back together

# Prokaryote normalization
prok.norm <- prok - prok_EAA$EAA
prok.means <- colMeans(prok.norm)
prok.sd <- apply(prok.norm, 2, sd)
prok.norm <- cbind(prok_IDs, prok.norm) # Put it back together

# Combine normalized Prokaryote and Eukaryote data
euk_prok <- rbind(euk.norm, prok.norm)

######################################
# EAA normalized Stetson C-CSIA data #
######################################

# seus_carbon %>%
#   dplyr::select(-c(Group.ID, Bulk, Year.AD, EAA, NEAA, RPP)) -> coral
seus_carbon %>% # REMEMBER: Lys was not included in average EAA calculation because of peak issue. However, I am plotting it.
  dplyr::select(Phe, Thr, Ile, Leu, Val, Lys, Asx, Glx, Pro, Ala, Ser, Gly) -> coral # re-ordering columns, multiple ways to do this
seus_carbon %>%
  dplyr::select(Group.ID) -> coral_IDs
coral.norm <- coral - seus_carbon$EAA

coral.means <- colMeans(coral.norm)
coral.sd <- apply(coral.norm, 2, sd)
coral.norm <- cbind(coral_IDs, coral.norm) # Put back together into data frame with IDs and normalized values
coral.norm %>%
  add_column(Group.ID2 = rep("Leiopathes", 10), .after = "Group.ID") -> coral.norm
  # so both normalized datasets have equal # columns

write.csv(coral.sd, "coral_sd.csv")
write.csv(euk.sd, "euk_sd.csv")
write.csv(prok.sd, "euk_sd.csv")

aa_norm <- rbind(euk_prok, coral.norm)

#' Now we have means for each: coral, euks and proks
#' We will attach the standard deviation values to them

means_sd_coral <- cbind(coral.means, coral.sd)
means_sd_euk <- cbind(euk.means, euk.sd)
means_sd_prok <- cbind(prok.means, prok.sd)

means_sd_all <- cbind(means_sd_coral, means_sd_euk, means_sd_prok)
names <- row.names(means_sd_all) # Kepp row names just in case you need them
means_sd_all <- as.data.frame(means_sd_all) 
data.table::setDT(means_sd_all, keep.rownames = TRUE) # Make row names into the first column
colnames(means_sd_all)[1] <- "AA" # Rename the first column if necessary
write.csv(means_sd_all, "means_sd_all_c-csia.csv")

#### Test plot below #####
test_plot <- ggplot(means_sd_all, aes(factor(AA), coral.means)) +
  geom_point(shape = 22, color = "black", fill = "#4575b4", size = 3) +
  geom_linerange(aes(ymin = coral.means -coral.sd, ymax = coral.means + coral.sd)) +
  geom_point(aes(AA, euk.means), shape = 22, color = "black", fill = "#d73027", size = 3) +
  geom_linerange(aes(ymin = euk.means - euk.sd, ymax = euk.means + euk.sd)) +
  geom_point(aes(AA, prok.means), shape = 23, color = "black", fill = "#fee090", size = 3) +
  geom_linerange(aes(ymin = prok.means - prok.sd, ymax = prok.means + prok.sd)) +
  scale_x_discrete(limits = names) +
  ylim(-10, 40) +
  ylab(c <- expression("Normalized "*delta^{13}*"C (\u2030)")) +
  xlab(NULL) +
  theme_bw()

test_plot

#' ---------------------------
#' Below is a regression of bulk d13C data from Stetson against 
#' d13C-AA data (EAA and NEAA)
#'
#' ---------------------------

# model <- lm(Bulk ~ Val, data = seus_carbon)
# summary(model)

########################
## C-CSIAA Mol % data ##
########################
mol <- read.csv('C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff_stetson aa mol percent.csv')
mol <- read.csv('~/Google Drive/projects/rproj/seus/data/schiff_stetson aa mol percent.csv')

mol$Total <- rowSums(mol[,2:13])
mol$Year.AD <- seus_carbon$Year.AD

mean(mol$Total)
sd(mol$Tota)

##################################
## C-CSIAA THAA Normalized data ##
##################################


##############################
## Suess effect correction  ##
##############################

#' Also relevant to bulk d13C data
#' McMahon (2015) used this method to correct
#' for the Suess effect in their data.
#' Two references: Francey et al, 1999; Quay et al, 2013
#' 
#' 0.16permil per decade since 1960
#' 0.05permil per decade between 1860 - 1960
#' 
#' One method: Make a for-loop to write a vector with corrected bulk d13C values
#' or just do it in Excel and re-import the corrected d13C values as a vector
#' 
#' For now, I'll reshape the data and just edit the output in Excel -- easier
#' 
library(reshape2)
View(seus_carbon)
seus_carbon %>% 
  melt(., "Year.AD") -> t.carbon
carbon_aa_to_suess_correct <- t.carbon[c(11:150),] # Get only Bulk and AAs
write.csv(carbon_aa_to_suess_correct, "carbon_aa_to_seus_correct.csv")

aa_suess_corrected <- read.csv("carbon_aa_seus_corrected.csv")
aa_suess_corrected %>%
  select(-c(value)) -> aa_suess_corrected
t.melt <- dcast(aa_seuss_corrected, Year.AD ~ variable)
##################
## UNUSED CODE  ##
##################