######################################
##                                  ##
## d13C-CSIAA Statistical analysis  ##
##                                  ##
######################################

#' Note: I do not recreate the figure from 
#' Larsen et al (2013), but I have the code to
#' do it in the script file: archive_schiff_carbon_csiaa.R
#' 
#' 
library(dplyr)
library(MASS)
library(vegan)

# Helpful link: https://rpubs.com/Nolan/298913

# Make sure data is loaded
head(carbon) # Loaded from carbon_csiaa script file

# Filter the data and normalize
carbon <- carbon[order(carbon$Group.ID2),] # Sort by name in Group 2 just to avoid any mixup later
carbon %>% # We are doing PCA on only the essential amino acids; McMahon et al (2014) uses four EAA
  dplyr::select(Group.ID2, Phe, Thr, Ile, Leu, Val) %>% # However, McMahon 2018 includes Valine (so 5 EAA)
  filter(Group.ID2 == "Prokaryote" 
         | Group.ID2 == "Isidella"
         | Group.ID2 == "Leiopathes"
         | Group.ID2 == "Sediment trap"
         | Group.ID2 == "Eukaryote") -> df_eaa # Now we need to normalize to the mean EAA for each sample
eaa_norm <- df_eaa[,-c(1)] - rowMeans(df_eaa[,-c(1)])
eaa_norm <- cbind(df_eaa[,1], eaa_norm)

carbon %>%
  filter(Group.ID2 == "Leiopathes") %>%
  dplyr::select(Group.ID3, Phe, Thr, Ile, Leu, Val) -> df_coral 

coral_norm <- df_coral[,-c(1)] - rowMeans(df_coral[,-c(1)])
coral_norm <- cbind(df_coral[,1], coral_norm)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#'
#'
#' Below are other kinds of filtered datasets.
#' I have decided to collate them here at the top
#' rather than throughout the file.
#' 

# Meant to mimic the idea in McMahon (2015)
carbon %>%
  dplyr::select(Group.ID, Group.ID2, Group.ID3, Source, Phe, Thr, Ile, Leu, Val) %>%
  filter(Group.ID3 == "N2 fixing" |
           Group.ID3 == "Non-N2 fixing" |
           Group.ID3 == "Eukaryotic algae" |
           Group.ID3 == "Heterotrophic bacteria" |
           Group.ID3 == "Leiopathes-post" |
           Group.ID3 == "Leiopathes-pre") -> pca_set

# Here, I filter Macroalgae, Microalgae, Bacteria, Seagrasses, Terrestrial plants
carbon %>%
  filter(Source == "Larsen 2013") %>%
  filter(Group.ID2 == "Bacteria"
         | Group.ID2 == "Macroalgae"
         | Group.ID2 == "Microalgae"
         | Group.ID2 == "Plants"
         | Group.ID2 == "Seagrasses") -> t.larsen

carbon %>%
  dplyr::select(Group.ID3, Phe, Thr, Ile, Leu, Val) %>%
  filter(Group.ID3 == "N2 fixing" |
           Group.ID3 == "Non-N2 fixing" |
           Group.ID3 == "Eukaryotic algae" |
           Group.ID3 == "Heterotrophic bacteria") -> training_set

training_norm <- training_set[,-c(1)] - rowMeans(training_set[,-c(1)])
training_set2 <- cbind(training_set[,1], training_norm)

# Group.ID3 == "Leiopathes-post" |
# Group.ID3 == "Leiopathes-pre") -> pca_set


####################################
##                                ##
##                                ##
## Principal Components Analysis  ## 03/07/2019
## Part III.                      ##
##                                ##
####################################
#' Here I am using PCA on the training sets
#' I create below for LDA. I may as well see how
#' they look.
#' 


#' This is for using biplot in vegan package.
#' I have mostly abandoned this for the flexibility 
#' with prcomp() function and ggplot2 for plotting
#' 


pca_set <- droplevels(pca_set)
pca_set2 <- pca_set[,c(5:9)] - rowMeans(pca_set[,c(5:9)])
t.pca <- rda(pca_set2, scale = TRUE)
df_pca <- prcomp(pca_set2, center = TRUE, scale. = TRUE)
summary(eigenvals(t.pca)) # Important for getting variance in vegan objects

par(pty = "s")
biplot(t.pca,
       display = c("species"),
       type = c("text"),
       col = c("#252525", "black"),
       xlab = "PC1 (58%)",
       ylab = "PC2 (23%)")
points(t.pca,
       pch = c(21,22,23,23,25)[as.numeric(pca_set$Group.ID3)],
       cex = 1.25,
       col = "black",
       bg = c("#2b83ba", "#abdda4", "#ffffbf", "#fdae61", "#d7191c")[as.numeric(pca_set$Group.ID3)])
legend("bottomleft",
       lty = NULL,
       bty = "n",
       pch = c(21,22,23,23,25),
       legend = levels(pca_set$Group.ID3),
       cex = 0.65,
       pt.cex = 0.95,
       pt.bg = c("#2b83ba", "#abdda4", "#ffffbf", "#fdae61", "#d7191c"))

####################################
##                                ##
##                                ##
## Principal Components Analysis  ## 03/08/2019; UPDATED 6/25/2019
## Part IV.                       ##
##                                ##
####################################
#' This is probably just going
#' to re-do what I have already done
#' but I don't want to erase old code (could still be useful).
#' Here, I filter Macroalgae, Microalgae, Bacteria, Seagrasses, Terrestrial plants

#' 6/25/2019 trying a new way to do PCA and then plotting it with ggplot
#' 
#' 
#' 
#' 
#' 
#' 

# pca_set <- droplevels(pca_set)
# pca_set2 <- pca_set[,c(5:9)] - rowMeans(pca_set[,c(5:9)])
# t.pca <- rda(pca_set2, scale = TRUE)
# df_pca <- prcomp(pca_set2, scale = TRUE)

# df_pca <- prcomp(t.norm, scale = TRUE)
df_pca <- prcomp(training_set2[,-c(1)], scale = TRUE)
df_out <- as.data.frame(df_pca$x)
df_out$group <- sapply(strsplit(as.character(row.names(training_set2)), "_"), "[[", 1)
head(df_out)

df_pca$rotation
df_pca$x
summary(df_pca)

# prediction of PCs for validation dataset
pred <- predict(df_pca, newdata=coral_norm[,2:6])
pred <- as.data.frame(pred)
legend_title <- NULL

p <- ggplot(df_out, aes(x=PC1,y=PC2, fill = training_set2[,c(1)], shape = training_set2[,c(1)]))
p + geom_point(size = 3.5, color = "black") +
  theme_classic() +
  # scale_color_manual(legend_title, values = cbPalette) +
  scale_shape_manual(legend_title, values = c(21:24)) +
  scale_fill_manual(legend_title, values = cbPalette) +
  # xlim(-2,2) +
  # ylim(-2,2) +
  xlab("PC1 (45%)") +
  ylab("PC2 (34%)") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  # geom_point(data = pred, aes(x=PC1, y=PC2), shape = 4, fill = "black", size = 4) +
  theme(axis.text.y   = element_text(size=12, color = "black"),
        axis.text.x   = element_text(size=12, color = "black"),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12))


# pca1 <- princomp(t.norm, cor = TRUE, scores = TRUE)
col8 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00') # 8 variables needing colors
col4 <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#2166ac", "#b2182b")


##################################
##                              ##
##                              ##
## Linear Discriminant Analysis ##
##                              ##
##                              ##
##################################
# Helpful link: https://rstudio-pubs-static.s3.amazonaws.com/35817_2552e05f1d4e4db8ba87b334101a43da.html
# Helpful link: https://datascienceplus.com/how-to-perform-logistic-regression-lda-qda-in-r/
# Use lda() function in MASS package
#' The few papers that do this all have the same
#' general idea, but the data varies a little bit.
#' For example, McMahon (2015) uses a modified dataset
#' combining Larsen (2013) and McCarthy lab culture data.
#' Vokhshoori (2014) uses Larsen data as they did in Larsen (2013)
#' 
#' For now, I will use Larsen data in a similar way that Vokhshoori (2014)
#' and Schiff (2014) used it. Kelton did the LDA for my Schiff (2014) paper.
#' 
t.schiff <- droplevels(t.schiff)
t.schiff %>%
  dplyr::select(Group.ID2, Phe, Thr, Ile, Leu, Val) -> t.schiff

carbon %>%
  filter(Source == "Larsen 2013") %>%
  filter(Group.ID2 == "Bacteria"
         | Group.ID2 == "Macroalgae"
         | Group.ID2 == "Microalgae"
         | Group.ID2 == "Plants"
         | Group.ID2 == "Seagrasses") -> t.larsen
t.larsen <- droplevels(t.larsen)
t.larsen %>%
  dplyr::select(Group.ID2, Phe, Thr, Ile, Leu, Val) -> larsen_eaa # dplyr may clash with MASS package
larsen_eaa_norm <
lda_larsen <- lda(Group.ID2 ~ Phe + Thr + Ile + Leu + Val, larsen_eaa, CV = FALSE)
lda_larsen_values <- predict(lda_larsen)
lda_predicted_schiff <- predict(lda_larsen, newdata = t.schiff)

# plot(lda_larsen, dimen = 2, abbrev = 1, ylim = c(-6,6), note 3/8/2019 -- keep working on this
#      panel = function(x, y, ...)
#      points(x, y, ...),
#      col = as.numeric(larsen_eaa$Group.ID2), pch = c(21,22,24))
par(pty = "s")
plot(lda_larsen, dimen = 2, abbrev = FALSE, ylim = c(-6,6), xlim = c(-6,6),
     pch = c(21,22,24)[as.numeric(larsen_eaa$Group.ID2)])
points(lda_predicted_schiff$x, pch = 23, bg = "#41b6c4", col = "black", cex = 1.25)
legend("topright",
       inset = c(0, 0),
       lty = NULL,
       bty = "n",
       pch = c(23,23,23,25,25),
       legend = levels(t.schiff$Group.ID2),
       cex = 0.65,
       pt.cex = 0.95,
       pt.bg = c("#41b6c4","#7fbf7b", "#ffffbf", "#fee090", "#fc8d59", "#d73027"))

plot(lda_predicted_schiff$x[,1],lda_predicted_schiff$x[,2], xlim = c(-6,6), ylim = c(-6,6)) # make a scatterplot
points(lda_larsen_values$x[,1],lda_larsen_values$x[,2])

table(lda_predicted_schiff$class, t.schiff$Group.ID2)



######################################
## We are now using the LDA method  ##
## on the Group.ID3 grouping 
## (based on McMahon 2015)
######################################

# carbon %>%
  # filter(Source == "Larsen 2013") %>%
  # filter(Group.ID3 == "N2 fixing" |
  #        Group.ID3 == "Non-N2 fixing" |
  #        Group.ID3 == "Euk microalgae" ) -> set1
  #        # Group.ID3 == "Euk macroalgae" |
  #        # Group.ID3 == "Het bacteria" |
  #        # Group.ID3 == "Pseudo" |
  #          # Group.ID3 == "Dinoflag") -> set1
  
carbon %>% 
  filter(Group.ID3 != "Isidella" & Group.ID3 != "Sediment trap" & Group.ID3 != "") -> set1
carbon %>% 
  filter(Group.ID3 != "" & Group.ID3 != "Leiopathes-pre" 
         & Group.ID3 != "Leiopathes-post" & Group.ID3 != "Isidella" & Group.ID3 != "Sediment trap") -> set1


carbon %>%
  filter(Source == "Schiff") %>%
  filter(Group.ID3 == "Leiopathes-post") -> set2
carbon %>%
  filter(Source == "Schiff") %>%
  filter(Group.ID3 == "Leiopathes-pre") -> set3

set2 <- droplevels(set2)
set2 %>%
  dplyr::select(Group.ID3, Phe, Thr, Ile, Leu) -> set2 # Can add Val back in
set2_norm <- set2[,c(2:5)] - rowMeans(set2[,c(2:5)])
set2_norm <- cbind(set2$Group.ID3, set2_norm)
colnames(set2_norm)[colnames(set2_norm)=="set2$Group.ID3"] <- "Group.ID3"

set3 <- droplevels(set3)
set3 %>%
  dplyr::select(Group.ID3, Phe, Thr, Ile, Leu) -> set3 # Can add Val back in
set3_norm <- set3[,c(2:5)] - rowMeans(set3[,c(2:5)])
set3_norm <- cbind(set3$Group.ID3, set3_norm)
colnames(set3_norm)[colnames(set3_norm)=="set3$Group.ID3"] <- "Group.ID3"


set1 <- droplevels(set1)
set1 %>%
  dplyr::select(Group.ID3, Phe, Thr, Ile, Leu) -> larsen_eaa # dplyr may clash with MASS package
larsen_eaa_norm <- larsen_eaa[,c(2:5)] - rowMeans(larsen_eaa[,c(2:5)])
larsen_eaa_norm <- cbind(larsen_eaa$Group.ID3, larsen_eaa_norm)
colnames(larsen_eaa_norm)[colnames(larsen_eaa_norm)=="larsen_eaa$Group.ID3"] <- "Group.ID3"

lda_larsen <- lda(Group.ID3 ~ Phe + Thr + Ile + Leu, larsen_eaa, CV = TRUE) # Doing this on normalized values gets an error: variables are collinear
lda_larsen_values <- predict(lda_larsen, dimen=2)$x
lda_predicted_set2 <- predict(lda_larsen, newdata = set2)
lda_predicted_set3 <- predict(lda_larsen, newdata = set3)

par(pty = "s")
plot(lda_larsen, dimen = 2, abbrev = TRUE, ylim = c(-6,6), xlim = c(-6,6),
     pch = c(21,22,24)[as.numeric(larsen_eaa$Group.ID3)])

eqscplot(lda_larsen_values, type = "p", pch = 22, cex = 1.25, bg = col4[as.numeric(larsen_eaa$Group.ID3)], xlab = "LD1", ylab = "LD2", tol = 0.25, las =1)
# legend("topright",
#        inset = c(0, 0),
#        lty = NULL,
#        bty = "n",
#        cex = 0.55,
#        # pch = c(22),
#        # col = "black",
#        fill = col4,
#        legend = levels(larsen_eaa$Group.ID3))

points(lda_predicted_set2$x, pch = 23, bg = "#2166ac", col = "black", cex = 1.25) # Set 2, post-1900
points(lda_predicted_set3$x, pch = 23, bg = "#b2182b",col = "black", cex = 1.25) # Set 3, pre-1900
legend("topleft",
       # inset = c(0, 0),
       lty = NULL,
       bty = "n",
       cex = 0.75,
       col = "black",
       pch = c(22,22,22,22,23,23),
       pt.bg = col4,
       pt.cex = 1,
       legend = c("Euk. algae", "Het. bacteria", "N2 fixing", "Non-N2 fixing", "post-1900", "pre-1900"))
       # legend = levels(larsen_eaa$Group.ID3))
       # cex = 0.65,
       # pt.cex = 0.95,
       # pt.bg = c("#41b6c4","#7fbf7b", "#ffffbf", "#fee090", "#fc8d59", "#d73027"))

plot(lda_predicted_set2$x[,1],lda_predicted_set2$x[,2], xlim = c(-6,6), ylim = c(-6,6)) # make a scatterplot
points(lda_larsen_values$x[,1],lda_larsen_values$x[,2])

table(lda_predicted_schiff$class, t.schiff$Group.ID2)




################################################
## Unrelated code for isotopic fingerprinting ##
################################################
# Helpful link: http://huboqiang.cn/2016/03/03/RscatterPlotPCA
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", '#003c30', '#8c510a')
test.set <- seus_carbon[,c(9,10,12,13)]
test.set <- test.set - rowMeans(test.set)
# test.pca <- rda(test.set, scale = TRUE)
test.pca <- prcomp(test.set, scale = TRUE)


# Using ggplot
plot(test.pca$x[,1], test.pca$x[,2])
df_out <- as.data.frame(test.pca$x)
df_out$group <- sapply(strsplit(as.character(row.names(df)), "_"), "[[", 1)
head(df_out)

library(ggplot2)
library(grid)
library(gridExtra)

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=seus_carbon$Group.ID))
p<-p+geom_point()
p

biplot(test.pca,
       display = c("species"),
       type = c("text"),
       col = c("#252525", "black"),
       xlab = "PC1",
       ylab = "PC2",
       xlim = c(-2.5,2.5),
       ylim = c(-2.5,2.5))
points(test.pca,
       pch = 22,
       cex = 1.25,
       # col = "black",
       bg = cbPalette[as.numeric(seus_carbon$Group.ID)])
legend("bottomleft",
       lty = NULL,
       bty = "n",
       pch = c(21,22,23,23,25),
       legend = levels(pca_set$Group.ID3),
       cex = 0.65,
       pt.cex = 0.95,
       pt.bg = c("#2b83ba", "#abdda4", "#ffffbf", "#fdae61", "#d7191c"))

#'
#'
#' UNUSED CODE
#' 
#' 

# carbon %>% 
#   filter(Group.ID5 != "" & Group.ID5 != "Leiopathes-pre" & Group.ID5 != "Leiopathes-post") -> t.larsen
# 
# carbon %>%
#   filter(Source == "Schiff") %>%
#   filter(Group.ID5 == "Leiopathes-pre" | Group.ID5 == "Leiopathes-post") -> t.schiff
# 
# t.compile <- rbind(t.larsen, t.schiff)
# write.csv(t.compile, "temp_csiaa.csv")
# t.compile <- read.csv("temp_csiaa.csv")
# t.compile %>%
#   dplyr::select(Group.ID5, Phe, Thr, Ile, Leu) -> t.compile
# t.compile <- na.omit(t.compile)
# t.compile <- droplevels(t.compile)
# t.norm <- t.compile[,-c(1)] - rowMeans(t.compile[,-c(1)])
# t.norm




####################################
##                                ##
##                                ##
## Principal Components Analysis  ##
##                                ##
##                                ##
####################################

# eaa_norm <- droplevels(eaa_norm) # Remove 'levels' not being used, since we subsetted the data
# eaa_pca <- rda(eaa_norm[,-c(1)], scale = TRUE)
# # eaa_pca <- prcomp(eaa_norm[,-c(1)], center = TRUE, scale. = TRUE)
# # eaa_pca <- princomp(eaa_norm[,-c(1)], cor = TRUE, score = TRUE)
# 
# coral_norm <- droplevels(coral_norm)
# coral_pca <- rda(coral_norm[,-c(1)], scale = TRUE)
# 
# # predicted <- predict(eaa_pca, newdata = coral_norm[,-c(1)]) # Not sure this works?
# 
# ##############
# ## Plotting ##
# ##############
# # If using rda() from vegan
# biplot(eaa_pca,
#        display = c("species"),
#        type = c("text"),
#        col = c("#252525", "black"),
#        xlab = "PC1",
#        ylab = "PC2")
# points(eaa_pca,
#        pch = c(21,22,23,24,25)[as.numeric(eaa_norm$`df_eaa[, 1]`)],
#        # col = c("red", "yellow", "blue", "green", "orange")[as.numeric(eaa_norm$`df_eaa[, 1]`)],
#        col = "black",
#        bg = c("red", "yellow", "blue", "darkgreen", "orange")[as.numeric(eaa_norm$`df_eaa[, 1]`)])
# legend("bottomleft",
#        lty = NULL,
#        bty = "n",
#        pch = c(21,22,23,24,25),
#        legend = levels(eaa_norm$`df_eaa[, 1]`),
#        cex = 0.45,
#        pt.cex = 0.95,
#        # col = "black",
#        pt.bg = c("red", "yellow", "blue", "darkgreen", "orange"))
# 
# biplot(coral_pca,
#        display = c("species"),
#        type = c("text"),
#        col = c("#252525", "black"),
#        xlab = "PC1",
#        ylab = "PC2")
# points(coral_pca,
#        pch = c(21,22,23,24,25)[as.numeric(coral_norm$`df_coral[, 1]`)],
#        # col = c("red", "yellow", "blue", "green", "orange")[as.numeric(eaa_norm$`df_eaa[, 1]`)],
#        col = "black",
#        bg = c("red", "yellow", "blue", "darkgreen", "orange")[as.numeric(coral_norm$`df_eaa[, 1]`)])
# legend("bottomleft",
#        lty = NULL,
#        bty = "n",
#        pch = c(21,22,23,24,25),
#        legend = levels(coral_norm$`df_coral[, 1]`),
#        cex = 0.45,
#        pt.cex = 0.95,
#        # col = "black",
#        pt.bg = c("red", "yellow", "blue", "darkgreen", "orange"))
# 
# # If using princomp() from base stats
# biplot(eaa_pca,
#        # col = c("#252525", "black"),
#        xlab = "PC1",
#        ylab = "PC2",
#        col=c(2,3), cex=c(1/2, 1.25))
# legend("bottomleft",
#        lty = NULL,
#        bty = "n",
#        legend = levels(eaa_norm$`df_eaa[, 1]`),
#        cex = 0.45,
#        pt.cex = 0.95)
#        # col = "black",
#        # pt.bg = c("red", "yellow", "blue", "darkgreen", "orange"))

####################################
##                                ##
##                                ##
## Principal Components Analysis  ## 03/07/2019
## Part II.                       ##
##                                ##
####################################

#' Here, I am doing PCA on the data
#' that was already normalized in the carbon_csiaa script.
#' 
# View(aa_norm)
# t.aa <- aa_norm %>% select(-c(Lys, Pro, Ser))
# t.aa <- droplevels(t.aa)
# aa_pca <- vegan::rda(t.aa[,-c(1,2)], scale = TRUE) # this includes all AA samples
# 
# ############
# # Plotting #
# ############
# biplot(aa_pca,
#        display = c("species"),
#        type = c("text"),
#        col = c("#252525", "black"),
#        xlab = "PC1",
#        ylab = "PC2")
# points(aa_pca,
#        pch = c(21,22,23,24,25)[as.numeric(t.aa$Group.ID2)],
#        # col = c("red", "yellow", "blue", "green", "orange")[as.numeric(eaa_norm$`df_eaa[, 1]`)],
#        col = "black",
#        bg = c("red", "yellow", "blue", "darkgreen", "orange")[as.numeric(t.aa$Group.ID2)])
# legend("bottomleft",
#        lty = NULL,
#        bty = "n",
#        pch = c(21,22,23,24,25),
#        legend = levels(t.aa$Group.ID2),
#        cex = 0.45,
#        pt.cex = 0.95,
#        # col = "black",
#        pt.bg = c("red", "yellow", "blue", "darkgreen", "orange"))




# pca1 <- rda(t.norm, scale = TRUE)
# par(pty = "s", bty = "L", xpd = FALSE)
# biplot(pca1,
#        display = c("species"),
#        type = c("text"),
#        col = c("#252525", "black"),
#        xlab = "PC1 (49%)",
#        ylab = "PC2 (28%)")
# points(pca1,
#        pch = c(22),
#        cex = 1.25,
#        col = "black",
#        bg = col8[as.numeric(t.compile$Group.ID5)])
# legend("topright",
#        inset = c(0, 0),
#        lty = NULL,
#        bty = "n",
#        pch = c(22),
#        legend = levels(t.compile$Group.ID5),
#        cex = 0.65,
#        pt.cex = 0.95,
#        pt.bg = c("#91bfdb","#7fbf7b", "brown", "#ffffbf", "#fee090", "#fc8d59", "#d73027"))


# plot(lda_larsen, dimen = 2, abbrev = 1, ylim = c(-6,6), note 3/8/2019 -- keep working on this
#      panel = function(x, y, ...)
#      points(x, y, ...),
#      col = as.numeric(larsen_eaa$Group.ID2), pch = c(21,22,24))

# x <- seq(-7, 5.5, 0.25) http link: https://stat.ethz.ch/pipermail/r-help/2009-August/402378.html
# y <- seq(-4.5, 6.5, 0.25)
# Xcon <- matrix(c(rep(x,length(y)),
#                  rep(y, rep(length(x), length(y)))),,2)


# carbon %>%
#   select(Group.ID, Group.ID2, Phe, Thr, Ile, Leu, Val) %>%
#   filter(Group.ID2 == "Macroalgae" |
#            Group.ID2 == "Microalgae" |
#            Group.ID2 == "Bacteria") -> larsen_train # Modified Larsen 2013 data
# lda_train <- na.omit(lda_train)
# 
# carbon %>%
#   select(Group.ID3, Source, Phe, Thr, Ile, Leu, Val) %>%
#   filter(Group.ID3 == "N2 fixing" |
#            Group.ID3 == "Non-N2 fixing" |
#            Group.ID3 == "Euk microalgae" |
#            Group.ID3 == "Euk macroalgae" |
#            Group.ID3 == "Het bacteria") -> mcmahon_train # Dataset used in McMahon 2015
# mcmahon_train <- droplevels(mcmahon_train)
# lda1 <- MASS::lda(Group.ID3 ~ Phe + Thr + Ile + Leu + Val, data = mcmahon_train)
# plot(lda1)
