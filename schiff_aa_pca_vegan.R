#' Building a principal components analysis
#' with d13C-EAA data
#' 
#' John Schiff
#' 

# ****************************
#
# Using vegan package 
# 
# As used in Larsen et al 2013
#
# ****************************



library(dplyr)
library(MASS)
library(vegan)

######################################
##                                  ##
## Principal Component Analysis     ##
##                                  ##
######################################

carbon %>%
  dplyr::select(Group.ID, Group.ID2, Feature, Source, Phe, Thr, Ile, Leu) %>%
  filter(Group.ID2 == "Microalgae" |
           Group.ID2 == "Macroalgae" |
           Group.ID2 == "Bacteria" |
           Group.ID2 == "Leiopathes" |
           Feature == "Leiopathes-post") -> pca_set

carbon %>%
  dplyr::select(Group.ID, Group.ID2, Group.ID3, Source, Phe, Thr, Ile, Leu) %>%
  filter(Group.ID2 == "Microalgae" |
           Group.ID2 == "Macroalgae" |
           Group.ID2 == "Bacteria") -> lda_set


pca_set <- droplevels(pca_set)
pca_set2 <- pca_set[,c(5:8)] - rowMeans(pca_set[,c(5:8)])
t.pca <- rda(pca_set2, scale = TRUE)
df_pca <- prcomp(pca_set2, center = TRUE, scale. = TRUE)
summary(eigenvals(t.pca)) # Important for getting variance in vegan objects

par(pty = "s")
biplot(t.pca,
       display = c("species"),
       type = c("text"),
       col = c("#252525", "black"),
       xlab = "PC1",
       ylab = "PC2")
points(t.pca,
       pch = c(21,22,23,23,25)[as.numeric(pca_set$Group.ID2)],
       cex = 1.25,
       col = "black",
       bg = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")[as.numeric(pca_set$Group.ID2)])
legend("topright",
       lty = NULL,
       bty = "n",
       pch = c(21,22,23,23,25),
       legend = levels(pca_set$Group.ID2),
       cex = 0.75,
       inset=c(-0.15,0),
       pt.cex = 0.95,
       pt.bg = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"))



##################################
##                              ##
##                              ##
## Linear Discriminant Analysis ##
##                              ##
##                              ##
##################################

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