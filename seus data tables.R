#' ------------------------------
#' Tables for SEUS Project
#' John Schiff
#' 02/07/2019
#' ------------------------------
#' 

library(stargazer)

#'
#' Radiocarbon data table
#' 

# \caption{$\delta$$^{13}$C$_{AA}$ for each Coral.} <- Latex syntax for special characters
path1 <- '~/Google Drive/projects/rproj/seus/data/schiff radiocarbon 02-05-2019.csv'
path2 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff radiocarbon 02-05-2019.csv'
radiocarbon <- read.csv(path1)
radiotable <- radiocarbon
levels(radiotable$Coral) <- c('Jacksonville-4684', 'Jacksonville-4907', 'Savannah Banks-4902','Stetson-old','Stetson Banks-4904')

View(levels(radiotable$Coral))
radiotable$Distance <- radiotable$Distance..um./1000
radiotable$Distance..um. <- NULL
radiotable <- radiotable %>% select(Coral, Sample, Distance, Fraction.modern, error, D14C, d13C, X14C.Age, error.1)
radiotable <- radiotable %>% dplyr::filter(Coral != 'Stetson-old')

# Output can be LaTeX or HTML. For copying into Word, do HTML
stargazer(radiotable, summary = FALSE, rownames = FALSE, digits = 1, digit.separator = "",
          type = "html", out = "radiocarbon.html")
write.table(radiotable, file = "radiocarbon.txt", sep = ",", quote = FALSE, row.names = FALSE) # I don't like doing this method with larger tables

#'
#' Carbon CSIAA data for each coral
#' 

seusDataC <- trainData %>% dplyr::filter(Group.ID2 == "Leiopathes")
seusDataC$Peel <- c(7, 22, 40, 109)
seusDataC$Coral <- c(rep("Stetson Banks-4904"))
seusDataC <- seusDataC %>% dplyr::select(Coral, Peel, Ala, Asx, Glx, Gly, Ile, Leu, Lys, Phe, Pro, Ser, Thr, Tyr, Val)

stargazer(seusDataC, summary = FALSE, rownames = FALSE, digits = 1, digit.separator = "",
          type = "html", out = "carbon csia.html")

#'
#' Nitrogen CSIAA data for each coral
#' 
# t.df <- read.csv("C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff ncsiaa 02-06-2019.csv")

t.df <- read.csv('~/Google Drive/projects/rproj/seus/data/schiff ncsiaa 02-06-2019.csv')
t.df <- t.df %>% select(Coral, sample.no., Ala, Asx, Glx, Gly, Ile, Leu, Phe, Pro, Ser, Val)
stargazer(t.df, summary = FALSE, rownames = FALSE, digits = 1, type = "html", out = "nitrogen csia.html")

#' 
#' Table showing specimens collected and where
#' Sample ID, Study Site, Lat/Long, Depth (m)
#' 

specimens <- data.frame("Coral" = c('JACK4907-BC1', 'JACK4684-BC1', 'JACK4686-BC1', 'JACK4907-BAM1','SAV4902-BC1','STET4904-BC1'),
                        "Site" = c("Jacksonville Lithoherms", "Jacksonville Lithoherms", "Jacksonville Lithoherms", "Jacksonville Lithoherms","Savannah Banks", "Stetson Banks"),
                        "Taxon" = c("Black coral", "Black coral", "Black coral", "Bamboo coral", "Black coral", "Black coral"),
                        "Year Collected" = c(2005, 2004, 2004, 2005, 2005, 2005),
                        "Latitude" = NA,
                        "Longitude" = NA,
                        "Depth (m)" = c(538, 561, "TBD", "TBD", 516, 685))
stargazer(specimens, summary = FALSE, rownames = FALSE, type = "html", digit.separator = "", out = "specimens.html")

#' 
#' Table showing inner and outer ages for each black coral specimen
#' 
#' 

ages <- data.frame("Coral" = c('JACK4907-BC1', 'JACK4684-BC1', 'JACK4686-BC1', 'JACK4907-BAM1','SAV4902-BC1','STET4904-BC1'),
                   "14C Inner" = NA,
                   "14C Outer" = NA,
                   "Life span" = NA,
                   "Linear growth rate" = NA)

#' ------------------------------------------------------------
#' Table showing the linear regression output for d15N-Bulk
#' vs. different d15N-AA parameters (SumV, Source AAs, Trophic
#' AAs, Trophic Position)
#' ------------------------------------------------------------
csiaa.all$Bulk.N <- csiaa$bulk.N # Make sure csiaa.all object is loaded from gom_ncsiaa.R file

lm.bulk.sumv <- lm(Bulk.N ~ Sum.V, data = csiaa.all)
lm.bulk.src <- lm(Avg.Src ~ Bulk.N, data = csiaa.all)
lm.bulk.tr <- lm(Avg.Tr ~ Bulk.N, data = csiaa.all)
lm.bulk.tp <- lm(Trophic.Position ~ Bulk.N, data = csiaa.all)
lm.bulk.phe <- lm(Bulk.N ~ Phe, data = csiaa.all)
lm.bulk.gly <- lm(Gly ~ Bulk.N, data = csiaa.all)

lm1 <- lm(Bulk.N ~ Ala + Asx + Glx + Ile + Leu +
            Pro + Val, data = csiaa.all)
plot(lm.bulk.sumv)

plot(Trophic.Position ~ Bulk.N, data = csiaa.all)
summary(lm1)

summary(lm(Bulk.N ~ Asx, data = csiaa.all))

##########################
## Specimens Data Table ##
##########################

coralSpec <- data.frame(
  # specimen.id = c('jack-4907-bc1-d3',
  #                 'jack-4907-bc1-d1',
  #                 'jack-4686-bc1-d1', # BC1 as far as I know and I am calling the disk I sampled Disk 1
  #                 'jack-4907-bam1',
  #                 'stet-4904-bc1-d2',
  #                 'sav-4902-bc1-unk',
  #                 'jack-4684-bc-unk'),
  # coral.id = c('jack-4907-bc1-d3',
  #                'jack-4907-bc1-d1',
  #                'jack-4686-bc1-d1',
  #                'jack-4907-bam1',
  #                'stet-4904-bc1-d2',
  #                'sav-4902-bc1-unk',
  #                'jack-4684-bc-unk'),
  coral.id = c('A',
               'B',
               'D',
               'E',
               'F',
               'G',
               'H'),
  location = c('Jacksonville Lithoherms',
               'Jacksonville Lithoherms',
               'Jacksonville Lithoherms',
               'Jacksonville Lithoherms',
               'Stetson Banks',
               'Savannah Banks',
               'Jacksonville Lithoherms'),
  group = c('Black Coral',
            'Black Coral',
            'Black Coral',
            'Bamboo Coral',
            'Black Coral',
            'Black Coral',
            'Black Coral'
  ),
  
  disk = c(3, 1, 1, 1, 2, NA, NA),
  c14.age.inner = c(3130, 3130, NA, NA, 1910, 2530, 940),
  error.inner = c(35, 40, NA, NA, 40, 40, 40),
  c14.age.outer = c(1010, 1010, NA, NA, NA, 1100, NA),
  error.outer = c(35, 35, NA, NA, NA, NA, 40),
  c14.no.samples = c(18, NA, NA, NA, 4, 7, 32),
  facility = c('NOSAMS',
               NA, 
               NA, 
               NA, 
               'NOSAMS/CAMS',
               'Beta Analytic',
               'NOSAMS/Beta Analytic')
  
)

View(coralSpec)

############################## 3/3/2019
## Lifespan vs Growth rate  ##
##############################
gr.table <- data.frame(coral = c("Jacksonville-4907 Disk 1",
                                 "Jacksonville-4684 Disk *",
                                 "Stetson-4904 Disk 5",
                                 "Stetson-4904 Disk *",
                                 "Savannah-4902 Disk *"),
                   lifespan = c(life.jack, life.jack2, life.stet, life.stet2, life.sav),
                   growth.rate = c(gr.jack, gr.jack2, gr.stet, gr.stet2, gr.sav),
                   growth.rate2 = NA,
                   bomb.growth.rate = NA)

################################
## Carbon C-CSIAA data table  ##
################################
dt <- carbon %>% filter(Group.ID2 == "Leiopathes") %>%
  select(Ala, Asx, Glx, Gly, Ile, Leu, Lys, Phe, Pro, Ser, Thr, Tyr, Val)

stargazer(dt, summary = FALSE, rownames = FALSE, type = "html", out = "carbontable.html")

####################################
## Bulk d15N and d13C data table  ##
## for all samples                ##
####################################

bulk.table <- df.bulk %>%
  filter(coral.id == "jack-4684-bc-unk" | coral.id == "jack-4907-bc1-d1" |
           coral.id == "jack-4907-bc1-d3" | coral.id == "sav-4902-bc1-unk" |
           coral.id == "stet-4904-bc1-d2")

bulk.table <- bulk.table %>% 
  select(coral.id, sample.no., distance, d15n, d13c, linear.ad)

stargazer(bulk.table, summary = FALSE, rownames = FALSE, type = "html", out = 'bulktable.html')

## --------------------
## Data table showing which age model
## was used for each disk
## --------------------

# models <- data.frame(coral = c("Jacksonville-4907 Disk 1",
#                                "Jacksonville-4684 Disk *",
#                                "Stetson-4904 Disk 5",
#                                "Stetson-4904 Disk *",
#                                "Savannah-4902 Disk *"),
#                      linear.age.model)
