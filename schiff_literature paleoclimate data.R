#' -----------------------------------------------------
#' Topic: Compiled paleoclimate data from other sources
#' John Schiff
#' 02/08/2019
#' -----------------------------------------------------
#' 
#' This script file is meant for importing, analyzing,
#' and combining paleoclimate data with the isotope data
#' from the black corals and bamboo corals. Different
#' datasets are imported and organized by the study
#' utilizing it.
#' 
#' Other datasets can be utilized. I got these from the
#' NOAA paleoclimate database.
#' 

# Import Saenger et al (2011) data from Paleoceanography
# Core 1: Core KNR140_2_59GGC
# Core 2: Core CH07_98_MC22
saenger.core1 <- '~/Google Drive/projects/rproj/seus/data/paleoclimate_data/saenger2011_core1.csv'
saenger.core2 <- '~/Google Drive/projects/rproj/seus/data/paleoclimate_data/saenger2011_core2.csv'
core1 <- read.csv(saenger.core1)
core2 <- read.csv(saenger.core2)

t.stet <- bulk.stet %>% select(year.ad, d15n.ra)
t.stet <- melt(t.stet, "year.ad")
t.core1 <- core1 %>% select(Year.AD., SST.Anand.)
t.core1 <- melt(t.core1, "Year.AD.")
t.core2 <- core2 %>% select(Year.AD., SST.Anand.)
t.core2 <- melt(t.core2, "Year.AD.")
t.core2$variable <- factor(t.core2$variable, label = "SST.A")
colnames(t.core1)[colnames(t.core1)=="Year.AD."] <- "year.ad"
colnames(t.core2)[colnames(t.core2)=="Year.AD."] <- "year.ad"
stet.core1 <- rbind(t.stet, t.core1, t.core2)
# proxynames <- c("delta^{15}*N (\u2030)",
#                `SST.Anand.` = "deg Celsius")
stet.core1$variable <- factor(stet.core1$variable, labels = c(expression(paste(delta^{15}*'N (\u2030)')),
                                                              expression("SST 59GGC " (degree~C)),
                                                              expression("SST MC22 " (degree~C))))
ggplot(data = stet.core1, aes(x = year.ad, y = value)) +
  geom_line() +
  facet_wrap(~ variable, scales = "free_y",
             strip.position = "left",
             nrow = 3,
             # labeller = as_labeller(proxynames)) +
             labeller = label_parsed) +
  xlab("Calendar Year (C.E.)") +
  ylab(NULL) +
  theme_classic() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  geom_vline(xintercept=600, color="black", linetype="dashed")

#' Throw in the d13C bulk data from Stetson as well
#' Goal: Figure with 4 charts stacked on top of each other
#' 
t.stet <- bulk.stet %>% select(year.ad, d15n.ra, d13c.ra)
t.stet <- melt(t.stet, "year.ad")
stet.cores2 <- rbind(t.stet, t.core1, t.core2)

stet.cores2$variable <- factor(stet.cores2$variable, labels = c(expression(paste(delta^{15}*'N (\u2030)')),
                                                                expression(paste(delta^{13}*'C (\u2030)')),
                                                                expression("SST 59GGC " (degree~C)),
                                                                expression("SST MC22 " (degree~C))))

#' -------------------------------------------------------------
#' Figure comparing Stetson d13C and d15N to cores are below
#' 
#' -------------------------------------------------------------

ggplot(data = stet.cores2, aes(x = year.ad, y = value)) +
  annotate("rect",xmin=1600,xmax=1850,
           ymin=-Inf,ymax=Inf,fill="#9ecae1",color=NA,size=0.25,alpha=0.25) +
  annotate("rect",xmin=950,xmax=1250,
           ymin=-Inf,ymax=Inf,fill="#e9a3c9",color=NA,size=0.25,alpha=0.25) +
  geom_line() +
  facet_wrap(~ variable, scales = "free_y",
             strip.position = "left",
             nrow = 4,
             # labeller = as_labeller(proxynames)) +
             labeller = label_parsed) +
  xlab("Calendar Year (C.E.)") +
  ylab(NULL) +
  theme_classic() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  geom_vline(xintercept=600, color="black", linetype="dashed")
ggsave("stetson_sst_cores.pdf", width=8, height=9)

#' -------------------------------------------------------------
#' Figure comparing Savannah d13C and d15N to cores, following
#' similar code to above with Stetson
#' 
#' -------------------------------------------------------------
#' 
t.sav <- bulk.sav %>% select(year.ad, d15n.ra, d13c.ra)
t.sav <- melt(t.sav, "year.ad")
sav.cores <- rbind(t.sav, t.core1, t.core2)
sav.cores$variable <- factor(sav.cores$variable, labels = c(expression(paste(delta^{15}*'N (\u2030)')),
                                                                expression(paste(delta^{13}*'C (\u2030)')),
                                                                expression("SST 59GGC " (degree~C)),
                                                                expression("SST MC22 " (degree~C))))
ggplot(data = sav.cores, aes(x = year.ad, y = value)) +
  annotate("rect",xmin=1600,xmax=1850,
           ymin=-Inf,ymax=Inf,fill="#9ecae1",color=NA,size=0.25,alpha=0.25) +
  annotate("rect",xmin=950,xmax=1250,
           ymin=-Inf,ymax=Inf,fill="#e9a3c9",color=NA,size=0.25,alpha=0.25) +
  geom_line() +
  facet_wrap(~ variable, scales = "free_y",
             strip.position = "left",
             nrow = 4,
             # labeller = as_labeller(proxynames)) +
             labeller = label_parsed) +
  xlab("Calendar Year (C.E.)") +
  ylab(NULL) +
  theme_classic() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  geom_vline(xintercept=600, color="black", linetype="dashed")
ggsave("sav_sst_cores.pdf", width=8, height=9)

#' -------------------------------------------------------------
#' Figure comparing Jacksonville-4684 d15N to cores, following
#' similar code to above with Stetson and Savannah
#' 
#' -------------------------------------------------------------
#'
t.youngjack <- bulk.youngjack %>% select(year.ad, d15n.ra)
t.youngjack <- melt(t.youngjack, "year.ad")
youngjack.cores <- rbind(t.youngjack, t.core1, t.core2)
youngjack.cores$variable <- factor(youngjack.cores$variable, labels = c(expression(paste(delta^{15}*'N (\u2030)')),
                                                                        expression("SST 59GGC " (degree~C)),
                                                                        expression("SST MC22 " (degree~C))))

ggplot(data = youngjack.cores, aes(x = year.ad, y = value)) +
  annotate("rect",xmin=1600,xmax=1850,
           ymin=-Inf,ymax=Inf,fill="#9ecae1",color=NA,size=0.25,alpha=0.25) +
  annotate("rect",xmin=950,xmax=1250,
           ymin=-Inf,ymax=Inf,fill="#e9a3c9",color=NA,size=0.25,alpha=0.25) +
  geom_line() +
  facet_wrap(~ variable, scales = "free_y",
             strip.position = "left",
             nrow = 4,
             # labeller = as_labeller(proxynames)) +
             labeller = label_parsed) +
  xlab("Calendar Year (C.E.)") +
  ylab(NULL) +
  theme_classic() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  geom_vline(xintercept=600, color="black", linetype="dashed")
ggsave("sav_sst_cores.pdf", width=8, height=9)