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

#' -------------------
#' Keiwgin 1996 data
#' -------------------
#' 

string1 <- '~/Documents/GitHub/data/paleoclimate_data/bc004a-tab.txt'
string2 <- '~/Documents/GitHub/data/paleoclimate_data/bc004d-tab.txt'
keigwin4a <- read.delim(string1, header = TRUE, comment.char = '#')
keigwin4d <- read.delim(string2, header = TRUE, comment.char = '#')

plot(d18Og.rub ~ yrBP, keigwin4a, type= "o")
keigwin4d %>%
  filter(d18Og.rub != '-999') %>%
plot(d18Og.rub ~ yrBP, ., type= "o")
keigwin4d %>%
  filter(carb. != '-999') %>%
  plot(carb. ~ yrBP, ., type= "o")
plot(d15n ~ bp, t.sav, type = "o", col = "gray")
lines(rollmean(d15n, 4, na.pad = TRUE) ~ bp, t.jack, col = "blue")
lines(d15n ~ bp, t.jack, type = "o", col = "blue")
lines(rollmean(d15n, 5, na.pad = TRUE) ~ bp, t.sav)

t.stet$d13c.norm <- t.stet$d13c - mean(t.stet$d13c, na.rm = TRUE)
t.stet$d15n.norm <- t.stet$d15n - mean(t.stet$d15n, na.rm = TRUE)

t.jack <- df.jack
t.jack$bp <- 1950 - t.jack$linear.ad
  
# Richey dataset
fisk <- read.csv('~/Documents/GitHub/data/paleoclimate_data/richey2009-fisk.csv', header = TRUE)
garrison <- read.csv('~/Documents/GitHub/data/paleoclimate_data/richey2009-garrison.csv', header = TRUE)


plot(SST ~ yrBP, fisk, type = "o", xlim = c(0,775))
plot(d13c ~ bp, t.stet, type = "o", xlim = c(0, 775), col = "gray")
lines(rollmean(d13c, 10, na.pad = TRUE) ~ bp, t.stet)

# Schmidt 2012 data
schmidt <- read.delim('~/Documents/GitHub/data/paleoclimate_data/schmidt2012.txt', header = TRUE, comment.char = '#')
plot(sst ~ age_calkaBP, schmidt, type = "o", xlim = c(0, 3))

plot(d13c ~ bp, t.sav, type = "o", xlim = c(0, 3000), col = "gray")
lines(rollmean(d13c, 10, na.pad = TRUE) ~ bp, t.sav)

t.stet$bp <- 1950 - t.stet$linear.ad

# Import Saenger et al (2011) data from Paleoceanography
# Core 1: Core KNR140_2_59GGC
# Core 2: Core CH07_98_MC22
saenger.core1 <- '~/Google Drive/projects/rproj/seus/data/paleoclimate_data/saenger2011_core1.csv'
saenger.core2 <- '~/Google Drive/projects/rproj/seus/data/paleoclimate_data/saenger2011_core2.csv'
core1 <- read.csv(saenger.core1)
core2 <- read.csv(saenger.core2)

t.stet <- df.stet %>% select(linear.ad, d15n)
t.stet <- melt(t.stet, "linear.ad")
t.core1 <- core1 %>% select(Year.AD., SST.Anand.)
t.core1 <- melt(t.core1, "Year.AD.")
t.core2 <- core2 %>% select(Year.AD., SST.Anand.)
t.core2 <- melt(t.core2, "Year.AD.")
t.core2$variable <- factor(t.core2$variable, label = "SST.A")
colnames(t.core1)[colnames(t.core1)=="Year.AD."] <- "linear.ad"
colnames(t.core2)[colnames(t.core2)=="Year.AD."] <- "linear.ad"
stet.core1 <- rbind(t.stet, t.core1)
# proxynames <- c("delta^{15}*N (\u2030)",
#                `SST.Anand.` = "deg Celsius")
stet.core1$variable <- factor(stet.core1$variable, labels = c(expression(paste(delta^{15}*'N (\u2030)')),
                                                              expression("SST 59GGC " (degree~C))))
                                                              #expression("SST MC22 " (degree~C))))
ggplot(data = stet.core1, aes(x = linear.ad, y = value)) +
  geom_line() +
  facet_wrap(~ variable, scales = "free_y",
             strip.position = "left",
             nrow = 3,
             # labeller = as_labeller(proxynames)) +
             labeller = label_parsed) +
  xlab("Year CE") +
  ylab(NULL) +
  theme_classic() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  geom_vline(xintercept=600, color="black", linetype="dashed")

#' Throw in the d13C bulk data from Stetson as well
#' Goal: Figure with 4 charts stacked on top of each other
#' 
t.stet <- bulk.stet %>% select(linear.ad, d15n, d13c)
t.stet <- melt(t.stet, "linear.ad")
stet.cores2 <- rbind(t.stet, t.core1, t.core2)

stet.cores2$variable <- factor(stet.cores2$variable, labels = c(expression(paste(delta^{15}*'N (\u2030)')),
                                                                expression(paste(delta^{13}*'C (\u2030)')),
                                                                expression("SST 59GGC " (degree~C))))
                                                               # expression("SST MC22 " (degree~C))))

#' -------------------------------------------------------------
#' Figure comparing Stetson d13C and d15N to cores are below
#' 
#' -------------------------------------------------------------

ggplot(data = stet.cores2, aes(x = linear.ad, y = value)) +
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
  xlab("Year CE") +
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
t.sav <- df.sav %>% select(linear.ad, d15n, d13c)
t.sav <- melt(t.sav, "linear.ad")
sav.cores <- rbind(t.sav, t.core1, t.core2)
sav.cores$variable <- factor(sav.cores$variable, labels = c(expression(paste(delta^{15}*'N (\u2030)')),
                                                                expression(paste(delta^{13}*'C (\u2030)')),
                                                                expression("SST 59GGC " (degree~C)),
                                                                expression("SST MC22 " (degree~C))))
ggplot(data = sav.cores, aes(x = linear.ad, y = value)) +
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
t.jack <- df.jack %>% select(linear.ad, d15n)
t.jack <- melt(t.jack, "linear.ad")
youngjack.cores <- rbind(t.jack, t.core1)
youngjack.cores$variable <- factor(youngjack.cores$variable, labels = c(expression(paste(delta^{15}*'N (\u2030)')),
                                                                        expression("SST 59GGC " (degree~C))))
                                                                        #expression("SST MC22 " (degree~C))))

ggplot(data = youngjack.cores, aes(x = linear.ad, y = value)) +
  annotate("rect",xmin=1600,xmax=1850,
           ymin=-Inf,ymax=Inf,fill="#9ecae1",color=NA,size=0.25,alpha=0.25) +
  annotate("rect",xmin=950,xmax=1250,
           ymin=-Inf,ymax=Inf,fill="#e9a3c9",color=NA,size=0.25,alpha=0.25) +
  geom_line(color = 'gray') +
  geom_line(aes(y=rollmean(value, 3, na.pad = TRUE)), color = 'black') +
  # geom_point(color = '#1f78b4', shape = 21) +
  facet_wrap(~ variable, scales = "free_y",
             strip.position = "left",
             nrow = 4,
             # labeller = as_labeller(proxynames)) +
             labeller = label_parsed) +
  xlab("Year CE") +
  ylab(NULL) +
  xlim(0,1300) +
  theme_classic() +
  theme(strip.background = element_blank(), strip.placement = "outside")
ggsave("sav_sst_cores.pdf", width=8, height=9)


pl1 <- t.jack %>%
  select(linear.ad, value) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = value), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.75) +
  geom_line(aes(y=rollmean(value, 3, na.pad = TRUE)), size = 0.85, color = '#cc4c02') +
  
  geom_vline(xintercept = 950, lty = 'dashed') +
  geom_vline(xintercept = 1250, lty = 'dashed') +
  geom_vline(xintercept = -900, lty = 'longdash') +
  geom_vline(xintercept = -300, lty = 'longdash') +
  geom_vline(xintercept = -250, lty = "dotted") +
  geom_vline(xintercept = 400, lty = "dotted") +
  
  annotate("text", x = 1100, y = 7, label = 'Medieval Warming', size = 2) +
  annotate("text", x = -600, y = 7, label = "Iron Age Cold Epoch", size = 4) +
  annotate("text", x = 75, y = 7, label = "Roman Warm Period", size = 2) +
  
  ylab(n) +
  xlab(NULL) +
  # xlim(-1200, 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0,1300)) +
  theme_classic() +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_blank(),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pl2 <- t.core1 %>%
  select(linear.ad, value) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = value), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.75) +
  geom_line(aes(y=rollmean(value, 3, na.pad = TRUE)), size = 0.85, color = '#fe9929') +
  
  geom_vline(xintercept = 950, lty = 'dashed') +
  geom_vline(xintercept = 1250, lty = 'dashed') +
  geom_vline(xintercept = -900, lty = 'longdash') +
  geom_vline(xintercept = -300, lty = 'longdash') +
  geom_vline(xintercept = -250, lty = "dotted") +
  geom_vline(xintercept = 400, lty = "dotted") +
  
  ylab(c) +
  theme_classic() +
  xlab(x) +
  # xlim(-1200, 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0,1300)) +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_text(size=10, color = "black"),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

grid.newpage()
grid.draw(rbind(ggplotGrob(pl1), ggplotGrob(pl2), size = "last"))
