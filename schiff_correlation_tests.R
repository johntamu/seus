# Required packages
library(maps)
library(ggplot2)
library(marmap)
library(oceanmap)
library(dplyr)
library(ggplot2)
library(lattice)
library(latticeExtra)
library(reshape2)
library(Hmisc)
library(grid)
library(zoo) 
library(forecast)

# path1 <- '~/Documents/GitHub/data/schiff_bulk_years_12-22-2019.csv'
path1 <- '~/Documents/GitHub/data/schiff_bulk_years_09-04-2019.csv'
# path2 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff_bulk_years_08-14-2019.csv'
# path3 <- '/home/john/Desktop/data/schiff_bulk_years_08-14-2019.csv'

path <- path1

df.bulk <- read.csv(path, header = TRUE)
# View(df.bulk)

colnames(df.bulk)[names(df.bulk) == "distance..mm."] <- "distance" # Rename some columns for easier coding
colnames(df.bulk)[names(df.bulk) == "d15n.vs.air"] <- "d15n"
colnames(df.bulk)[names(df.bulk) == "d13c.vs.vpdb"] <- "d13c"

df.jack <- df.bulk %>% filter(coral.id == 'jack-4907-bc1-d3')
df.sav <- df.bulk %>% filter(coral.id == 'sav-4902-bc1-unk')
df.stet <- df.bulk %>% filter(coral.id == 'stet-4904-bc1-d2')
df.jack2 <- df.bulk %>% filter(coral.id == 'jack-4907-bc1-d1')
df.jack4684 <- df.bulk %>% filter(coral.id == 'jack-4684-bc-unk')
df.jack4686 <- df.bulk %>% filter(coral.id == 'jack-4686-bc-d1-t1')


df.stet2 <- df.stet %>% 
  filter(linear.ad < 1850)
df.stet3 <- df.stet %>% 
  filter(linear.ad > 1850)

test.jack4684 <- df.jack4684 %>%
  filter(linear.ad > 1850)


df.jack <- df.jack %>%
  filter(linear.ad > -251)
  
mod1 <- lm(d15n ~ d13c, df.jack)
mod2 <- lm(d15n ~ d13c, df.sav)
mod3 <- lm(d15n ~ d13c, df.stet)

summary(mod1)
summary(mod2)
summary(mod3)
summary(lm(d15n ~ d13c, df.stet2))
summary(lm(d15n ~ d13c, df.stet3))

summary(lm(d15n ~ linear.ad, test.jack4684))
summary(lm(d13c ~ linear.ad, test.jack4684))
summary(lm(d15n ~ linear.ad, df.stet3))
summary(lm(d13c ~ linear.ad, df.stet3))

plot(d15n ~ linear.ad, test.jack4684)
abline(lm(d15n ~ linear.ad, test.jack4684))

plot(d15n ~ linear.ad, df.stet3)
abline(lm(d15n ~ linear.ad, df.stet3))

plot(d13c ~ linear.ad, df.stet3)
abline(lm(d13c ~ linear.ad, df.stet3))

jack4686 <- bulkdata %>% filter(coral.id == 'jack-4686-bc-d1-t1')
jack4686 <- jack4686 %>% filter(linear.ad > 1849)

plot(d15n.vs.air ~ linear.ad, jack4686)
summary(lm(d15n.vs.air ~ linear.ad, jack4686))
summary(lm(d13c.vs.vpdb ~ linear.ad, jack4686))

par(mfrow = c(1,3), bty = 's')
plot(d15n ~ d13c, df.jack,
     pch = 21, bg = 'white', col = "black", 
     ylab = n, xlab = c)
abline(mod1)
plot(d15n ~ d13c, df.sav,
     pch = 21, bg = 'white', col = "black", 
     ylab = n, xlab = c)
abline(mod2)
plot(d15n ~ d13c, df.stet,
     pch = 21, bg = 'white', col = "black",
     ylab = n, xlab = c)
abline(mod3)

plot(d15n ~ d13c, df.stet3,
     pch = 21, bg = 'white', col = "black",
     ylab = n, xlab = c)
abline(lm(d15n ~ d13c, df.stet3))

testmodel <- lm(X14C.Age ~ Distance.microns, r.jack)
summary(testmodel)
plot(X14C.Age ~ Distance.microns, r.jack)
abline(lm(X14C.Age ~ Distance.microns, r.jack))
t.jack <- df.jack
t.jack$Distance.microns <- t.jack$distance*1000
predicted <- predict(testmodel, t.jack)
t.jack$predicted <- predicted
plot(d15n ~ predicted, t.jack)
tt.jack <- r.jack

calibrated <- read.csv('~/Documents/GitHub/data/jack4907_calibages.csv', header = TRUE)
tt.jack$calibrated <- calibrated$median
testmodel <- lm(calibrated ~ Distance.microns, tt.jack)
summary(testmodel)
plot(calibrated ~ Distance.microns, tt.jack)
abline(lm(calibrated ~ Distance.microns, tt.jack))
predicted <- predict(testmodel, t.jack, interval = "confidence")
t.jack <-cbind(t.jack, predicted)
plot(d15n ~ lwr, t.jack)

summary(lm(X14C.Age ~ Distance.microns, r.stet))
plot(X14C.Age ~ Distance.microns, r.stet)
abline(lm(X14C.Age ~ Distance.microns, r.stet))

summary(lm(X14C.Age ~ Distance.microns, r.stet2))
plot(X14C.Age ~ Distance.microns, r.stet2)
abline(lm(X14C.Age ~ Distance.microns, r.stet2))

# Figure
library(Hmisc)
c14 <- expression(Delta^{14}*"C (\u2030)")
dis <- expression(paste("Distance from edge", " (",mu,"m)"))
par(mfrow = c(2,2))
pch = 22
plot(X14C.Age ~ Distance.microns, r.stet, type = "o", pch = pch, col = 'black', bg = "black", tck = -0.0275, ylab = "Age (years BP)", xlab = dis)
minor.tick(nx=4, ny=4, tick.ratio=0.25)
abline(lm(X14C.Age ~ Distance.microns, r.stet))
# plot(X14C.Age ~ Distance.microns, r.stet2, type = "o", pch = pch, col = 'black', bg = "black", tck = -0.0275, ylab = "Age (years BP)", xlab = dis)
# minor.tick(nx=4, ny=4, tick.ratio=0.25)
# abline(lm(X14C.Age ~ Distance.microns, r.stet2))
plot(X14C.Age ~ Distance.microns, r.jack, type = "o", pch = pch, col = 'black', bg = "black", tck = -0.0275, ylab = "Age (years BP)", xlab = dis)
minor.tick(nx=4, ny=4, tick.ratio=0.25)
abline(lm(X14C.Age ~ Distance.microns, r.jack))
r.jack2 %>% filter(Distance.microns > 600) -> filtered
plot(X14C.Age ~ Distance.microns, filtered, type = "o", pch = pch, col = 'black', bg = "black", tck = -0.0275, ylab = "Age (years BP)", xlab = dis)
minor.tick(nx=4, ny=4, tick.ratio=0.25)
abline(lm(X14C.Age ~ Distance.microns, filtered))
plot(X14C.Age ~ Distance.microns, r.sav, type = "o", pch = pch, col = 'black', bg = "black", tck = -0.0275, ylab = "Age (years BP)", xlab = dis)
minor.tick(nx=4, ny=4, tick.ratio=0.25)
abline(lm(X14C.Age ~ Distance.microns, r.sav))

filt <- radiocarbon %>% filter(X14C.Age > 0)

library(dplyr)
first <- r.jack %>%
  filter(Distance.microns < 900)
second <- r.jack %>%
  filter(Distance.microns > 900 & Distance.microns < 3500)
third <- r.jack %>%
  filter(Distance.microns > 3500)

summary(lm(X14C.Age ~ Distance.microns, first))
plot(X14C.Age ~ Distance.microns, first)
abline(lm(X14C.Age ~ Distance.microns, first))

summary(lm(X14C.Age ~ Distance.microns, second))
summary(lm(X14C.Age ~ Distance.microns, third))

summary(lm(X14C.Age ~ Distance.microns, r.sav))


carbon %>%
  # dplyr::select(Group.ID, Group.ID2, Group.ID3, Source, Phe, Thr, Ile, Leu) %>%
  filter(Group.ID2 == "Leiopathes") -> leio_aa

summary(lm(Leu ~ Bulk, leio_aa))
plot(Leu ~ Bulk, leio_aa)

# Filter the C-CSIAA data and then create a table for use in dissertation
carbon %>%
  dplyr::select(Group.ID, Group.ID2, Ala, Asx, Glx, Gly, Ser, Pro, Ile, Leu, Lys, Phe, Thr, Tyr, Val) -> newtable
write.csv(newtable, "newtable.csv")


summary(lm(Phe ~ SrcAA, seus))
plot(Phe ~ SrcAA, seus)
abline(lm(Phe~SrcAA,seus))


seus$TP <- (((seus$Glx - seus$Phe)-3.4)/7.6) + 1
seus %>%
  filter(Coral == "Jacksonville Lithoherms-4907" | Coral == "Jacksonville Lithoherms-4684" | Coral == "Savannah Banks-4902") -> corals



phe.aov <- aov(Phe ~ Coral, corals) # One-way ANOVA with Phe
sraa.aov <- aov(SrcAA ~ Coral, corals) # One-way ANOVA with avg Sr-AA
tp.aov <- aov(TP ~ Coral, corals)
sumv.aov <- aov(Sum.V ~ Coral, corals)

summary(phe.aov)
summary(sraa.aov)
summary(tp.aov)
summary(sumv.aov)

gom %>% filter(Disk == 1 | Disk == 7) -> gom
gom$SrcAA <- ((gom$Ser + gom$Phe + gom$Gly)/3)
summary(lm(Phe ~ SrcAA, gom))

plot(Phe ~ SrcAA, gom)
abline(lm(Phe ~ SrcAA, gom))

carbon %>%
  dplyr::select(Group.ID2, Phe, Thr, Ile, Leu) %>%
  filter(Group.ID2 == "Microalgae" |
           Group.ID2 == "Macroalgae" |
           Group.ID2 == "Bacteria" |
           Group.ID2 == "Leiopathes") -> test_table

phe <- test_table$Phe
thr <- test_table$Thr
ile <- test_table$Ile
leu <- test_table$Leu

res.man <- manova(cbind(phe, thr, ile, leu) ~ Group.ID2, data = test_table)
summary(res.man)
summary.aov(res.man)

sepl <- iris$Sepal.Length
petl <- iris$Petal.Length

bind <- cbind(sepl, petl)

# Looks good, but will need to make sure the samples are normalized to the means first
print(bind)


# CLUSTER ANALYSIS #

carbon %>%
  dplyr::select(Group.ID, Group.ID2, Phe, Thr, Ile, Leu) %>%
  filter(Group.ID2 == "Microalgae" |
           Group.ID2 == "Macroalgae" |
           Group.ID2 == "Bacteria" |
           Group.ID2 == "Leiopathes") -> test_table
mydata <- test_table
mydata <- data.frame(mydata[,-1], row.names=mydata[,1])
mydata <- mydata %>% dplyr::select(Phe, Thr, Ile, Leu)
mydata <- scale(mydata)
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
fit <- kmeans(mydata, 4) # 4 cluster solution
# get cluster means
aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(mydata, fit$cluster)

df <- USArrests
df <- na.omit(df)
df <- scale(df)

library(cluster)
clusplot(mydata, fit$cluster, color=TRUE, shade=TRUE,
         labels=2, lines=0)

mydata2 <- cbind(mydata, test_table$Group.ID2)

head(iris)
df <- iris
df <- as.data.frame(iris)
row.names(df) <- paste(df$Species, row.names(df), sep="_") 
df$Species <- NULL
head(df)


# t-test before and after 1900 d15N

t.table <- read.csv('~/Documents/GitHub/data/schiff_bulk_years_06-17-2020.csv')
library(dplyr)

t.table <- t.table %>% filter(coral.id == "jack-4684-bc-unk" | coral.id == 'stet-4904-bc1-d2')
t.table1 <- t.table %>% filter(coral.id == "jack-4684-bc-unk")
t.table2 <- t.table %>% filter(coral.id == "stet-4904-bc1-d2")

t.jack <- t.table %>% filter(coral.id == "jack-4907-bc1-d3")
t.jack <- t.jack %>% filter(linear.ad > -750 & linear.ad < -450)
lm1 <- lm(d15n.vs.air ~ d13c.vs.vpdb,data=t.jack)
summary(lm(lm1))
plot(d15n.vs.air ~ d13c.vs.vpdb,data=t.jack)

# statistical tests
t.table1 <- t.table %>% filter(coral.id == "jack-4686-bc-d1-t1")
t.table2 <- t.table %>% filter(coral.id == "jack-4907-bc1-d3")
t.table3 <- t.table %>% filter(coral.id == "stet-4904-bc1-d2")
t.table4 <- t.table %>% filter(coral.id == "sav-4902-bc1-unk")

lm1 <- lm(d15n.vs.air ~ d13c.vs.vpdb,data=t.table1)
lm2 <- lm(d15n.vs.air ~ d13c.vs.vpdb,data=t.table2)
lm3 <- lm(d15n.vs.air ~ d13c.vs.vpdb,data=t.table3)
lm4 <- lm(d15n.vs.air ~ d13c.vs.vpdb,data=t.table4)

summary(lm1)
plot(d15n.vs.air ~ d13c.vs.vpdb,data=t.table1)
abline(lm1)
summary(lm2)
plot(d15n.vs.air ~ d13c.vs.vpdb,data=t.table2) # jack4907
abline(lm2)
summary(lm3)
plot(d15n.vs.air ~ d13c.vs.vpdb,data=t.table3) # stetson4904
abline(lm3)
summary(lm4)
plot(d15n.vs.air ~ d13c.vs.vpdb,data=t.table4) # sav4902
abline(lm4)


stet1 <- t.table2 %>% filter(linear.ad > 1900)
stet2 <- t.table2 %>% filter(linear.ad < 1900 & linear.ad > 1500)

stet3 <- t.table2 %>% filter(linear.ad > 700 & linear.ad < 900)

t.test(stet1$d15n.vs.air, stet2$d15n.vs.air)

j1 <- t.table1 %>% filter(linear.ad > 1900)
j2 <- t.table1 %>% filter(linear.ad < 1900 & linear.ad > 1500)

t.test(j1$d15n.vs.air, j2$d15n.vs.air)








# Mol percent calculation and figure with error bar
# mol <- read.csv('C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff_stetson aa mol percent.csv')
mol <- read.csv('~/Documents/GitHub/data/schiff_stetson aa mol percent.csv')

mol$Total <- rowSums(mol[,2:13])
# mol$Year.AD <- seus_carbon$Year.AD

mean(mol$Total)
sd(mol$Tota)

data <- data.frame(
  name=letters[1:5],
  value=sample(seq(4,15),5),
  sd=c(1,0.2,3,2,4))

ggplot(data) +
  geom_bar( aes(x=name, y=value), stat="identity", fill="skyblue", alpha=0.5) +
  geom_linerange( aes(x=name, ymin=value-sd, ymax=value+sd), colour="orange", alpha=0.9, size=1.3)

library(dplyr)
library(reshape2)
library(lattice)

mol %>%
  melt(., "Sample") -> t.data
write.csv(t.data, "data.csv")
levels(t.data$variable)
mol %>%
  melt(., "Sample") %>%
  barchart(value ~ variable,
           data = .,
           horiz = FALSE,
           ylab = 'Amino acids (mol, %)',
           col = "darkgrey")

t.mol <- mol
rownames(t.mol) <- t.mol[,1]
t.mol$Sample <- NULL
mean <- colMeans(t.mol)
sd <- apply(t.mol, 2, sd)
mean <- as.data.frame(mean)
mean$sd <- sd
library(data.table)
mean <- setDT(mean, keep.rownames = "aa")[]
ggplot(mean) +
  geom_bar( aes(x=aa, y=mean), stat="identity", fill="#1b9e77", alpha=1, color = "black", size = 0.25) +
  geom_linerange( aes(x=aa, ymin=mean-sd, ymax=mean+sd), colour="#000000", alpha=1, size=0.5) +
  theme_linedraw() +
  theme(panel.grid = element_blank(), axis.ticks.x = element_line(size = 0.25, color = "black"),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) +
  ylab("Amino acid (Mol%)") +
  xlab(NULL)
ggsave("mol_percent.png",width = 8, height = 5)


library(dplyr)
ndata <- read.csv('~/Documents/GitHub/data/cleaned_ndata_06-17-2020.csv')
ndata %>% filter(Region == 'SEUS') -> table
table %>% filter(Phe < 13) -> table
ancient <- table %>% filter(Year.CE < 1500)
ancient <- table %>% filter(Year.CE > 1000)
plot(Phe ~ Year.CE, table)
abline(lm(Phe ~ Year.CE, table))
plot(Phe ~ Bulk, ancient)
summary(lm(Phe ~ Year.CE, table))
summary(lm(Phe ~ Year.CE, ancient))
summary(lm(Phe ~ Bulk, ancient))
summary(lm(Glu ~ Year.CE, table))

phe.aov <- aov(Phe ~ Sample.ID2, table) # One-way ANOVA with Phe
gly.aov <- aov(Gly ~ Sample.ID2, table)
ser.aov <- aov(Ser ~ Sample.ID2, table)
glu.aov <- aov(Glu ~ Sample.ID2, table)
sraa.aov <- aov(SrcAA ~ Coral, corals) # One-way ANOVA with avg Sr-AA
tp.aov <- aov(TP ~ Coral, corals)
sumv.aov <- aov(Sum.V ~ Coral, corals)

bulkdata <- read.csv('~/Documents/GitHub/data/schiff_bulk_years_06-21-2020.csv')
colnames(bulkdata)[names(bulkdata) == "distance..mm."] <- "distance" # Rename some columns for easier coding
colnames(bulkdata)[names(bulkdata) == "d15n.vs.air"] <- "d15n"
colnames(bulkdata)[names(bulkdata) == "d13c.vs.vpdb"] <- "d13c"
colnames(bulkdata)[names(bulkdata) == "linear.ad"] <- "yrAD"

df.jack <- bulkdata %>% filter(coral.id == 'jack-4907-bc1-d3')
df.sav <- bulkdata %>% filter(coral.id == 'sav-4902-bc1-unk')
df.stet <- bulkdata %>% filter(coral.id == 'stet-4904-bc1-d2')
df.jack2 <- bulkdata %>% filter(coral.id == 'jack-4907-bc1-d1')
df.jack4684 <- bulkdata %>% filter(coral.id == 'jack-4684-bc-unk')
df.jack4686 <- bulkdata %>% filter(coral.id == 'jack-4686-bc-d1-t1')

plot(d15n ~ d13c, df.jack)
abline(lm(d15n ~ d13c, df.jack2))
plot(d15n ~ d13c, df.jack2)
plot(d15n ~ d13c, df.jack4686)
plot(d15n ~ d13c, df.stet)
summary(lm(d15n ~ d13c, df.jack))
summary(lm(d15n ~ d13c, df.jack4686))
summary(lm(d15n ~ d13c, df.jack2))
summary(lm(d15n ~ d13c, df.stet))

df.jack <- df.jack %>% filter(yrAD < 300)
df.jack <- df.jack %>% filter(yrAD > -250)
df.jack <- df.jack %>% filter(yrAD > 400 & yrAD < 1500 )
df.sav <- df.sav %>% filter(yrAD > 400 & yrAD < 1500 )
plot(d15n ~ yrAD, df.jack) 
summary(lm(d15n ~ yrAD, df.jack))   
summary(lm(d15n ~ yrAD, df.sav))     

plot(d15n ~ d13c, df.jack)
ggplot(table, aes(x=Year.CE, y=Phe))+
  geom_point()+
  geom_smooth(method=lm, se=TRUE)
ggplot(df.jack, aes(x=d13c, y=d15n))+
  geom_point()+
  geom_smooth(method=lm, se=TRUE)
