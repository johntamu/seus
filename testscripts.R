xyplot(c.stet1$d13c.vs.vpdb ~ c.stet1$sample.no.,
       groups = NULL,
       labels = c.stet1$sample.no,
       auto.key=list(space="right"),
       panel = function(x, y, groups, subscripts, labels, do.labels, ...) {
         panel.xyplot(x, y, groups = groups, subscripts = subscripts, ...)
                      labels <- labels[subscripts]
                      do.labels <- do.labels[subscripts]
                      panel.text(x[do.labels], y[do.labels],
                                 labels = labels[do.labels],
                                 pos = 1, cex=0.5, ...)
       })


plot(c.stet1$d13c.vs.vpdb ~ c.stet1$sample.no., xlim=c(200,400))
with(c.stet1, 
     text(c.stet1$d13c.vs.vpdb ~ c.stet1$sample.no., 
          labels = c.stet1$sample.no, pos = 4))

plot(c.stet1$d15n.vs.air ~ c.stet1$sample.no., xlim=c(300,400))
with(c.stet1, 
     text(c.stet1$d15n.vs.air ~ c.stet1$sample.no., 
          labels = c.stet1$sample.no, pos = 4))

plot(c.jack1$d13c.vs.vpdb ~ c.jack1$clam.ad, type = "l")
with(c.jack1, 
     text(c.jack1$d13c.vs.vpdb ~ c.jack1$clam.ad, 
          labels = c.jack1$sample.no, pos = 4))


plot(c.jack1$d13c.vs.vpdb ~ c.jack1$clam.ad,
     xlim=c(-1000, -400))
with(c.jack1, 
     text(c.jack1$d13c.vs.vpdb ~ c.jack1$clam.ad, 
          labels = c.jack1$sample.no, pos = 4))
# symbols(x=-800, circles=rep(1,2), add=T, inches=F)

library(ggplot2)
sample <- c.jack1$sample.no.
ggplot(data=c.jack1, aes(clam.ad, d13c.vs.vpdb)) +
  theme_bw() +
  geom_point() +
  geom_point(data=c.jack1[sample == "232",],
             pch=21, fill=NA, size=4, colour="red", stroke=1) +
  geom_point(data=c.jack1[sample == "122",],
             pch=21, fill=NA, size=4, colour="red", stroke=1) +
  geom_point(data=c.jack1[sample == "123",],
             pch=21, fill=NA, size=4, colour="red", stroke=1) +
  geom_point(data=c.jack1[sample == "142",],
             pch=21, fill=NA, size=4, colour="red", stroke=1) +
  geom_point(data=c.jack1[sample == "143",],
             pch=21, fill=NA, size=4, colour="red", stroke=1) +
  geom_point(data=c.jack1[sample == "19",],
             pch=21, fill=NA, size=4, colour="red", stroke=1) +
  geom_point(data=c.jack1[sample == "9",],
             pch=21, fill=NA, size=4, colour="red", stroke=1)

sample.n <- c.stet1$sample.no.
ggplot(data=c.stet1, aes(clam.ad, d15n.vs.air)) +
  theme_bw() +
  geom_point(shape=21) +
  geom_text(data=c.stet1, aes(label=sample.no.), hjust=0, vjust=-1) +
  geom_line() +
  xlim(500, 2000) +
  geom_point(data=c.stet1[sample.n == "343" | sample.n == "344"
                          | sample.n == "345"
                          | sample.n == "336" | sample.n == "337" | sample.n == "338"
                          | sample.n == "366" | sample.n == "367" | sample.n == "368"
                          | sample.n == "379"
                          | sample.n == "267",],
             pch=21, fill=NA, size=4, colour="red", stroke=1)

# These are the samples I sent to Stephanie Christensen last year
ggplot(data=c.stet1, aes(clam.ad, d13c.vs.vpdb)) +
  theme_bw() +
  geom_point(shape=21) +
  geom_text(data=c.stet1, aes(label=sample.no.), hjust=0, vjust=-1) +
  geom_line() +
  xlim(1400, 2000) +
  geom_point(data=c.stet1[sample.n == "7" | sample.n == "8"
                          | sample.n == "22"
                          | sample.n == "25" | sample.n == "40" | sample.n == "46"
                          | sample.n == "103" | sample.n == "109" | sample.n == "135",],
             pch=21, fill=NA, size=4, colour="red", stroke=2) +
  labs(x="Calendar Year (C.E.)",
       y=expression(delta^{13}*"C (\u2030)"))


pro.time <- ggplot(csiaa, aes(x=year.ad, y=Pro)) +
  annotate("rect",xmin=950,xmax=1250,
           ymin=-Inf,
           ymax=Inf,fill='#e9a3c9',color=NA,size=0.25,alpha=0.25) +
  # Little Ice Age
  annotate("rect",xmin=1550,xmax=1850,
           ymin=-Inf,
           ymax=Inf,fill="#9ecae1",color=NA,size=0.25,alpha=0.25) +
  # Bond event 1
  annotate("rect",xmin=595,xmax=605,
           ymin=-Inf,
           ymax=Inf,fill="black",color="black",size=0.25,alpha=0.25) +
  # label for MCA
  annotate("text", x=1100, y=13.5, label="MCA") +
  # label for Little Ice Age
  annotate("text", x=1725, y=13.5, label="LIA") +
  # label for Bond Event
  annotate("text", x=210, y=13, label="Bond event 1") +
  geom_point(aes(color=coral.id), 
             size=4.5) +
  geom_point(shape = 21, color="black", size=4.5) +
  scale_color_brewer(palette = "Dark2",
                     name = "Specimen",
                     breaks = c("4907bc1-d1",
                                "4902bc1-d?",
                                "4684bc1-d?"),
                     labels = c(" Jacksonville (4907)",
                                " Savannah (4902)",
                                " Jacksonville (4684)")) +
  theme_linedraw() +
  # theme(panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank()) +
  labs(x="Calendar Year (C.E.)",
       y=expression(delta^{15}*"N Phe")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # Medieval Climate Anomaly, warming
  geom_segment(x=200, y=13.5,xend=585, yend=13.75) +
  ylim(20,35) +
  scale_x_continuous(limits=c(-750, 2005),
                     breaks=seq(-750, 2005, by=550))

pro.time

# traa.time <- csiaa.long %>% filter(variable == 'traa')
traa.time <- ggplot(csiaa, aes(x=year.ad, y=traa)) +
  annotate("rect",xmin=950,xmax=1250,
           ymin=-Inf,
           ymax=Inf,fill='#e9a3c9',color=NA,size=0.25,alpha=0.25) +
  # Little Ice Age
  annotate("rect",xmin=1550,xmax=1850,
           ymin=-Inf,
           ymax=Inf,fill="#9ecae1",color=NA,size=0.25,alpha=0.25) +
  # Bond event 1
  annotate("rect",xmin=595,xmax=605,
           ymin=-Inf,
           ymax=Inf,fill="black",color="black",size=0.25,alpha=0.25) +
  # label for MCA
  annotate("text", x=1100, y=13.5, label="MCA") +
  # label for Little Ice Age
  annotate("text", x=1725, y=13.5, label="LIA") +
  # label for Bond Event
  annotate("text", x=210, y=13, label="Bond event 1") +
  geom_point(aes(color=coral.id), 
             size=4.5) +
  geom_point(shape = 21, color="black", size=4.5) +
  scale_color_brewer(palette = "Dark2",
                     name = "Specimen",
                     breaks = c("4907bc1-d1",
                                "4902bc1-d?",
                                "4684bc1-d?"),
                     labels = c(" Jacksonville (4907)",
                                " Savannah (4902)",
                                " Jacksonville (4684)")) +
  theme_linedraw() +
  # theme(panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank()) +
  labs(x="Calendar Year (C.E.)",
       y=expression(delta^{15}*"N Average Tr-AA")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # Medieval Climate Anomaly, warming
  geom_segment(x=200, y=13.5,xend=585, yend=13.75) +
  ylim(15,25) +
  scale_x_continuous(limits=c(-750, 2005),
                     breaks=seq(-750, 2005, by=550))

traa.time

jack.4686 <- c.jack4
jackplot <- ggplot(jack.4686, aes(distance..mm., y=d15n.vs.air)) +
  # geom_point(aes(color=coral.id), 
  #            size=4.5) +
  # geom_point(shape = 21, color="black", size=4.5) +
  geom_line(aes(color=coral.id), size = 0.85) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2",
                     name = "Specimen",
                     breaks = c("4907bc1-d1",
                                "4902bc1-d?",
                                "4684bc1-d?",
                                "jack-4686-bc-d1-t1"),
                     labels = c(" Jacksonville (4907)",
                                " Savannah (4902)",
                                " Jacksonville (4684)",
                                " Jacksonville (4686)")) +
  labs(x="Distance (mm)",
       y=expression(delta^{15}*"N (\u2030)"))

jackplot

jack.4684 <- c.jack3
jackplot <- ggplot(jack.4684, aes(clam.year.bp, y=d15n.vs.air)) +
  # geom_point(aes(color=coral.id), 
  #            size=4.5) +
  # geom_point(shape = 21, color="black", size=4.5) +
  geom_line(aes(color=coral.id), size = 0.85) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2",
                     name = "Specimen",
                     breaks = c("4907bc1-d1",
                                "4902bc1-d?",
                                "jack-4684-bc-unk",
                                "jack-4686-bc-d1-t1"),
                     labels = c(" Jacksonville (4907)",
                                " Savannah (4902)",
                                " Jacksonville (4684)",
                                " Jacksonville (4686)")) +
  labs(x="Years Before Present (1950 AD)",
       y=expression(delta^{15}*"N (\u2030)"))

jackplot

# ----------------------------------------------------------------------------------
# Coming up with a way to quickly calculate SumV across a lot of different samples
# Each sample will have d15N data from Trophic AA: Ala, Asp, Glu, Ile, Leu, Pro, Val)
# SumV = (1/7)*Sum(Abs(ChiAA))
# (1/7) <-- assuming calculating from seven Trophic AAs (if 6, use 1/6 etc.; some pubs remove Valine)
# ChiAA <- deviation of each Trophic AA --> d15N-AA - Avg d15N-Trophic AA (from all seven)
# Take the absolute value of each ChiAA
# Sum them all and multiply by (1/7)
# ----------------------------------------------------------------------------------

# First, remove all non-Trophic AAs
# gom.tr <- gom %>% select((!!as.name(column)) == 'Ala' | 'Asp' | 'Glu' | 'Ile' | 'Leu' | 'Pro' | 'Val')
gom.tr <- gom %>% select(Ala, Asp, Glu, Ile, Leu, Pro, Val, Avg.Tr)
# Make a column of Ala - Avg Tr and go from there!
gom.tr$Dev.Ala <- abs(gom.tr$Ala - gom.tr$Avg.Tr)
gom.tr$Dev.Asp <- abs(gom.tr$Asp - gom.tr$Avg.Tr)
gom.tr$Dev.Glu <- abs(gom.tr$Glu - gom.tr$Avg.Tr)
gom.tr$Dev.Ile <- abs(gom.tr$Ile - gom.tr$Avg.Tr)
gom.tr$Dev.Leu <- abs(gom.tr$Leu - gom.tr$Avg.Tr)
gom.tr$Dev.Pro <- abs(gom.tr$Pro - gom.tr$Avg.Tr)
gom.tr$Dev.Val <- abs(gom.tr$Val - gom.tr$Avg.Tr)
gom.tr$Sum.Chi <- rowSums(gom.tr[, c(9,11)]) # Make sure you are selecting the correct columns

# Making a loop to do the above
# Still rusty with loops, unfortunately

# Nitrogen
xyplot(d15n.vs.air ~ distance..mm.,
       jack.4686,
       type = "o",
       pch = 20,
       cex = 1.5,
       xlab = "Distance (mm)",
       ylab = expression(delta^{15}*"N"),
       main = "Jacksonville-4686 BC1 (No Age Model)")

# Carbon
xyplot(d13c.vs.vpdb ~ distance..mm.,
       jack.4686,
       type = "o",
       pch = 20,
       cex = 1.5,
       xlab = "Distance (mm)",
       ylab = expression(delta^{13}*"C"),
       main = "Jacksonville-4686 BC1")


xyplot(d13c.ra ~ clam.ad | coral.id,
       bind2,
       type = "l",
       lty = 1,
       lwd = 1,
       col.line = "black")

xyplot()



t <- seq(0, 1, len = 100)
x <- sin(2*pi*t*2.3) + 0.25*rnorm(length(t))
z <- filter(bf, x)
plot(t, x, type="l")
lines(t,z, col = "red")

stet <- c.stet1
stet <- stet %>% select(d15n.vs.air, clam.ad)
stet <- stet %>% select(d15n.vs.air)
stet <- as.numeric(unlist(stet))
stet <- stet[1:376]

library(signal)

bf <- butter(3, c(0.05, 1)) # 0.1 means 10 Hz low pass filter
z <- filter(bf, stet)

xyplot(z)

plot(stet, type = "l")
lines(z, col = "red")

filt <- fftfilt(rep(1, 10)/10, stet)
plot(filt,
     type="l")

y <- filtfilt(bf, stet)
plot(y)
u <- fftfilt(rep(1, 10)/10,stet, n = 385)
View(bf)
k <- Ma()
iodine <- iodine.stet
iodine <- iodine %>% select(Total.Br)
iodine <- as.numeric(unlist(iodine))

z <- filter(filt = bf$b, x = stet, a = bf$a)
z <- filtfilt(filt = bf$b, x = stet, a = bf$a)
plot(z)
plot(u, type ="l")
# Melting the CSIAA data to put in a new column for "trophic" and "source" amino acids; same for NEAA and EAA
# Ultimately so we can make a nice chart

# v <- seus[]
v <- melt(seus, 'Sample.ID')
library(data.table)
v <- data.table(v)
# v[, Group.1 := ifelse(variable == c('Phe','Gly','Ser','Lys','Tyr'), "Source", NA)]
# v[, Group.1 := ifelse(variable %in% c('Phe','Gly','Ser','Lys','Tyr'), "Source", NA)]
# v[, Group.1 := ifelse(variable %in% c('Glu','Asp','Ala','Ile','Leu','Pro','Val'), "Trophic", NA)]
v[, Group.1 := ifelse(variable %in% c('Phe','Gly','Ser','Lys','Tyr'), "Source",
                      ifelse(variable %in% c('Glu','Asp','Ala','Ile','Leu','Pro','Val'), "Trophic", NA))]


t.bulk <- bulk %>% select(d15n.vs.air, d13c.vs.vpdb, clam.ad, coral.id)
t.bulk <- t.bulk %>% filter(coral.id == "stet-4904-bc1-d2" | coral.id == "jack-4684-bc-unk")

xyplot(d13c.vs.vpdb ~ clam.ad,
       t.bulk,
       pch = 21,
       col = c("red", "blue"),
       group=coral.id,
       type="l")


xyplot(value ~ variable1,
       seus.means.melt,
       panel = function(x,  y, ...){
         panel.xyplot(x, y, ...)
         panel.abline(h = 9, lty = 1)
         panel.abline(h = 9.25, lty = 2)
         panel.abline(h = 8.75, lty = 2)
       },
       cex=1.5,
       ylim=c(-20, 40), # ylim and xlim are "first class" parameters and don't need to be in scales=list()
       group = Sample.ID2,
       xlab = NULL,
       ylab = expression({delta}^15*"N (\u2030)"),
       # auto.key=list(columns=2, cex = 0.75),
       auto.key = FALSE,
       par.settings = list(superpose.symbol = list(pch = c(0,2:6), cex = 1.5,
                                                                  col = c('#c7e9b4',
                                                                    '#7fcdbb',
                                                                    '#41b6c4',
                                                                    '#1d91c0',
                                                                    '#225ea8',
                                                                    '#0c2c84'),
                                                   lwd = 3)))


train.means <- aggregate(train.norm[,3:13], list(Sample.ID.2=train.norm$Sample.ID.2), mean, na.rm = TRUE)
train.means2 <- aggregate(t.train[,3:12], list(Sample.ID.2=t.train$Sample.ID.2), mean, na.rm = TRUE)


gg <- droplevels(train.norm2)
str(gg)

biplot(train.pca)
biplot(train.pca2)

summary(train.pca2, loadings = TRUE)

View(USairpollution)
View(cor(USairpollution[,-1]))

library(flexclust)
data(nutrient, package="flexclust")

head(nutrient)
row.names(nutrient) <- tolower(row.names(nutrient))
nutrient.scaled <- scale(nutrient)
d <- dist(nutrient.scaled)
as.matrix(d)[1:4,1:4]

fit.average <- hclust(train.dropped, method="average")
plot(fit.average, hang=-1, cex=0.8, main="Average Linkage Clustering")

rownames(train.dropped) <- train.dropped[,1]
train.scaled <- scale(train.dropped[,-c(1,2)])
binded <- cbind(train.dropped[,c(1,2)],train.scaled)

t.pca <- rda(binded[,-c(1,2)])
biplot(t.pca)

d <- dist(train.dropped)
binded.fit <- hclust(d, method = "ward")

plot(binded.fit, hang=-2, cex=0.5, main="Ward's Linkage Clustering")

training_sample <- sample(c(TRUE, FALSE), nrow(larsenData), replace = T, prob = c(0.6,0.4))
train <- larsenData[training_sample, ]
test <- larsenData[!training_sample, ]

lda.l <- lda(Group.ID2 ~ Ile + Leu + Phe + Thr + Val, train)
lda.l

plot(lda.l, dimen = 2, pch = as.numeric(train$Group.ID2))

lda.train <- predict(lda.l)
train$lda <- lda.train$class
table(train$lda, train$Group.ID2)

lda.test <- predict(lda.l,test)
test$lda <- lda.test$class
table(test$lda, test$Group.ID2)

lda.isidella <- predict(lda.l, isidella)

isidella <- trainData %>% filter(Group.ID2 == "Isidella")
isidella <- isidella %>% select(Group.ID, Group.ID2, Ile, Leu, Phe, Thr, Val)
isidella <- droplevels(isidella)
isidella <- na.omit(isidella)

isidella$lda <- lda.isidella$class


col_names <- colnames(isidella)
for(i in c(1:5)) {
  nth <- isidella[, col_names[i]]
  View(nth)
}



stargazer(stetlm, slm, jlm1, title = "Results", align = FALSE,
          dep.var.labels = c("Stetson Banks", "Savannah Banks", "Jacksonville Lithoherms"))

tempfile <- stargazer(ndata[,-c(1:2)], omit.summary.stat = c("p25", "p75"))

stargazer(seusDataC, summary = FALSE, rownames = FALSE)

t.df <- read.csv("C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff ncsiaa 02-06-2019.csv")
t.df <- t.df %>% select(Coral, sample.no., Ala, Asx, Glx, Gly, Ile, Leu, Phe, Pro, Ser, Val)
stargazer(t.df, summary = FALSE, rownames = FALSE, digits = 1)

radiotable <- radiocarbon
radiotable$Coral[radiotable$Coral == 'jack-4684-bc1'] <- 'Jacksonville-4684'
radiotable$Coral[radiotable$Coral == 'jack-4907-bc1-d1'] <- 'Jacksonville-4907'
radiotable$Coral[radiotable$Coral == 'sav-4902-bc1'] <- 'Savannah Banks-4902'
radiotable$Coral[radiotable$Coral == 'stet-4904-bc1-d5'] <- 'Stetson Banks-4904'

levels(radiotable$Coral) <- c('Jacksonville-4684', 'Jacksonville-4907', 'Savannah Banks-4902','Stetson-old','Stetson Banks-4904')

View(levels(radiotable$Coral))
radiotable$Distance <- radiotable$Distance..um./1000
radiotable$Distance..um. <- NULL
radiotable <- radiotable %>% select(Coral, Sample, Distance, Fraction.modern, error, D14C, d13C, X14C.Age, error.1)
radiotable <- radiotable %>% dplyr::filter(Coral != 'Stetson-old')


stargazer(radiotable, summary = FALSE, rownames = FALSE, digits = 1, digit.separator = "")

seusDataC[,3:ncol(seusDataC)] <- round(seusDataC[,3:ncol(seusDataC)], digits = 1)

frame <- c.stet1

mod <- lm(d15n.vs.air ~ d13c.vs.vpdb, data = frame)
summary(mod)

plot(d15n.vs.air ~ d13c.vs.vpdb, data = frame)
plot(d15n.vs.air ~ distance..mm., data = frame)

plot(d15n.vs.air ~ d13c.vs.vpdb, data = frame[207:nrow(frame),])
mod <- lm(d15n.vs.air ~ d13c.vs.vpdb, data = frame[100:nrow(frame),])
summary(mod)


testmod <- glm(Bulk.N ~ Sum.V, data = csiaa.all)
summary(testmod)

testmod2 <- glm(d15n.vs.air ~ d13c.vs.vpdb, data = frame)
summary(testmod2)

ndata <- 3000
tt <- seq(0, 9, length=ndata)
xt <- sin(pi * tt)
plot(tt)
plot(xt)
extrema(xt)

c <- ssa(co2, L = 120)
plot(c, type = "wcor")
plot(c, type = "vectors")
plot(c$length)

lst <- grouping.auto.wcor(c, group = 1:6, nclust=3)
plot(lst$F2)
r <- reconstruct(c, groups = list(c(1), c(2, 3, 4), c(5, 6)))
plot(pgram(s, groups = list(c(1), c(2,3, 4), c(5, 6)))) 
plot(r)
#' every 500 microns
#' 8 microns per year
#' every 60 years?

f <- TTR::SMA(iodine, n = 500)
library(TTR)
plot(f, type = "l")
plot(iodine, type = "l")
lines(f, col = "blue")

y.smooth <- loess(iodine ~ stetdata$Distance, span = 0.005) # 0.5% smoothing span
plot(y.smooth$fitted, type = "l")
peak <- findpeaks(y.smooth$fitted)
length(peak)
de <- detrend(y.smooth$fitted)
plot(de, type = "l")

plot(stetson$BSE, type = "l")
length(findpeaks(stetson$BSE, minpeakheight = 50))

f <- TTR::SMA(na.omit(stetson$Total.Br), n = 250)
plot(f, type= "l")

y <- stet.table$d15n.ra
y1 <- y+0.3
y2 <- y-0.3
x <- stet.table$year.ad
plot(y ~ x,
       type = "l",
       col = "black",
     ylim = c(6,10))
axis(3,labels=F,tcl=-0.5) 
polygon(c(x,rev(x)),c(y1,rev(y2)),col=alpha("skyblue", alpha = 0.4),
        border = FALSE)

# lines(d15n.vs.air ~ year.ad, data = stet.table,
#       col = alpha("black", alpha = 0.4))
lines(d15n.ra+0.3 ~ year.ad, data = stet.table,
      type = "l",
      lty = 2,
      col = alpha("black", alpha = 0.5))
lines(d15n.ra-0.3 ~ year.ad, data = stet.table,
      type = "l",
      lty = 2,
      col = alpha("black", alpha = 0.5))
lines(d15n.ra ~ min.ad, data = stet.table,
      type = "l",
      lty = 2,
      col = alpha("black", alpha = 0.5))
y1 <- d15n.ra + 0.03
y2 <- d15n.ra - 0.03
polygon(c(year.ad, rev(year.ad)), 
        c(y2, rev(y1)),
        col = "grey30", border = NA)

xyplot(d15n.ra + d15n.vs.air ~ year.ad, data = stet.table,
       type = "l")

v <- stats::spectrum(scale(bromine), method = "pgram")
plot(v, xlim=c(0.2, 0.25))
f <- stats::fft((bromine/length(bromine)), inverse = FALSE)
plot(f)
p <- TSA::periodogram(bromine, plot = TRUE)
plot(scale(bromine), type = "l")

convert.fft <- function(cs, sample.rate=1) {
  cs <- cs / length(cs) # normalize
  
  distance.center <- function(c)signif( Mod(c),        4)
  angle           <- function(c)signif( 180*Arg(c)/pi, 3)
  
  df <- data.frame(cycle    = 0:(length(cs)-1),
                   freq     = 0:(length(cs)-1) * sample.rate / length(cs),
                   strength = sapply(cs, distance.center),
                   delay    = sapply(cs, angle))
  df
}

plot(convert.fft(fft(bromine)))
filt <- spectral::filter.fft(bromine, x = stetdata$Distance)
filt

x <- stetdata$BSE.1
r <- EMD::emd(x, tt = stetdata$Distance, boundar = "periodic")
# Keep 5 components -- you may need more, or less.
y <- apply( r$imf[,5:10], 1, sum ) + mean(r$residue)
plot(x, type="l", col="grey",
     xlim = c(2000,2100))
lines( y, type="l", lwd=2)
n <- length(y)
i <- y[2:(n-1)] > y[1:(n-2)] & y[2:(n-1)] > y[3:n]
points( which(i), y[i], pch=15)

plot(bromine, type = "l")

ndata <- 3000
tt <- seq(0, 9, length=ndata)
xt <- sin(pi * tt)

plot(tt) # tt is just the x-axis
View(tt)
plot(tt, xt)
View(EMD::extrema(xt))

### Generating a signal
ndata <- 3000
par(mfrow=c(1,1), mar=c(1,1,1,1))
tt2 <- seq(0, 9, length=ndata)
xt2 <- sin(pi * tt2) + sin(2* pi * tt2) +
  + sin(6 * pi * tt2) + 0.5 * tt2
plot(tt2, xt2, xlab="", ylab="", type="l",
       + axes=FALSE); box()

### Extracting the first IMF by sifting process
tryimf <- EMD::extractimf(xt2, tt2, check=TRUE)

try <- EMD::emd(x, tt = stetdata$Distance, boundary = "periodic")
y <- apply( try$imf[,5:10], 1, sum ) + mean(try$residue)
plot(x, type="o", col="grey",
     xlim = c(2200,2300))
lines( y, type="l", lwd=2)
n <- length(y)
i <- y[2:(n-1)] > y[1:(n-2)] & y[2:(n-1)] > y[3:n]
points( which(i), y[i], pch=15)
abline(h = 50)

av <- TTR::SMA(x, n = 16)
plot(pracma::detrend(av), x = stetdata$Distance, type = "l")
plot(x ~ stetdata$Distance, 
     xlim = c(2200, 2300),
     type = "l")

x == pracma::detrend(x)

plot(pracma::detrend(x), type = "l")
plot(x)


plot(y ~ x,
     type = "l",
     col = "black",
     ylim = c(6,10))
axis(3,labels=F,tcl=-0.5) 
lines(d15n.ra ~ min.ad, data = stet.table,
      type = "l",
      lty = 2,
      col = alpha("black", alpha = 0.5))
lines(d15n.ra ~ max.ad, data = stet.table,
      type = "l",
      lty = 2,
      col = alpha("black", alpha = 0.5))

##################################
## Jack-4684 bomb spike figure  ##
##################################

jc <- radiocarbon %>% filter(Coral == "jack-4684-bc1")
par(pty = "s")
plot(D14C ~ Distance..um., data = jc,
     type="o",
     xlim = c(0, 2500),
     ylim = c(-75, 150))
lines(rise ~ Distance..um., data = rising, col = "#41b6c4", lwd = 2) # predicted based on lm.rising
lines(desc ~ Distance..um., data = descending, col = "#253494", lwd = 2) # predicted based on lm.descending
# abline(lm.rising, col = "#41b6c4", lwd = 2) # from lm() models below
# abline(lm.desc, col = "#253494", lwd = 2)
rising <- jc %>% filter(Distance..um. < 700 & Distance..um. > 269)
descending <- jc %>% filter(Distance..um. < 300)
lm.rising <- lm(D14C ~ Distance..um., data = rising)
lm.desc <- lm(D14C ~ Distance..um., data = descending)
rise <- predict(lm.rising, newdata = rising)
desc <- predict(lm.desc, newdata = descending)

# Decompose 'co2' series with default parameters
s <- ssa(co2)
# Reconstruct the series, grouping elementary series.
r <- reconstruct(s, groups = list(Trend = c(1, 4), Season1 = c(2,3), Season2 = c(5, 6)))
plot(r)
# 'groups' argument might contain duplicate entries as well
r <- reconstruct(s, groups = list(1, 1:4, 1:6))
plot(r)

p <- findpeaks(n, minpeakheight = 0.025)
plot(n, col = ifelse(n %in% n[p], "red", "black"), type = "o", xlim = c(0,50))

plot(ssa(n))


par(mfrow=c(2,1))
astsa::tsplot(bam3ts, ylab='Sr/Ca 12pt MA', col=4, main='ALV 3808-3', cex.main=1.5,
              xlim=c(1925,2000),
              ylim = c(-1,1))
abline(h = 0, lty = "dashed")
astsa::tsplot(sma.ts, ylab='Sr/Ca 12pt MA', col=4, main='ALV 3808-4', cex.main=1.5,
              xlim=c(1925,2000),
              ylim=c(-1,1))
abline(h = 0, lty = "dashed")

# Radiocarbon figure for Brendan
par(mfrow=c(2,2), pty = "s")
plot(y1 ~ x1,
     main = "Jacksonville BC1 (Dive #4907)",
     ylab = "Calibrated Years BP",
     xlab = "Distance (microns)",
     cex = 1.25)
abline(linearMod1, col="darkblue", lty = "dashed")
abline(h = 1350)
plot(y2 ~ x2,
     main = "Savannah Banks BC1 (Dive #4902)",
     ylab = "Calibrated Years BP",
     xlab = "Distance (microns)",
     cex = 1.25)
abline(linearMod2, col="darkblue", lty = "dashed")
abline(h = 1350)
plot(y4 ~ x4,
     main = "Stetson Banks (Dive #4904)",
     ylab = "Calibrated Years BP",
     xlab = "Distance (microns)",
     cex = 1.25)
abline(linearMod4, col="darkblue", lty = "dashed")
abline(h = 1350)




library(Rssa)
c <- ssa(co2, L = 120)
plot(c)
recon <- reconstruct(c, groups = list(c(1,4), c(2, 3), c(5, 6)))
plot(recon, plot.method = "xyplot", type = "cumsum")


savyrs <- sav %>% filter(sample.no. == 2 | sample.no. == 7 | sample.no. == 15 | sample.no. ==  33
                         | sample.no. == 61 | sample.no. == 66 | sample.no. == 115)
jackyrs <- jack2 %>% filter(sample.no. == 23 | sample.no. == 55 | sample.no. == 56 | sample.no. ==  134)

stetyrs <- df.stet %>% filter(sample.no. == 7 | sample.no. == 8 | sample.no. == 22 | sample.no. ==  25
                              | sample.no. ==  40 | sample.no. ==  46 | sample.no. ==  103 | sample.no. ==  103
                              | sample.no. ==  109 | sample.no. ==  135 | sample.no. ==  170)

plot <- ggplot(means_sd_all, aes(AA, coral.means)) +
  geom_dotplot(binaxis = "y", stackdir = "center", shape = 21) +
  geom_dotplot(aes(AA, euk.means),
               binaxis = "y", stackdir = "center", shape = 22) + 
  geom_dotplot(aes(AA, prok.means),
               binaxis = "y", stackdir = "center", shape = 23)

plot2 <- ggplot(means_sd_all, aes(factor(AA), coral.means)) +
  geom_point() +
  geom_linerange(aes(ymin = coral.means -coral.sd, ymax = coral.means + coral.sd)) +
  geom_point(aes(AA, euk.means), shape = 22) +
  geom_linerange(aes(ymin = euk.means - euk.sd, ymax = euk.means + euk.sd)) +
  geom_point(aes(AA, prok.means), shape = 23) +
  geom_linerange(aes(ymin = prok.means - prok.sd, ymax = prok.means + prok.sd)) +
  scale_x_discrete(limits = names) +
  ylim(-10, 40)

plot2

xyplot(Thr ~ Year.AD, data = ndata,
       group = Sample.ID2)
library(lattice)

par(pty ="s")
plot(d13c.3pt ~ linear.ad, data = df.jack, type = "l", col = "blue", lwd = 2, xlim = c(0, 1500), 
     xlab = x,
     ylab = c)
lines(d13c.3pt ~ linear.ad, data = df.stet, col= "red", lwd = 1.5)
# lines(d13c.3pt ~ linear.ad, data = df.sav, col= "purple", lwd = 2)
legend("bottomright",
       bty = "n",
       col = c("blue", "red"),
       legend = c("Jacksonville (4907)", "Stetson (4904)"),
       lty = 1, cex = 0.75)

par(pty = "s")
plot(d15n ~ distance, df.jack4686, type = "l", col = "purple", lwd = 2.5,
     xlab = "Distance (mm)",
     ylab = n)
legend("topright",
       bty = "n",
       col = c("purple"),
       legend = c("Jacksonville (4686) Live"),
       lty = 1, cex = 0.75)

par(pty ="s")
plot(d15n.3pt ~ linear.ad, data = df.jack, type = "l", ylim = c(6, 10), col = "blue", lwd = 2,
     xlab = x,
     ylab = n)
lines(d15n.3pt ~ linear.ad, data = df.stet, col= "red", lwd = 2)
legend("bottomleft",
       bty = "n",
       col = c("red", "blue"),
       legend = c("Stetson (4904)", "Jacksonville (4907)"),
       lty = 1, cex = 0.75)

par(pty ="s")
plot(d15n.3pt ~ linear.ad, data = df.jack, type = "l", xlim = c(-1010, 1500), ylim = c(7, 11), col = "blue", lwd = 2,
     xlab = x,
     ylab = n)
lines(d15n.3pt ~ linear.ad, data = df.sav, col= "red", lwd = 2)
legend("topleft",
       bty = "n",
       col = c("red", "blue"),
       legend = c("Savannah (4902)", "Jacksonville (4907)"),
       lty = 1, cex = 0.75)

ggplot(data=df.stet, aes(linear.ad, d13c)) +
  theme_bw() +
  geom_point(shape=21) +
  geom_text(data=df.stet, aes(label=df.stet$sample.no.), hjust=0, vjust=-1) +
  geom_line() +
  xlim(1200, 2000) +
  geom_point(data=df.stet[df.stet$sample.no. == "3" | df.stet$sample.no. == "5"
                          | df.stet$sample.no. == "6"
                          | df.stet$sample.no. == "8" | df.stet$sample.no. == "20" | df.stet$sample.no. == "21"
                          | df.stet$sample.no. == "36" | df.stet$sample.no. == "49" 
                          | df.stet$sample.no. == "83"
                          | df.stet$sample.no. == "93",],
             pch=21, fill=NA, size=4, colour="red", stroke=2)


# Doing basic stat calculations for Results section
# Need means, ranges, etc.

## R plot showing Stetson-4904 and Jack-4684 through time

df.bulk %>%
  filter(coral.id == "jack-4684-bc-unk"  | coral.id == "stet-4904-bc1-d2") %>%
  select(coral.id, d13c, d15n, linear.ad) %>%
  na.omit(.) -> tdf

library(lattice)
xyplot(d15n ~ linear.ad,
       group = coral.id,
       xlim = c(1400,2005),
       type = "l",
       lwd = 2,
       ylab = n,
       xlab = x,
       data = tdf)

p <- ggplot(tdf, aes(linear.ad, d15n))
p + geom_line(aes(color = factor(coral.id)), size = 1.5) +
  xlim(1400, 2005) +
  ylab(n) +
  xlab(x) +
  theme_bw() +
  theme(aspect.ratio=1)

##### Calculate N* based on Gruber and Sarmiento, 1997
# N* = (NO3 + NO2) - 16P + 2.9 umol per kg

n_star <- function(nitrate, nitrite, phosphorus) {
  product <- ((nitrate + nitrite) - 16*phosphorus + 2.9) # What are the concentrations?
  return(product)
}

n_star(16,17,100)

n_star <- function(n_plus_n, phosphorus) {
  product <- (n_plus_n - 16*phosphorus + 2.9) # What are the concentrations?
  return(product)
}

library(ggplot2)
library(grid)
library(dplyr)
library(lubridate)

#' Create some data to play with. Two time series with the same timestamp.
df <- data.frame(DateTime = ymd("2010-07-01") + c(0:8760) * hours(2), series1 = rnorm(8761), series2 = rnorm(8761, 100))

#' Create the two plots.
plot1 <- df %>%
  select(DateTime, series1) %>%
  na.omit() %>%
  ggplot() +
  geom_point(aes(x = DateTime, y = series1), size = 0.5, alpha = 0.75) +
  ylab("Red dots / m") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

plot2 <- df %>%
  select(DateTime, series2) %>%
  na.omit() %>%
  ggplot() +
  geom_point(aes(x = DateTime, y = series2), size = 0.5, alpha = 0.75) +
  ylab("Blue drops / L") +
  theme_minimal() +
  theme(axis.title.x = element_blank())

grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))


plot(d15n ~ d13c, df.jack)
plot(d15n ~ d13c, df.sav)
mod <- lm(d15n ~ d13c, df.jack)
mod2 <- lm(d15n ~ d13c, df.sav)

summary(mod)
summary(mod2)

plot(Bulk ~ EAA, data = seus_carbon)
t.model <- lm(Bulk ~ EAA, data = seus_carbon)
summary(t.model)


#' I will now test out the findpeaks() function from the pracma library
#' on a vector of BSE values from the scanned Stetson coral. There
#' should be approximately 1500 peaks if BSE peaks are annual.
#' However, I will use the function on an SSA vector so the peaks are
#' 'cleaner'.
#' 




brom <- read.delim('~/Desktop/bromine_ssavector.txt')
brom$distance.um <- 1:11909

bsevec <- read.delim('~/Desktop/BSEvector.txt', header = FALSE)
bsevec$distance.um <- 1:11910
bsevec$linear.ad <- 2005 - ((bsevec$distance.um)/7.46)

v1 <- r.stet$Distance..um.
v2 <- r.stet$X14C.Age
m1 <- cbind(v1, v2)
m1 <- as.data.frame(m1)
m1 <- m1[-c(4,22,24),]

bsevec %>% 
  filter(distance.um %in% m1$v1) -> new
new$c.age <- m1$v2
new$linear.bp <- 1950 - new$linear.ad

iodinevec <- read.delim('~/Desktop/iodinevector.txt', header = FALSE)
iodinevec$distance.um <- 1:11910
iodinevec$linear.ad <- 2005 - ((iodinevec$distance.um)/7.5)
iodinevec %>% filter(V1 > 0) -> t
plot(V1 ~ linear.ad, data = iodinevec, type = "l")

ggplot(iodinevec, aes(x = linear.ad, color = distance.um)) +
  geom_histogram(binwidth = 5, color = "blue")
  # geom_density(alpha = 0.2, fill = "red")

test <- melt(df.bulk, id = 'coral.id')
ggplot(df.bulk, aes (linear.ad, d15n)) +

View(stet.post)
