library(lattice)
library(latticeExtra)
library(dplyr)

View(df.bulk)

df.bulk %>%
  filter(coral.id == "jack-4684-bc-unk" | coral.id == "jack-4907-bc1-d1" |
           coral.id == "jack-4907-bc1-d3" | coral.id == "sav-4902-bc1-unk" |
           coral.id == "stet-4904-bc1-d2") -> bulk2

xyplot() # plan is to create an xyplot of each of these corals to see what the "overall" record looks like

View(corals)
plot(TP ~ bp, corals)
max(corals$TP)
min(corals$TP)
mean(corals$TP)
sd(corals$TP)

corals %>%
  filter(Sample.ID2 == "Jacksonville-4684") -> t.coral
corals %>%
  filter(Sample.ID2 == "Jacksonville-4907" | Sample.ID2 == "Savannah Banks-4902") -> t.coral2
corals %>%
  filter(Sample.ID2 == "Savannah Banks-4902") -> t.coral3

mean(t.coral$TP)
sd(t.coral$TP)
mean(t.coral2$TP)
sd(t.coral2$TP)

max(corals$Sum.V)
min(corals$Sum.V)
mean(corals$Sum.V)
sd(corals$Sum.V)

plot(Sum.V ~ bp, corals)
corals <- droplevels(corals)
xyplot(Sum.V ~ bp, data = corals,
       group = Sample.ID2, auto.key = TRUE)
