#' ------------------------
#' Analysis of nutrient data 
#' from Nancy Prouty/USGS
#' 08/27/2019
#' ------------------------

# Libraries
library(lattice)
library(latticeExtra)
library(ggplot2)

# Read data
nitrate_d15n <- read.csv('~/Documents/GitHub/data/nutrients_d15n.csv')
head(nitrate_d15n)

# Plot

# Lattice
xyplot(depth ~ d15N,
       group = site,
       data = nitrate_d15n,
       xlab = n, # From figures script file
       ylab = 'Depth (m)',
       ylim = rev(c(0,2500)),
       xlim = c(0,8))
       
# ggplot (easier)
colors <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99')
p <- ggplot(nitrate_d15n, aes(x=d15N, y=depth))
p + geom_point(aes(fill= factor(site)), shape = 23, size = 3) +
  scale_y_continuous(trans = "reverse") +
  scale_x_continuous(position = "top", limits = c(0,8)) +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(panel.border = element_rect(color = 'black'),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.25, color = "black"),
        axis.line = element_blank()) +
  ylab('Depth (m)') +
  xlab(n)
