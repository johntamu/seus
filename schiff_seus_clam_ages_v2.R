#'-----------------------------------
#' Age-depth models using clam in R
#' Also includes Bacon age-depth modeling
#'-----------------------------------
#'
library(clam)
library(rbacon)

core.dir <- "~/Documents/GitHub/data/clam/Cores"

# Stetson-4904 BC1 
clam(core="stetson", type = 4, smooth = 0.8, prob = 0.95, its = 1000,
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, thickness = 0.1, est = 2)

Bacon(core = "stetson", thick = 0.1, depth.unit = 'mm')

t.stet %>%
  select(distance, d15n, d13c) -> proxies
write.csv(proxies, '~/Documents/GitHub/data/clam/Cores/stetson/stetson_proxies.csv')
stetdepths <- read.delim('~/Documents/GitHub/data/clam/Cores/stetson/stetson_smooth_spline_ages.txt')

# Jacksonville-4907 BC1
clam(core="jack4907", type = 4, smooth = 0.6, prob = 0.95, its = 1000,
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, proxies = TRUE, thickness = 0.1, youngest = c(1), est = 1)
jackdepths <- read.delim('~/Documents/GitHub/data/clam/Cores/jack4907/jack4907_smooth_spline_ages.txt')
# Savannah-4902 BC1
clam(core="sav4902", type = 4, smooth = 0.6, prob = 0.95, its = 1000,
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, proxies = TRUE)
t.sav %>% 
  select(distance) -> depths
write.csv(depths, '~/Documents/GitHub/data/clam/Cores/sav4902/sav4902_depths.csv')
