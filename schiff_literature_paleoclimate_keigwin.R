#' 
#' This is just for importing and playing around
#' with the paleoclimate data from Keigwin,
#' since it has more variables
#' 
#' 11/15/2019
#' 

# Keiwgin 1996; Bermuda Rise
string1 <- '~/Documents/GitHub/data/paleoclimate_data/bc004a-tab.txt'
string2 <- '~/Documents/GitHub/data/paleoclimate_data/bc004d-tab.txt'
keigwin4a <- read.delim(string1, header = TRUE, comment.char = '#')
keigwin4d <- read.delim(string2, header = TRUE, comment.char = '#')

keigwin4d %>%
  filter(carb. == -999) -> rub
keigwin4d %>%
  filter(d18Og.rub == -999) -> carb
  
plot(d18Og.rub ~ yrBP, rub,
    type = "o",
    bty = "l",
    col = alpha("black", 0.99),
    xlab = "Years BP",
    ylab = "d18O, Core 4D")

plot(carb. ~ yrBP, carb,
     type = "o",
     bty = "l",
     col = alpha("black", 0.99),
     xlab = "Years BP",
     ylab = "% carbonate, Core 4D")

plot(d18Og.rub ~ yrBP, keigwin4a,
     type = "o",
     bty = "l",
     col = alpha("black", 0.99),
     xlab = "Years BP",
     ylab = "d18O, Core 4D")

plot(carb. ~ yrBP, keigwin4a,
     type = "o",
     bty = "l",
     col = alpha("black", 0.99),
     xlab = "Years BP",
     ylab = "% carbonate, Core 4D")

lims <- c(0,3050)
plot(d15n ~ bp, t.stet, type = "o", col = alpha("black", 0.99), xlab = "Years BP", ylab = n, xlim = lims)
plot(d18Og.rub ~ yrBP, rub,
     type = "o",
     bty = "l",
     col = alpha("black", 0.99),
     xlim = lims,
     xlab = "Years BP")
lines(d18Og.rub ~ yrBP, keigwin4a, type = "o", col = alpha("black", 0.3), xlab = "Years BP", ylab = n, xlim = lims)


obj2 <- xyplot(rollmean(d15n, 5, fill = NA, align = "center") ~ bp, data = t.jack, type = "l", lwd = 2, ylab = n, col = "#56B4E9")
obj1 <- xyplot(d18Og.rub ~ yrBP, data = rub, type= "o", lwd = 2, col = "#009E73")
doubleYScale(obj2, obj1, add.ylab2 = TRUE, style1= NULL, style2 = NULL)
