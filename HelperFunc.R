
library(ggplot2)
library(ggfortify)
library(xtable)

fancyHist <- function(x,xlabel, binwid) {
  ggplot(data = as.data.frame(x), 
         aes(x)) + geom_histogram(aes(y = ..density..),binwidth = binwid, alpha = 0.8, fill = "white",
                                  color = "black") +
    labs(x = xlabel, y = "Density") +
   geom_density(fill = "#FF6666", alpha = 0.2)
  # Uncomment if you want a smooth estimate
}
#fancyBoxPlot <- function(x,ylabel) {
#  ggplot(data = as.data.frame(x),aes( y = x,x = "")) + stat_boxplot(geom = "errorbar") +
#    geom_boxplot(outlier.alpha = 0.8) +
#    labs(x = "" , y = ylabel )
#}

tableFactory <- function(x,fileN) {
  print(xtable(x,type = "latex"),file = fileN,booktabs =TRUE)
}