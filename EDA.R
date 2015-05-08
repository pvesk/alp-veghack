
### Exploratory data analysis on Alps air photos.

library(GGally)
ggpairs(veg[,c(3:6,8:14)])
pairs(veg[,c(3:6,8:14)])
library(car)
scatterplot.matrix(~ Tr.hit+Sh.hit+Gr.hit+V6+ 
                     Burn.or.clearT1..WD + Burn.or.clearT2..WD + wildfires.between.T1T2 + wildfires.between.T2T3 + Altitude.m + AnnPrecip + MeanAnnTemp|ord
                  , data=veg,
                   main="Three Periods")
by(veg[,c(3:6,8:14)], veg[,15], pairs)
by(veg[,c(3:6)], veg[,1], pairs)
by(veg[,c(3:6)], veg[,1], cor)
cor.tab <-as.data.frame(by(veg[,c(3:6)], veg[,1], cor))

cor.tab <- with(veg,
            by(veg[,c(3:6)], veg[,1],
               cor))
sapply(cor.tab, "[[", 2)
sapply(cor.tab, "[[", 3)
sapply(cor.tab, "[[", 4)
sapply(cor.tab, "[[", 7)
sapply(cor.tab, "[[", 8)
sapply(cor.tab, "[[", 12)


summary(sapply(cor.tab, "[[", 2))
summary(sapply(cor.tab, "[[", 3))
summary(sapply(cor.tab, "[[", 4))
summary(sapply(cor.tab, "[[", 7))
summary(sapply(cor.tab, "[[", 8))
summary(sapply(cor.tab, "[[", 12))


cor.spa <-with(veg,by(veg[,c(3:6,8:14)], veg[,15], cor))
cor.spa[1]

tve correl through time between trees and shrubs
and bw Gr and wetland/other

perfect negtive correlation bw Shrub and Grass
next strongest neg correl bwt Trees and Wetland,
then less strong: 
 other negative correl b/w Shrub and wetland

names(cor.tab)
tmp <- rep(NA, times=23*4*4)
dim(tmp) <- c(4,4,23)
for (i in 1:23) {
  tmp[,,i] <- as.matrix(cor.tab[i])
}
sapply(cor.tab, "[[", 2)
summary(sapply(cor.tab, "[[", 2))
tmp<- as.data.frame(cor(veg[,c(3:6)]))
str(tmp)
summary(tmp)
