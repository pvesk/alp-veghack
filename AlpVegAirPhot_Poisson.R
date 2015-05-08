#Alpine air photo veg change analysis
library(R2jags)
#source('~/Documents/research/coextinction/R_Stirlings/code/doJAGS.R') 
setwd('~/Documents/research/alpine air photo veg change')

veg <- read.csv('Alpvegcov.csv')   # read in from data folder
veg <- veg[with(veg,order(veg$Site, veg$date)), ]
veg[,3:5] <- round(veg[,3:5]) 
sites <- read.csv('Alpsites.csv')   # read in from data folder
x <- cbind(veg$Tr.hit,veg$Sh.hit,veg$Gr.hit)
n1 <- as.integer(floor(rowSums(x)))
#x <- cbind(veg$Sh.hit,veg$Gr.hit)          
#n2 <- as.integer(floor(rowSums(x)))
rm(x)
veg <- cbind(veg,n1,n2)
yrs <- veg$date-1970
veg <- merge(veg,sites)
veg <- veg[with(veg,order(veg$Site, veg$date)), ]
ysg <- veg$date - veg$last.grazed
ysg <- ifelse(ysg<0,0,ysg)

#hits = as.integer(floor(veg$Tr.hit))
hits = as.integer(floor(veg$Sh.hit))
wf1 <-veg$wildfires.between.T1T2
wf2 <-veg$wildfires.between.T2T3

yswf1 <- ifelse(yrs < -15,0,1)
yswf2 <- ifelse(yrs > 15,1,0)

c.b <- c(50,0,0,  50,0,0,  75,0,0,   0,0,0,   0,75,0,  90,0,0,  20,20,0, 98,0,0,  
         20, 0,0,  90,0,0,   0,60,0,  50,90,0, 98,0,0, 0, 0,0,   0,90,0,  10,0,0,   
         0,0,0,   0,0,0,   0,0,0,  25,25,0,  80,80,0, 0,0,0,  10,0,0)

#c.b <- c(50,0,NA,  50,0,NA,  75,0,NA,   0,0,NA,   0,75,NA,  90,0,NA,  20,20,NA, 98,0,NA,  
 #        20, 0,NA,  90,0,NA,   0,60,NA,  50,90,NA, 98,0,NA, 0, 0,NA,   0,90,NA,  10,0,NA,   
  #       0,0,NA,   0,0,NA,   0,0,NA,  25,25,NA,  80,80,NA, 0,0,NA,  10,0,NA)

#c.b <- ifelse(c.b>0,1,0)
#c.b <- log((c.b+1)/(100-(c.b+1)))
#c.b <-  1-(c.b ) / 100

#define data
veg.data <- list( 
  Nobs=69, 
  Nsites=23,
  Nveg = 3,
  site = as.numeric(as.factor(veg$Site)),
  yrs = yrs,
  Y = veg[,3:5], 
  ysg = ysg,
  wf1=wf1,
  wf2=wf2,
  c.b=c.b,
  yswf1=yswf1,
  yswf2=yswf2
)
# THIS IS THE DATA LIST

veg.inits <- function() list(
  mn.a = c(0,0,0),
  a = structure( c(rep(0, times=69)), .Dim = as.integer(c(23,3))), 
  b =rnorm(3),
  g.ef=rnorm(3),
  f.ef=rnorm(3),
  cb.ef=rnorm(3),
  sd = runif(1,0,5) 
) 
#THIS IS THE INITIALS LIST

########multinomial-poisson transform
Model <- function(){
  for (i in 1:Nobs) {
    for (k in 1:Nveg) {
      Y[i,k] ~ dpois(mu[i,k])
      log(mu[i,k]) <- a[site[i],k]                     # the av count for veg k, with random effect for sites.
                    + b[k] * yrs[i]/10                     # effect of time
                    + g.ef[k] * ysg[i]/10 
                    + f.ef[k] * (wf1[i]* ( yswf1[i] - yswf2[i] ) + wf2[i]* yswf2[i]) /10
                    + cb.ef[k] * c.b[i] #+ log(max(Y[site[i],1]))) 
      
  ## predict for each obs    
      pred[i,k] <- exp( a[site[i],k]                     # the av count for veg k, with random effect for sites.
      + b[k] * yrs[i]/10                     # effect of time
      + g.ef[k] * ysg[i]/10 
      + f.ef[k] * (wf1[i]* ( yswf1[i] - yswf2[i] ) + wf2[i]* yswf2[i]) /10
      + cb.ef[k] * c.b[i] )
    }    
  }
  
  for (j in 1:Nsites) {
    for (k in 1:Nveg) {
      a[j, k] ~ dnorm(mn.a[k], sd^-2)
    }
  }
  
  for (t in 1:60) {   #predictions
    for (k in 1:Nveg) {
      pHit[t,k] <- exp(mn.a[k] + b[k] * (t-30)/10 ) 
      ppHit[t,k] <- pHit[t,k] / 955
      #predict effect of 60% clearing WD in 1955
      pHit.cb[t,k] <- exp(mn.a[k] + b[k] * (t-30)/10 + cb.ef[k] * 0.6 * step(t-15)) 
      ppHit.cb[t,k] <- pHit.cb[t,k] / 955
      #predict for grazing removal in 1964 (median)
      pHit.gr[t,k] <- exp(mn.a[k] + b[k] * (t-30)/10 + g.ef[k] * step(t-24)  /10) 
      ppHit.gr[t,k] <- pHit.gr[t,k] / 955
#predict effect of wildifre in 1975
      pHit.wf[t,k] <- exp(mn.a[k] + b[k] * (t-30)/10 + f.ef[k] * step(t-35) /10) 
      ppHit.wf[t,k] <- pHit.wf[t,k] / 955
    }
  }
 
  for (k in 1:Nveg) {   
      mn.a[k] ~ dnorm(0,1.0E-6)
      #  cover1970[k] <- ilogit(mn.a[k])
      b[k] ~ dnorm(0,1.0E-6)
      g.ef[k] ~ dnorm(0,1.0E-6) # for grazing removal
      f.ef[k] ~ dnorm(0,1.0E-6) # for fe for wildfire
      cb.ef[k] ~ dnorm(0,1.0E-6)  # for clearing & burn
      rate[k] <- exp(b[k])
      rate.g[k] <-exp(b[k]+g.ef[k]) # swa g for f for grazing
  }  
  sd ~ dunif(0,10)
}  # END MODEL

###fit model
AlpVegP <- jags(veg.data, veg.inits, 
                c('mn.a', 'b', 'sd',  'g.ef', 'f.ef', 'rate', 'rate.g', 'cb.ef', 'pred', # 'a',
                  'pHit',  'pHit.cb', 'pHit.gr', 'pHit.wf' ,   
                  'ppHit', 'ppHit.cb', 'ppHit.gr', 'ppHit.wf'), 
                Model, 3, 8000, n.burnin = 3000) 



print(AlpVegP)
par(mfrow = c(4, 2))
traceplot(AlpVegP, mfrow=c(4, 2), varname = c('a', 'mn.a', 'b', 'sd','g.ef','f.ef', 'cb.ef', 'lambda'), ask=FALSE)
attach.jags(AlpVegP)
par(mfrow = c(4,2))
plot(density(mn.a[,1]))
plot(density(mn.a[,2]))
plot(density(mn.a[,3]))
plot(density(a[,1,1]))
plot(density(a[,1,2]))
plot(density(a[,1,3]))
plot(density(b[,1]))
plot(density(b[,2]))
plot(density(b[,3]))
plot(density(g.ef[,1]))
plot(density(g.ef[,2]))
plot(density(g.ef[,3]))
plot(density(f.ef[,1]))
plot(density(f.ef[,2]))
plot(density(f.ef[,3]))
plot(density(cb.ef[,1]))
plot(density(cb.ef[,2]))
plot(density(cb.ef[,3]))

plot(density(sd[, 1]))
plot(density(sd[, 2]))


#AlpVeg.upd <- update(AlpVeg, n.iter=10000)
save(AlpVegP, file="AlpVegP.Rdata")
rm(AlpVeg.upd.upd)
#only want $BUGS.output and $sims.list within that for density plots etc
load('AlpVegP.Rdata')
#traceplot
plot( AlpVegP$BUGSoutput$sims.list$mn.a[,1], type='l')
#denisty plot
plot(density(AlpVegP$BUGSoutput$sims.list$mn.a[,1]))

#plotting for predicted vs obs for trees first
preds <- as.data.frame(matrix(rep(NA,times=69*3*3), nrow=69))
colnames(preds) <- c('wd.med', 'wd.lb','wd.ub', 
                     'sh.med', 'sh.lb','sh.ub',
                     'gr.med', 'gr.lb','gr.ub')
preds <- cbind(preds, veg[,1:5])
preds[,1] <- apply(AlpVegP$BUGSoutput$sims.list$pred[,,1],2,function(x) quantile(x,prob=0.5))
preds[,2] <- apply(AlpVegP$BUGSoutput$sims.list$pred[,,1],2,function(x) quantile(x,prob=0.025))
preds[,3] <- apply(AlpVegP$BUGSoutput$sims.list$pred[,,1],2,function(x) quantile(x,prob=0.975))


plot(preds[,1], preds[,12], xlab = "predicted woodland", ylab = "observed woodland", main = "Woodland")
abline(a = 0, b = 1, col = 2)

#shrubs
preds[,4] <- apply(AlpVegP$BUGSoutput$sims.list$pred[,,2],2,function(x) quantile(x,prob=0.5))
preds[,5] <- apply(AlpVegP$BUGSoutput$sims.list$pred[,,2],2,function(x) quantile(x,prob=0.025))
preds[,6] <- apply(AlpVegP$BUGSoutput$sims.list$pred[,,2],2,function(x) quantile(x,prob=0.975))

plot(preds[,4], preds[,13], xlab = "predicted shrubland", ylab = "observed shrubland", main = "Shrubland")
abline(a = 0, b = 1, col = 2)
#grassland
preds[,7] <- apply(AlpVegP$BUGSoutput$sims.list$pred[,,3],2,function(x) quantile(x,prob=0.5))
preds[,8] <- apply(AlpVegP$BUGSoutput$sims.list$pred[,,3],2,function(x) quantile(x,prob=0.025))
preds[,9] <- apply(AlpVegP$BUGSoutput$sims.list$pred[,,3],2,function(x) quantile(x,prob=0.975))

plot(preds[,7], preds[,14], xlab = "predicted grassland", ylab = "observed grassland",main = "Grassland")
abline(a = 0, b = 1, col = 2)


### now for continuous predictions 
### set up the dataframe to collect the predictions
p.hits <- as.data.frame(matrix(rep(NA,times=60*(36+1)), nrow=60))
colnames(p.hits) <- c( "yr", 'wd.med', 'wd.lb','wd.ub', 
                     'sh.med', 'sh.lb','sh.ub',
                     'gr.med', 'gr.lb','gr.ub',
                     
                     'wd.med.cb', 'wd.lb.cb','wd.ub.cb', 
                     'sh.med.cb', 'sh.lb.cb','sh.ub.cb',
                     'gr.med.cb', 'gr.lb.cb','gr.ub.cb',
                     
                     'wd.med.gr', 'wd.lb.gr','wd.ub.gr', 
                     'sh.med.gr', 'sh.lb.gr','sh.ub.gr',
                     'gr.med.gr', 'gr.lb.gr','gr.ub.gr',
                     
                     'wd.med.wf', 'wd.lb.wf','wd.ub.wf', 
                     'sh.med.wf', 'sh.lb.wf','sh.ub.wf',
                     'gr.med.wf', 'gr.lb.wf','gr.ub.wf')

p.hits[,1] <- 1940:1999
#trees
p.hits[,2] <- apply(AlpVegP$BUGSoutput$sims.list$pHit[,,1],2,function(x) quantile(x,prob=0.5))
p.hits[,3] <- apply(AlpVegP$BUGSoutput$sims.list$pHit[,,1],2,function(x) quantile(x,prob=0.025))
p.hits[,4] <- apply(AlpVegP$BUGSoutput$sims.list$pHit[,,1],2,function(x) quantile(x,prob=0.975))

#shrubs
p.hits[,5] <- apply(AlpVegP$BUGSoutput$sims.list$pHit[,,2],2,function(x) quantile(x,prob=0.5))
p.hits[,6] <- apply(AlpVegP$BUGSoutput$sims.list$pHit[,,2],2,function(x) quantile(x,prob=0.025))
p.hits[,7] <- apply(AlpVegP$BUGSoutput$sims.list$pHit[,,2],2,function(x) quantile(x,prob=0.975))

#grassland
p.hits[,8] <- apply(AlpVegP$BUGSoutput$sims.list$pHit[,,3],2,function(x) quantile(x,prob=0.5))
p.hits[,9] <- apply(AlpVegP$BUGSoutput$sims.list$pHit[,,3],2,function(x) quantile(x,prob=0.025))
p.hits[,10] <- apply(AlpVegP$BUGSoutput$sims.list$pHit[,,3],2,function(x) quantile(x,prob=0.975))

## now for the cleared and burnt
#trees
p.hits[,11] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.cb[,,1],2,function(x) quantile(x,prob=0.5))
p.hits[,12] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.cb[,,1],2,function(x) quantile(x,prob=0.025))
p.hits[,13] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.cb[,,1],2,function(x) quantile(x,prob=0.975))

#shrubs
p.hits[,14] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.cb[,,2],2,function(x) quantile(x,prob=0.5))
p.hits[,15] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.cb[,,2],2,function(x) quantile(x,prob=0.025))
p.hits[,16] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.cb[,,2],2,function(x) quantile(x,prob=0.975))

#grassland
p.hits[,17] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.cb[,,3],2,function(x) quantile(x,prob=0.5))
p.hits[,18] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.cb[,,3],2,function(x) quantile(x,prob=0.025))
p.hits[,19] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.cb[,,3],2,function(x) quantile(x,prob=0.975))

##now for the grazing removal:
#trees
p.hits[,20] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.gr[,,1],2,function(x) quantile(x,prob=0.5))
p.hits[,21] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.gr[,,1],2,function(x) quantile(x,prob=0.025))
p.hits[,22] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.gr[,,1],2,function(x) quantile(x,prob=0.975))

#shrubs
p.hits[,23] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.gr[,,2],2,function(x) quantile(x,prob=0.5))
p.hits[,24] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.gr[,,2],2,function(x) quantile(x,prob=0.025))
p.hits[,25] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.gr[,,2],2,function(x) quantile(x,prob=0.975))

#grassland
p.hits[,26] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.gr[,,3],2,function(x) quantile(x,prob=0.5))
p.hits[,27] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.gr[,,3],2,function(x) quantile(x,prob=0.025))
p.hits[,28] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.gr[,,3],2,function(x) quantile(x,prob=0.975))

### now for the wildfire:
#trees
p.hits[,29] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.wf[,,1],2,function(x) quantile(x,prob=0.5))
p.hits[,30] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.wf[,,1],2,function(x) quantile(x,prob=0.025))
p.hits[,31] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.wf[,,1],2,function(x) quantile(x,prob=0.975))

#shrubs
p.hits[,32] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.wf[,,2],2,function(x) quantile(x,prob=0.5))
p.hits[,33] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.wf[,,2],2,function(x) quantile(x,prob=0.025))
p.hits[,34] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.wf[,,2],2,function(x) quantile(x,prob=0.975))

#grassland
p.hits[,35] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.wf[,,3],2,function(x) quantile(x,prob=0.5))
p.hits[,36] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.wf[,,3],2,function(x) quantile(x,prob=0.025))
p.hits[,37] <- apply(AlpVegP$BUGSoutput$sims.list$pHit.wf[,,3],2,function(x) quantile(x,prob=0.975))


## PLOT: Here for woodland amount
obs <- veg$Tr.hit

obs.pts <- cbind(veg,obs)
obs.pts <- subset(obs.pts, select = c(Site, date, obs))
obs.pts <- obs.pts[with(obs.pts,order(obs.pts$Site, obs.pts$date)), ]
ord <- rep(1:3,times=23)
obs.pts <- cbind(obs.pts,ord)
obs.pts.w <- reshape(obs.pts, timevar="ord", idvar="Site", direction="wide")

dow <- obs.pts.w[,c(3,5,7)]
dyw <- obs.pts.w[,c(2,4,6)]

plot.new()
plot.window(xlim=c(-5,65),ylim=c(0,800))
box()
axis(1, at=(0:6)*10,labels=((0:6)*10)+1940)
axis(2, las=1)
title(main = "Woodland", xlab = 'Year', ylab = 'Points')


#observed
for (i in 1:23) { 
  points(x=dyw[i,]-1940, y=dow[i,], type='l', col = 'grey80', lwd=2) 
  points(x=dyw[i,1]-1940, y=dow[i,1], cex = veg.w[i,8]/30, pch = 19)
  points(x=dyw[i,2]-1940, y=dow[i,2], cex = veg.w[i,9]/30, pch = 19)
}

plot.new()
plot.window(xlim=c(-5,65),ylim=c(0,800))
box()
axis(1, at=(0:6)*10,labels=((0:6)*10)+1940)
axis(2, las=1)
title(main = "Woodland prediction", xlab = 'Year', ylab = 'Points')

preds.a <- preds[,1]
dim(preds.a) = c(3,23)
preds.a <- t(preds.a)

#predicted 
for (i in 1:23) { 
  points(x=dyw[i,]-1940, y=preds.a[i,1:3], type='l', col = 'grey80', lwd=2)   
  points(x=dyw[i,1]-1940, y=preds.a[i,1], cex = veg.w[i,8]/30, pch = 19)
  points(x=dyw[i,2]-1940, y=preds.a[i,2], cex = veg.w[i,9]/30, pch = 19)
}

#now for mean predictions over observations,
plot.new()
plot.window(xlim=c(-5,65),ylim=c(0,800))
box()
axis(1, at=(0:6)*10,labels=((0:6)*10)+1940)
axis(2, las=1)
title(main = "Woodland", xlab = 'Year', ylab = 'Points')
#observed
for (i in 1:23) { 
  points(x=dyw[i,]-1940, y=dow[i,], type='l', col = 'grey80', lwd=2) 
  points(x=dyw[i,1]-1940, y=dow[i,1], cex = veg.w[i,8]/30, pch = 19)
  points(x=dyw[i,2]-1940, y=dow[i,2], cex = veg.w[i,9]/30, pch = 19)
}

## undisturbed
points(1:60,p.hits[,2], type='l', lwd=2)
points(1:60,p.hits[,3], type='l',lty=2, lwd=2)
points(1:60,p.hits[,4], type='l',lty=2, lwd=2)

##clear & burn
points(1:60,p.hits[,11], type='l', col = 'blue', lwd=2)  
points(1:60,p.hits[,12], type='l',lty=2, col = 'blue', lwd=2)
points(1:60,p.hits[,13], type='l',lty=2, col = 'blue', lwd=2)

###grazing removal
points(1:60,p.hits[,20], type='l', col = 'green', lwd=2)  
points(1:60,p.hits[,21], type='l',lty=2, col = 'green', lwd=2)
points(1:60,p.hits[,22], type='l',lty=2, col = 'green', lwd=2)

##wildfire
points(1:60,p.hits[,29], type='l', col = 'red', lwd=2)
points(1:60,p.hits[,30], type='l',lty=2, col = 'red', lwd=2)
points(1:60,p.hits[,31], type='l',lty=2, col = 'red', lwd=2)


##### for shrubs
obs <- veg$Sh.hit
obs.pts <- cbind(veg,obs)
obs.pts <- subset(obs.pts, select = c(Site, date, obs))
obs.pts <- obs.pts[with(obs.pts,order(obs.pts$Site, obs.pts$date)), ]
ord <- rep(1:3,times=23)
obs.pts <- cbind(obs.pts,ord)
obs.pts.w <- reshape(obs.pts, timevar="ord", idvar="Site", direction="wide")

dow <- obs.pts.w[,c(3,5,7)]
dyw <- obs.pts.w[,c(2,4,6)]

plot.new()
plot.window(xlim=c(-5,65),ylim=c(0,1000))
box()
axis(1, at=(0:6)*10,labels=((0:6)*10)+1940)
axis(2, las=1)
title(main = "Shrubland", xlab = 'Year', ylab = 'Points')

for (i in 1:23) { 
  points(x=dyw[i,]-1940, y=dow[i,], type='l', col = 'grey80', lwd=2) 
}

#Undisturbed
points(1:60,p.hits[,5], type='l', lwd=2)
points(1:60,p.hits[,6], type='l',lty=2, lwd=2)
points(1:60,p.hits[,7], type='l',lty=2, lwd=2)


##clear & burn
points(1:60,p.hits[,14], type='l', col = 'blue', lwd=2)  
points(1:60,p.hits[,15], type='l',lty=2, col = 'blue', lwd=2)
points(1:60,p.hits[,16], type='l',lty=2, col = 'blue', lwd=2)

###grazing removal
points(1:60,p.hits[,23], type='l', col = 'green', lwd=2)  
points(1:60,p.hits[,24], type='l',lty=2, col = 'green', lwd=2)
points(1:60,p.hits[,25], type='l',lty=2, col = 'green', lwd=2)

##wildfire
points(1:60,p.hits[,32], type='l', col = 'red', lwd=2)
points(1:60,p.hits[,33], type='l',lty=2, col = 'red', lwd=2)
points(1:60,p.hits[,34], type='l',lty=2, col = 'red', lwd=2)


############# for grassland

obs <- veg$Gr.hit
obs.pts <- cbind(veg,obs)
obs.pts <- subset(obs.pts, select = c(Site, date, obs))
obs.pts <- obs.pts[with(obs.pts,order(obs.pts$Site, obs.pts$date)), ]
ord <- rep(1:3,times=23)
obs.pts <- cbind(obs.pts,ord)
obs.pts.w <- reshape(obs.pts, timevar="ord", idvar="Site", direction="wide")

dow <- obs.pts.w[,c(3,5,7)]
dyw <- obs.pts.w[,c(2,4,6)]

plot.new()

plot.window(xlim=c(-5,65),ylim=c(0,1400))
box()
axis(1, at=(0:6)*10,labels=((0:6)*10)+1940)
axis(2, las=1)
title(main = "Grassland", xlab = 'Year', ylab = 'Points')
for (i in 1:23) { 
  points(x=dyw[i,]-1940, y=dow[i,], type='l', col = 'grey80', lwd=2) 
}

#Undisturbed
points(1:60,p.hits[,8], type='l', lwd=2)
points(1:60,p.hits[,9], type='l',lty=2, lwd=2)
points(1:60,p.hits[,10], type='l',lty=2, lwd=2)

##clear & burn
points(1:60,p.hits[,17], type='l', col = 'blue', lwd=2)  
points(1:60,p.hits[,18], type='l',lty=2, col = 'blue', lwd=2)
points(1:60,p.hits[,19], type='l',lty=2, col = 'blue', lwd=2)

###grazing removal
points(1:60,p.hits[,26], type='l', col = 'green', lwd=2)  
points(1:60,p.hits[,27], type='l',lty=2, col = 'green', lwd=2)
points(1:60,p.hits[,28], type='l',lty=2, col = 'green', lwd=2)

##wildfire
points(1:60,p.hits[,35], type='l', col = 'red', lwd=2)
points(1:60,p.hits[,36], type='l',lty=2, col = 'red', lwd=2)
points(1:60,p.hits[,37], type='l',lty=2, col = 'red', lwd=2)



#######################################
#### as proportions  #################
pp.hits <- as.data.frame(matrix(rep(NA,times=60*(36+1)), nrow=60))
colnames(pp.hits) <- c( "yr", 'wd.med', 'wd.lb','wd.ub', 
                       'sh.med', 'sh.lb','sh.ub',
                       'gr.med', 'gr.lb','gr.ub',
                       
                       'wd.med.cb', 'wd.lb.cb','wd.ub.cb', 
                       'sh.med.cb', 'sh.lb.cb','sh.ub.cb',
                       'gr.med.cb', 'gr.lb.cb','gr.ub.cb',
                       
                       'wd.med.gr', 'wd.lb.gr','wd.ub.gr', 
                       'sh.med.gr', 'sh.lb.gr','sh.ub.gr',
                       'gr.med.gr', 'gr.lb.gr','gr.ub.gr',
                       
                       'wd.med.wf', 'wd.lb.wf','wd.ub.wf', 
                       'sh.med.wf', 'sh.lb.wf','sh.ub.wf',
                       'gr.med.wf', 'gr.lb.wf','gr.ub.wf')

pp.hits[,1] <- 1940:1999


#trees
pp.hits[,2] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit[,,1],2,function(x) quantile(x,prob=0.5))
pp.hits[,3] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit[,,1],2,function(x) quantile(x,prob=0.025))
pp.hits[,4] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit[,,1],2,function(x) quantile(x,prob=0.975))

#shrubs
pp.hits[,5] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit[,,2],2,function(x) quantile(x,prob=0.5))
pp.hits[,6] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit[,,2],2,function(x) quantile(x,prob=0.025))
pp.hits[,7] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit[,,2],2,function(x) quantile(x,prob=0.975))

#grassland
pp.hits[,8] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit[,,3],2,function(x) quantile(x,prob=0.5))
pp.hits[,9] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit[,,3],2,function(x) quantile(x,prob=0.025))
pp.hits[,10] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit[,,3],2,function(x) quantile(x,prob=0.975))


## now for the cleared and burnt
#trees
pp.hits[,11] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.cb[,,1],2,function(x) quantile(x,prob=0.5))
pp.hits[,12] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.cb[,,1],2,function(x) quantile(x,prob=0.025))
pp.hits[,13] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.cb[,,1],2,function(x) quantile(x,prob=0.975))

#shrubs
pp.hits[,14] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.cb[,,2],2,function(x) quantile(x,prob=0.5))
pp.hits[,15] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.cb[,,2],2,function(x) quantile(x,prob=0.025))
pp.hits[,16] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.cb[,,2],2,function(x) quantile(x,prob=0.975))

#grassland
pp.hits[,17] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.cb[,,3],2,function(x) quantile(x,prob=0.5))
pp.hits[,18] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.cb[,,3],2,function(x) quantile(x,prob=0.025))
pp.hits[,19] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.cb[,,3],2,function(x) quantile(x,prob=0.975))

##now for the grazing removal:
#trees
pp.hits[,20] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.gr[,,1],2,function(x) quantile(x,prob=0.5))
pp.hits[,21] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.gr[,,1],2,function(x) quantile(x,prob=0.025))
pp.hits[,22] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.gr[,,1],2,function(x) quantile(x,prob=0.975))

#shrubs
pp.hits[,23] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.gr[,,2],2,function(x) quantile(x,prob=0.5))
pp.hits[,24] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.gr[,,2],2,function(x) quantile(x,prob=0.025))
pp.hits[,25] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.gr[,,2],2,function(x) quantile(x,prob=0.975))

#grassland
pp.hits[,26] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.gr[,,3],2,function(x) quantile(x,prob=0.5))
pp.hits[,27] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.gr[,,3],2,function(x) quantile(x,prob=0.025))
pp.hits[,28] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.gr[,,3],2,function(x) quantile(x,prob=0.975))

### now for the wildfire:
#trees
pp.hits[,29] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.wf[,,1],2,function(x) quantile(x,prob=0.5))
pp.hits[,30] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.wf[,,1],2,function(x) quantile(x,prob=0.025))
pp.hits[,31] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.wf[,,1],2,function(x) quantile(x,prob=0.975))

#shrubs
pp.hits[,32] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.wf[,,2],2,function(x) quantile(x,prob=0.5))
pp.hits[,33] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.wf[,,2],2,function(x) quantile(x,prob=0.025))
pp.hits[,34] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.wf[,,2],2,function(x) quantile(x,prob=0.975))

#grassland
pp.hits[,35] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.wf[,,3],2,function(x) quantile(x,prob=0.5))
pp.hits[,36] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.wf[,,3],2,function(x) quantile(x,prob=0.025))
pp.hits[,37] <- apply(AlpVegP$BUGSoutput$sims.list$ppHit.wf[,,3],2,function(x) quantile(x,prob=0.975))

#######################################
#### as cumulative proportions  #################
###########################################
cum.hits <- pp.hits

#### cumulative  for trees  
#grassland at the bottom
# as was
#shrubs
cum.hits[,5] <- cum.hits[,5] + cum.hits[,8]
cum.hits[,6] <- cum.hits[,6] + cum.hits[,9]
cum.hits[,7] <- cum.hits[,7] + cum.hits[,10]
#woodland
cum.hits[,2] <- cum.hits[,5] + cum.hits[,2]
cum.hits[,3] <- cum.hits[,6] + cum.hits[,3]
cum.hits[,4] <- cum.hits[,7] + cum.hits[,14]

####################

obs <- veg$Gr.hit
obs.pts <- cbind(veg,obs)
obs.pts <- subset(obs.pts, select = c(Site, date, obs))
obs.pts <- obs.pts[with(obs.pts,order(obs.pts$Site, obs.pts$date)), ]
ord <- rep(1:3,times=23)
obs.pts <- cbind(obs.pts,ord)
obs.pts.w <- reshape(obs.pts, timevar="ord", idvar="Site", direction="wide")

dow.g <- obs.pts.w[,c(3,5,7)]

dyw <- obs.pts.w[,c(2,4,6)]
obs <- veg$Sh.hit
obs.pts <- cbind(veg,obs)
obs.pts <- subset(obs.pts, select = c(Site, date, obs))
obs.pts <- obs.pts[with(obs.pts,order(obs.pts$Site, obs.pts$date)), ]
ord <- rep(1:3,times=23)
obs.pts <- cbind(obs.pts,ord)
obs.pts.w <- reshape(obs.pts, timevar="ord", idvar="Site", direction="wide")

dow.t <- obs.pts.w[,c(3,5,7)]

obs <- veg$Sh.hit
obs.pts <- cbind(veg,obs)
obs.pts <- subset(obs.pts, select = c(Site, date, obs))
obs.pts <- obs.pts[with(obs.pts,order(obs.pts$Site, obs.pts$date)), ]
ord <- rep(1:3,times=23)
obs.pts <- cbind(obs.pts,ord)
obs.pts.w <- reshape(obs.pts, timevar="ord", idvar="Site", direction="wide")

dow.s <- obs.pts.w[,c(3,5,7)]

cum.obs <- dyw 
cum.obs <- cbind(cum.obs, dow.g, dow.s, dow.t)
####grassland cover, as was
####shrubland + grassland cover
cum.obs[,7] <- cum.obs[,4] + cum.obs[,7]  
cum.obs[,8] <- cum.obs[,5] + cum.obs[,8]  
cum.obs[,9] <- cum.obs[,6] + cum.obs[,9]  
####shrubland + grassland + woodland cover
cum.obs[,10] <- cum.obs[,10] + cum.obs[,7]  
cum.obs[,11] <- cum.obs[,11] + cum.obs[,8]  
cum.obs[,12] <- cum.obs[,12] + cum.obs[,9]  
## grassland as proportion
cum.obs[,4]  <- cum.obs[,4] / cum.obs[,10]                  
cum.obs[,5]  <- cum.obs[,5] / cum.obs[,11] 
cum.obs[,6]  <- cum.obs[,6] / cum.obs[,12] 
## grassland & shrubaland as proportion
cum.obs[,7]  <- cum.obs[,7] / cum.obs[,10]                  
cum.obs[,8]  <- cum.obs[,8] / cum.obs[,11] 
cum.obs[,9]  <- cum.obs[,9] / cum.obs[,12] 

#### build the plots for proportional cover

#windows()
plot.new()
plot.window(xlim=c(-5,65),ylim=c(0,1.0))
box()
axis(1, at=(0:6)*10,labels=((0:6)*10)+1940)
axis(2, las=1)

for (i in 1:23) { 
  points(x=cum.obs[i,1:3]-1940, y=cum.obs[i,4:6], type='l', col = 'grey80', lwd=2) 
  points(x=cum.obs[i,1:3]-1940, y=cum.obs[i,7:9], type='l', col = 'black', lwd=2) 
}
#grassland is grey
#shrub+grass is Black, remainder iw woodland


points(1:60,p.hits[,2], type='l', lwd=2)
points(1:60,p.hits[,3], type='l',lty=2, lwd=2)
points(1:60,p.hits[,4], type='l',lty=2, lwd=2)

points(1:60,mn.pg, type='l', col = 'red', lwd=2)
points(1:60,lb.pg, type='l',lty=2, col = 'red', lwd=2)
points(1:60,ub.pg, type='l',lty=2, col = 'red', lwd=2)

points(1:60,mn.pcb, type='l', col = 'blue', lwd=2)
points(1:60,lb.pcb, type='l',lty=2, col = 'blue',lwd=2)
points(1:60,ub.pcb, type='l',lty=2, col = 'blue',lwd=2)

mtext('Years',1, line=3)
mtext('Proportion',2, line=3)

dev.off()


