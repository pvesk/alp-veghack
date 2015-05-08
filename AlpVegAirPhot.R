#Alpine air photo veg change analysis
library(R2jags)
#source('~/Documents/research/coextinction/R_Stirlings/code/doJAGS.R') 
setwd('~/Documents/research/alpine air photo veg change')

veg <- read.csv('Alpvegcov.csv')   # read in from data folder
veg <- veg[with(veg,order(veg$Site, veg$date)), ]
sites <- read.csv('Alpsites.csv')   # read in from data folder
x <- cbind(veg$Tr.hit,veg$Sh.hit,veg$Gr.hit)
n1 <- as.integer(floor(rowSums(x)))
x <- cbind(veg$Sh.hit,veg$Gr.hit)          
n2 <- as.integer(floor(rowSums(x)))
rm(x)
veg <- cbind(veg,n1,n2)
yrs <- veg$date-1970
veg <- merge(veg,sites)
veg <- veg[with(veg,order(veg$Site, veg$date)), ]
ysg <- veg$date-veg$last.grazed
ysg <- ifelse(ysg<0,0,ysg)
#hits = as.integer(floor(veg$Tr.hit))
hits = as.integer(floor(veg$Sh.hit))
wf1 <-veg$wildfires.between.T1T2
wf2 <-veg$wildfires.between.T2T3
#obs <- veg$Tr.hit/veg$n1
obs <- veg$Sh.hit/veg$n2
obs.pts <- cbind(veg,obs)
obs.pts <- subset(obs.pts, select = c(Site, date, obs))
obs.pts <- obs.pts[with(obs.pts,order(obs.pts$Site, obs.pts$date)), ]
ord <- rep(1:3,times=23)
obs.pts <- cbind(obs.pts,ord)
obs.pts.w <- reshape(obs.pts, timevar="ord", idvar="Site", direction="wide")

dow <- obs.pts.w[,c(3,5,7)]
dyw <- obs.pts.w[,c(2,4,6)]
yswf1 <- ifelse(yrs < -15,0,1)
yswf2 <- ifelse(yrs > 15,1,0)

c.b <- c(50,0,0,  50,0,0,  75,0,0,   0,0,0,   0,75,0,  90,0,0,  20,20,0, 98,0,0,  
        20, 0,0,  90,0,0,   0,60,0,  50,90,0, 98,0,0, 0, 0,0,   0,90,0,  10,0,0,   
        0,0,0,   0,0,0,   0,0,0,  25,25,0,  80,80,0, 0,0,0,  10,0,0)
#c.b <- ifelse(c.b>0,1,0)
#c.b <- log((c.b+1)/(100-(c.b+1)))

                                    
#define data
veg.data <- list(
  Nobs=69, 
  Nsites=23,                 
  site = as.numeric(as.factor(veg$Site)),
  yrs = yrs,
  hits= hits,
  Npts = n1,
  ysg =ysg,
  wf1=wf1,
  wf2=wf2,
  c.b=c.b,
  yswf1=yswf1,
  yswf2=yswf2
  )
# THIS IS THE DATA LIST

veg.inits <- function() list(
  int = rnorm(1),
  b =rnorm(1),
  ge=rnorm(1),
  fe=rnorm(1),
  cb.ef=rnorm(1),
  sd = runif(2,0,5),
  Y = rnorm(69)) 
#THIS IS THE INITIALS LIST

#model
Model <- function(){
  for (i in 1:Nobs) {
    hits[i] ~ dbin(p[i],Npts[i])
    logit(p[i]) <- Y[i]
    Y[i] ~ dnorm(mu[i], prec[1])
 #   mu[i] <- lambda[i] + cb.ef * c.b[i]
  #  lambda[i] <- 
    mu[i] <- a[site[i]]  
    + b * yrs[i] 
    + g.ef * ysg[i] 
    + f.ef * (wf1[i]* ( yswf1[i] - yswf2[i] ) + wf2[i]* yswf2[i]) 
    #short-term effect  #b * yrs[i]   
    # fe * (wf1*step(yrs[i]+15) + wf2*step(yrs[i]-15)) #persistent effect
        
    l.pred[i] <- a[site[i]] +  b * yrs[i] + g.ef * ysg[i]  ## 
    #l.pred.g[i] <- a[site[i]] +  b * yrs[i] + g.ef*ysg[i]
    # + fe * (wf1[i]* ( yswf1[i] - yswf2[i] ) + wf2[i]* yswf2[i]) 
    pred[i] <- ilogit(l.pred[i])
    dv[i] <- c.b[i] + ysg[i] + wf1[i] + wf2[i] +yswf2[1]
  }
  for (j in 1:Nsites) {
    a[j] ~ dnorm(int, prec[2])
  }
  for (k in 1:60) {
    l.p[k] <- int + b * (k-30)
    p.veg[k] <- ilogit(l.p[k])
    l.pg[k] <- int + (b + g.ef * step(k-36)) * (k-30) #predict for grazing removal in 1964 (median)
  #  l.pcb[k] <- int + b * (k-30) + cb.ef *step(k-20) #predict effect of clearing in 1960
    logit(pg[k]) <- l.pg[k]
   # logit(pcb[k]) <- l.pcb[k]
  }
  
  b ~ dnorm(0,1.0E-6)
  g.ef ~ dnorm(0,1.0E-6) # or ge for grazing
  f.ef ~ dnorm(0,1.0E-6) # or fe for wildfire
  cb.ef ~ dnorm(0,1.0E-6) #
  logit(rate) <- b
 # rate.g <-exp(b+ge) # swa g for f for grazing
  int ~ dnorm(0,1.0E-6)
  cover1970 <- ilogit(int)
  
  for (m in 1:2) {
  prec[m] <- 1/sd[m]/sd[m]	
  sd[m] ~ dunif(0,100)
  }
  
} #'

AlpVeg <- jags(veg.data, veg.inits, 
               c('int', 'b',  'sd',  'g.ef','pg','f.ef',
                'p.veg', 'rate'), 
                 Model, 3, 1000) 

 

print(AlpVeg)
par(mfrow = c(4, 2))
traceplot(AlpVeg, mfrow=c(4, 2), varname = c('int', 'b', 'sd','g.ef','f.ef', 'cb.ef'), ask=FALSE)
attach.jags(AlpVeg)
par(mfrow = c(4,2))
plot(density(int))
plot(density(b))
plot(density(g.ef))
plot(density(f.ef))
plot(density(cb.ef))
plot(density(sd[, 1]))
plot(density(sd[, 2]))


AlpVeg.upd <- update(AlpVeg, n.iter=10000)
save(AlpVeg, file="AlpVeg.Rdata")
rm(AlpVeg.upd.upd)
#only want $BUGS.output and $sims.list within that for density plots etc
load('AlpVeg.upd.Rdata')
#traceplot
plot( AlpVeg.upd$BUGSoutput$sims.list$mn[,1], type='l')
#denisty plot
plot(density(AlpVeg$BUGSoutput$sims.list$mn[,1]))

lb.p.veg <- apply(AlpVeg$BUGSoutput$sims.list$p.veg,2,function(x) quantile(x,prob=0.025))
ub.p.veg <- apply(AlpVeg$BUGSoutput$sims.list$p.veg,2,function(x) quantile(x,prob=0.975))
mn.p.veg <- colMeans(AlpVeg$BUGSoutput$sims.list$p.veg)
#mn.p.veg <- apply(AlpVeg$BUGSoutput$sims.list$p.veg,2,mean)

lb.pg <- apply(AlpVeg$BUGSoutput$sims.list$pg,2,function(x) quantile(x,prob=0.025))
ub.pg <- apply(AlpVeg$BUGSoutput$sims.list$pg,2,function(x) quantile(x,prob=0.975))
mn.pg <- colMeans(AlpVeg$BUGSoutput$sims.list$pg)

lb.pcb <- apply(AlpVeg$BUGSoutput$sims.list$pcb,2,function(x) quantile(x,prob=0.025))
ub.pcb <- apply(AlpVeg$BUGSoutput$sims.list$pcb,2,function(x) quantile(x,prob=0.975))
mn.pcb <- colMeans(AlpVeg$BUGSoutput$sims.list$pcb)

mn.p <- colMeans(AlpVeg$BUGSoutput$sims.list$p)

windows()
plot.new()
plot.window(xlim=c(-5,65),ylim=c(0,0.8))
box()
axis(1, at=(0:6)*10,labels=((0:6)*10)+1940)
axis(2, las=1)

for (i in 1:23) { 
  points(x=dyw[i,]-1940, y=dow[i,], type='l', col = 'grey80', lwd=2) 
}



points(1:60,mn.p.veg, type='l', lwd=2)
points(1:60,lb.p.veg, type='l',lty=2, lwd=2)
points(1:60,ub.p.veg, type='l',lty=2, lwd=2)

points(1:60,mn.pg, type='l', col = 'red', lwd=2)
points(1:60,lb.pg, type='l',lty=2, col = 'red', lwd=2)
points(1:60,ub.pg, type='l',lty=2, col = 'red', lwd=2)

points(1:60,mn.pcb, type='l', col = 'blue', lwd=2)
points(1:60,lb.pcb, type='l',lty=2, col = 'blue',lwd=2)
points(1:60,ub.pcb, type='l',lty=2, col = 'blue',lwd=2)

mtext('Years',1, line=3)
mtext('Proportion',2, line=3)

dev.off()
   
   
   