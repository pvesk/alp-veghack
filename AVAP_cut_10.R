###set up a dataframe for the observations
setwd('~/Documents/research/alpine air photo veg change')

veg <- read.csv('Alpvegcov.csv')   # read in from data folder
veg <- veg[with(veg,order(veg$Site, veg$date)), ]
veg[,3:5] <- round(veg[,3:5]) 
veg[,6] <- 1600 - veg[,3] - veg[,4] - veg[,5]
veg[,2] <- floor((veg$date-1926)/10)
################this recodes dates into 10yr intervals###################

names(veg[,6]) <- "other.hits"
sites <- read.csv('Alpsites.csv')   # read in from data folder
#yrs <- veg$date-1970
hd <- floor((veg$date-1926)/10)
veg <- merge(veg,sites)
veg <- veg[with(veg,order(veg$Site, veg$date)), ]
veg[,7] <- floor((veg[,7]-1926)/10)
first <- min(veg[,2])
last <- max(veg[,2])

ord <- rep(1:3,times=23)
veg <- cbind(veg,ord)
veg.w <- reshape(veg, timevar="ord", idvar="Site", direction="wide")

all.yrs.a <- rep(NA,times = (23) * (last - first + 1) * 4)
dim(all.yrs.a) <- c(23,(last - first + 1),4)

#years of observations
y.o.obs <- rep(NA, times=23*3) #empty vectors
dim(y.o.obs) <- c(23,3)
y.o.obs[,1] <- veg.w[,2]#this is the year of the first observation
y.o.obs[,2] <- veg.w[,15]
y.o.obs[,3] <- veg.w[,28]
yoo1 <- y.o.obs[,1]


all.yrs.a <- rep(NA,times = (23) * (last - first +1) * 4)
dim(all.yrs.a) <- c(23,(last - first + 1),4)
o<- rep(NA, times=23)

for (k in 1:4) {
for (i in 1:23) {
 o[i] <- veg.w[i, k+2]  #this is the 1st observation for veg type k
 all.yrs.a[i,  y.o.obs[i,1] - first + 1, k] <- o[i]
}
for (i in 1:23) {
  o[i] <- veg.w[i, k+15]  #this is the 2nd observation
  all.yrs.a[i,  y.o.obs[i,2] - first + 1, k] <- o[i]  
}
for (i in 1:23) {
  o[i] <- veg.w[i, k+28]  #this is the 3rd observation
  all.yrs.a[i, y.o.obs[i,3] - first + 1, k] <- o[i] 
}
}

# setup covariates for time since grazing, wildfires and clearing
tsg <- floor((sites$last.grazed-1926)/10)
t <- all.yrs.a[,,1]
for (i in 1:23){
  for (j in 1 : (last - first + 1) ) {
    t[i, j] <- ifelse( (j + first - 2) > tsg[i] , 1, 0) 
  }
} 
# now for CB
cb.t <- t
cb.t[,1:(last - first + 1)] <- 0
for (i in 1:23) {
  o[i] <- ifelse(sites[i, 3] > 0, (sites[i, 3]/100), 0)  #this is the observed redcution
#  cb.t[i, y.o.obs[i,1]:(y.o.obs[i,1]+1)] <- o[i]  # 2 period effect of C&B
  cb.t[i, y.o.obs[i,1]] <- o[i]
}
#cb1 <- cb.t
#cb.t[,2:69] <- 1
for (i in 1:23) {
  o[i] <- ifelse(sites[i, 4] > 0, (sites[i, 4]/100), 0)  #this is the 2nd observation
#  cb.t[i, y.o.obs[i,2]:(y.o.obs[i,2]+1)] <- o[i]  
  cb.t[i, y.o.obs[i,2]] <- o[i]  
}


# for wildfire, lets just assume 15 y before photo
fy.t <- t
fy.t[ , 1:(last - first + 1)] <- 0
for (i in 1:23) {
  o[i] <- ifelse(sites[i, 5] > 0, 1, 0)  #this is the 2nd observation
  fy.t[i, y.o.obs[i,2] - 1] <- o[i]  
}

#now for 2nd interval

for (i in 1:23) {
  o[i] <- ifelse(sites[i, 6] > 0, 1, 0)  #this is the 3rd observation
  fy.t[i, y.o.obs[i,3] - 1] <- o[i]  
}

library(R2jags)
#library(runjags)
#source('~/insects model/doJAGS.R')



KK.init <- matrix(runif(92,0,4), nrow=23)
for (i in 1:23) {   KK.init[i,1] <- log(max(veg.w[i,c(3,16,29)])+2)  }
for (i in 1:23) {   KK.init[i,2] <- log(max(veg.w[i,c(4,17,30)])+2)  }
for (i in 1:23) {   KK.init[i,3] <- log(max(veg.w[i,c(5,18,31)])+2)  }
for (i in 1:23) {   KK.init[i,4] <- log(max(veg.w[i,c(6,19,32)])+2)  }

interval <- veg.w[,15] - veg.w[,2] 
interval2 <- veg.w[,28] -  veg.w[,15] 

dif <- veg.w[,16] - veg.w[,3] 
dif <- cbind(dif, veg.w[,17] - veg.w[,4] ) 
dif <- cbind(dif, veg.w[,18] - veg.w[,5] ) 
dif <- cbind(dif, veg.w[,19] - veg.w[,6] ) 
dif <- cbind(dif, veg.w[,29] - veg.w[,16] ) 
dif <- cbind(dif, veg.w[,30] - veg.w[,17] ) 
dif <- cbind(dif, veg.w[,31] - veg.w[,18] ) 
dif <- cbind(dif, veg.w[,32] - veg.w[,19] ) 

N.init <- all.yrs.a
for (i in 1:23) {
  for (j in 1:(last - first + 1)) {
    for (k in 1:4) {
      N.init[i,j,k] <- ifelse(j < y.o.obs[i,1]+1, NA, 
                              (ifelse(is.na(N.init[i,j,k]), max(round(veg.w[i, k+2] + dif[i,k] / interval[i] * (j - y.o.obs[i,1]) ),1), N.init[i,j,k] ) ) )
      N.init[i,j,k] <- ifelse(j < y.o.obs[i,2], N.init[i,j,k],  
                              max(round(veg.w[i, k+15] + dif[i,k+3] / interval2[i] * (j - y.o.obs[i,2])),1))
      N.init[i,j,k] <- ifelse(j == y.o.obs[i,2], NA, N.init[i,j,k])
      N.init[i,j,k] <- ifelse(j == y.o.obs[i,3], NA, N.init[i,j,k])      
      N.init[i,j,4] <- ifelse(is.na(N.init[i,j,4]), NA, round(N.init[i,j,4]*.7))
  }}}


veg.inits <- function() list(
 # mn.K = c(3,3,3,3),
  alpha.100 = c(5,9,0.5,1.8), #rnorm(3, mean = 0.05, sd = 0.01),
  N = N.init,
 # KK = KK.init, #runif(23, 500, 1500),
  sd = runif(2,0.5,1), 
  gr.ef = rnorm(4, mean = 0.5, sd = 0.1),
 KT = rnorm(4, mean = 0.5, sd = 0.1),
  cb.ef = rnorm(4, mean = 0.5, sd = 0.1),
  f.ef = rnorm(4, mean = 0.0, sd = 0.1),
 a = c(3,2.5,3.5, 3),
 b = c(0,1,-1, 0)
)   #THIS IS THE INITIALS LIST

#define data
veg.data <- list( 
  Nsites = 23,
 # alt = sites[,7],
  precip = sites[,8],
  Nveg = 4,
  T = last - first,
  yoo1 = yoo1,
  N = all.yrs.a,
  sg = as.matrix(t[ , 1:(last - first + 1)]),
  cb.t = as.matrix(cb.t[ , 1:(last - first + 1)]),
  fy.t = as.matrix(fy.t[ , 1:(last - first + 1)])) ## , 
# THIS IS THE DATA LIST

###lets try a Ricker function
cat(' model{
for (j in 1:Nsites) {
  
    for (t in yoo1[j] : T) { 
      
      for (k in 1:Nveg) {
        Er[j, t, k] <- exp( alpha.100[k] / 100 * ( 1 - N[j, t, k] / K[j, t, k])  ) 

        N[j, t + 1, k] ~ dpois(lambda[j, t, k])
        log(K[j, t, k]) <- log.K[j, t, k]
        log.K[j, t, k] ~ dnorm(mu.K[j,t,k], sd1^-2)T(0,7.4) 
# spatio-temporal variation aorund the linear predicted K (grazing, fire)
    mu.K[j,t,k]  <-  KK[j, k] + KT[k] * (t-4) + gr.ef[k] * sg[j, t]
    } 

for ( k in 1: 2) {
        lambda[j, t, k] <- (N[j, t, k] + 0.5) * ( Er[j, t, k] ) * exp(cb.ef[k] * cb.t[j,t])
}
for ( k in 3: 4) {
        lambda[j, t, k] <- (N[j, t, k] + 0.5) * ( Er[j, t, k] ) 
}

 #   mu.K[j,t,1]  <-  KK[j, 1] + gr.ef[1] * sg[j, t]  # + f.ef[1] * fy.t[j, t] +  cb.ef[1] * cb.t[j, t] 
  #     mu.K[j,t,2]  <-  KK[j, 2]  + gr.ef[2] * sg[j, t]#  + f.ef[2] * fy.t[j, t]   +  cb.ef[2] * cb.t[j, t]
   #    mu.K[j,t,3]  <-  KK[j, 3]  + gr.ef[3] * sg[j, t] #  + f.ef[3] * fy.t[j, t]   + cb.ef[3] * cb.t[j, t]
    #   mu.K[j,t,4]  <-  KK[j, 4]  + 0 * fy.t[j, t] + 0 *gr.ef[4] * sg[j, t]  + 0 *cb.ef[4] * cb.t[j, t] 
#log.K[j,t,4] <- log(1600 - exp(log.K[j,t,1]) - exp(log.K[j,t,2]) - exp(log.K[j,t,3]))
  }
}

for (k in 1:Nveg) { 
  alpha.100[k] ~ dnorm(0.05,0.7)T(0, )
  gr.ef[k] ~ dnorm(0, 0.05) # for grazing removal 
  cb.ef[k] ~ dnorm(0, 0.05)T(-1,1)
#  f.ef[k] ~ dnorm(0, 0.05)

  a[k] ~ dnorm(4.5,.7)T(0,7.38)
  b[k] ~ dnorm(0,0.01)
KT[k] ~ dnorm(0,0.01)
}

  for (j in 1: Nsites) { 
    for (k in 1: 4) {
      mn.K[j,k] <- a[k] + b[k] * (precip[j] - 1465)/216 # (alt[j]-1433)/185  #standardize
      KK[j,k]  ~ dnorm(mn.K[j,k], sd2^-2)T(0,7.38)   #site level departures from predicted K from Precip
}

}
sd1 ~ dt(0, .016, 1)T(0, 7)
sd2 ~ dt(0, .016, 1)T(0, 7)
}  ## End model
'  , file=(modelfile <- tempfile()))



###fit model
AlpVeg10 <- jags(veg.data, veg.inits, 
                c('KT','sd1','sd2', 'sd', 'alpha.100', 'gr.ef', 'cb.ef', 'f.ef', 'KK', 'a', 'b'),
                modelfile, n.chains=3, n.iter=20000) #, n.burnin=5000)  



AlpVeg10.upd <- update(AlpVeg10, 
                      params= c('KT','sd1','sd2', 'alpha', 'gr.ef', 'cb.ef', 'f.ef', 'KK', 'a', 'b'),
                      n.iter=10000, n.thin=10)
print(AlpVeg10)
par(mfrow = c(4, 2))
traceplot(AlpVeg10, mfrow=c(4, 2), varname = c( 'K', 'sd',  'gr.ef', 'cb.ef', 'f.ef', 'a', 'b'), ask=FALSE)
attach.jags(AlpVeg10)
par(mfrow = c(4,2))
plot(density(a[,1]))
plot(density(a[,2]))
plot(density(a[,3]))

plot(KK[,19,1], KK[,19,2])
plot(KK[,19,3], KK[,19,2])
plot(KK[,19,1], KK[,19,3])

plot(KT[,1], gr.ef[,1])
abline(h=0,v=0)
plot(KT[,2], gr.ef[,2])
abline(h=0,v=0)
plot(KT[,3], gr.ef[,3])
abline(h=0,v=0)
plot(KT[,4], gr.ef[,4])
 abline(h=0,v=0)


## cb vs fire
plot(cb.ef[,1], f.ef[,1])
abline(h=0,v=0)


plot(a[,1], b[,1])
plot(a[,1], alpha.100[,1])
plot(a[,1], cb.ef[,1])
plot(a[,1], gr.ef[,1])
plot(a[,1], f.ef[,1])
plot(cb.ef[,1], f.ef[,1])
plot(gr.ef[,1], f.ef[,1])
plot(alpha.100[,1], b[,1])
plot(alpha.100[,1], cb.ef[,1])
plot(alpha.100[,1], gr.ef[,1])
plot(alpha.100[,1], f.ef[,1])
plot(a[,1], KK[,19,1])
plot(b[,1], KK[,19,1])
plot(alpha.100[,1], KK[,19,1])
plot(cb.ef[,1], KK[,19,1])
plot(gr.ef[,1], KK[,19,1])
plot(f.ef[,1], KK[,19,1])


plot(a[,2], b[,2])
plot(a[,2], alpha.100[,2])
plot(a[,2], cb.ef[,2])
plot(a[,2], gr.ef[,2])
plot(a[,2], f.ef[,2])
plot(cb.ef[,2], f.ef[,2])
plot(gr.ef[,2], f.ef[,2])
plot(alpha.100[,2], b[,2])
plot(alpha.100[,2], cb.ef[,2])
plot(alpha.100[,2], gr.ef[,2])
plot(alpha.100[,2], f.ef[,2])
plot(a[,2], KK[,19,2])
plot(b[,2], KK[,19,2])
plot(alpha.100[,2], KK[,19,2])
plot(cb.ef[,2], KK[,19,2])
plot(gr.ef[,2], KK[,19,2])
plot(f.ef[,2], KK[,19,2])


plot(a[,3], b[,3])
plot(a[,3], alpha.100[,3])
plot(a[,3], cb.ef[,3])
plot(a[,3], gr.ef[,3])
plot(a[,3], f.ef[,3])
plot(cb.ef[,3], f.ef[,3])
plot(gr.ef[,3], f.ef[,3])
plot(alpha.100[,3], cb.ef[,3])
plot(alpha.100[,3], gr.ef[,3])
plot(alpha.100[,3], f.ef[,3])


plot(a[,3], KK[,19,3])
plot(b[,3], KK[,19,3])
plot(alpha.100[,3], KK[,19,3])
plot(cb.ef[,3], KK[,19,3])
plot(gr.ef[,3], KK[,19,3])
plot(f.ef[,3], KK[,19,3])

plot(density(sd[, 1]))
plot(density(alpha.100[, 1]))
plot(density(gr.ef[, 1]))
plot(density(cb.ef[, 1]))


save(AlpVegRL.upd, file="AlpVegRL.Rdata")
rm(AlpVegRL.upd)
#only want $BUGS.output and $sims.list within that for density plots etc
load('AlpVegRL.Rdata')












AlpVegRicker$BUGSoutput$sims.list$pred[,,1],2
plot(AlpVegRicker$BUGSoutput$sims.list$alpha, AlpVegRicker$BUGSoutput$sims.list$gr.ef)

## or with parallel
AlpVegRickerP <- run.jags(model, n.iter=10000, method = "rjparallel", inits = veg.inits, 
                          data=veg.data,  monitor = c('mn.K', 'K', 'sd', 'alpha', 'gr.ef', 'cb.ef1', 'cb.ef2') ) 

#with do.jags
AlpVegRickerP <- do.jags(.modelfile, n.iter=2000, parallel= TRUE, inits= veg.inits, n.chains=3,
                         data=veg.data,  parameters.to.save =c('mn.K', 'K', 'sd', 'alpha', 'gr.ef', 'cb.ef1', 'cb.ef2') ) 
