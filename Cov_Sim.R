##This script is set to run remotely via a Linux bash script
## Hence the "CommandArgs" lines

args <- commandArgs(TRUE)
Start <- as.numeric(args[1])
Stop <- as.numeric(args[2])

#sim.Burrs <- dget(file = 'Sim_data_05112018.txt')
sim.Burrs <- dget(file = 'Sim_data_05072018.txt')
for (i in Start:Stop) {
  source('RJMCMC.R')
  require(Rdistance)
  library(runjags)
  library(rjags)
All.Data <- data.frame(sim.Burrs[[i]]$all.burr)
Burrs <- data.frame(dist = All.Data$dist, size = (All.Data$diam)/10, Found = All.Data$found, Trans = All.Data$trans)
Found <- subset(Burrs, Burrs$Found == 1)
n.transects <- sim.Burrs[[i]]$n.lines
n <- nrow(Burrs)

require(Rdistance)
detection.data <- data.frame(siteID = as.factor(Found$Trans), groupsize = 1, dist = Found$dist)
transect.data <- data.frame(siteID = as.factor(c(1:n.transects)), length = sim.Burrs[[i]]$len)
area = 10000
est <- F.automated.CDA(detection.data, transect.data, w.lo=0, w.hi=max(Found$dist),
	likelihoods=c("halfnorm", "uniform", "Gamma"),
	series=c("cosine", "hermite", "simple"), expansions=0:3, warn=TRUE,
	area=area, ci=0.95, R=500, by.id=FALSE, plot.bs=FALSE, plot=FALSE)

pt.est.null <- est$n.hat
ci.est.null <- unname(est$ci)
truth <- (n/100 <= ci.est.null[2] & n/100 >= ci.est.null[1])

#######################
#bayesian method
modelstring.Foo = "
  model
{
  for (i in 1:(nind +nz)) {
  w[i] ~ dbern(psi)					##augmentation 
  x[i] ~ dunif(0,Bx)					#distance from line for the missed ones; Bx = max(distances) 
  
  z[i] ~ dgamma(a[i],c.beta[i])T(5,50)		                
  a[i] <- shape[clust[i]]		## a depends on which cluster the size comes from
  c.beta[i] <- betaClust[clust[i]]	#precision is allowed to vary between clusters
  clust[i] ~ dcat( pClust[1:Nclust] )	## clusters are defined as categories, assuming Nclust clusters
  
  
  sigma[i] <- exp(sigma.int+sigma.beta*z[i])	#log link
  logp[i] <- -((x[i]*x[i])/(2*sigma[i]*sigma[i]))	#This is the half normal distribution with an adjustment for covariates
  p[i] <- exp(logp[i])*xi[i]				# probabilty of seeing the thing, regardless of us actually seeing it; scaled by xi
  xi[i] <- ifelse(z[i] < b.point, m*z[i]+intercept, 1) #if less than b.point, probability is linear; if larger than b.point, perfect detection on line
  
  mu[i] <- w[i]*p[i] 					## probabilty of seeing it for all the ones we DID see (where w[i] = 1)
  y[i] ~ dbern(mu[i])
  
  o[i]~ dbin(o2[i], 1)
  logit(o2[i]) <- o.int + z.beta*z[i]
  }			 		
  
  p.online ~dunif(0.2, 0.7)		#estimated detection on line for 4 cm objects
  m <- (1-p.online)/(b.point-4) 	#slope for detection on the line for smaller objects		
  intercept <- p.online-(4*m)	## finding intercept via the detection of the 4cm object 
  
  
  sigma.int~ dnorm(0,s.int.tau)
  s.int.tau <- 1/(s.int.sd*s.int.sd)
  s.int.sd ~ dunif(.00001,5)	
  sigma.beta~ dnorm(0,s.beta.tau)
  s.beta.tau <- 1/(s.beta.sd*s.beta.sd)
  s.beta.sd ~ dunif(.00001,5)
  o.int~ dnorm(0,o.int.tau)
  o.int.tau <- 1/(o.int.sd*o.int.sd)
  o.int.sd ~ dunif(.00001,5)		
  z.beta~ dnorm(0,z.beta.tau)
  z.beta.tau <- 1/(z.beta.sd*z.beta.sd)
  z.beta.sd ~ dunif(.00001,5)			
  
  for (clustIdx in 1: Nclust) {
  shape[clustIdx] ~ dunif(1,100)
  betaClust[clustIdx] ~ dunif(.2,2)
  }
  pClust[1:Nclust] ~ ddirch(psizes) ## probability of each cluster = the probability of that category in the ddirch distribution
  
  
  psi~ dunif(0,1)			#exists or not		
  b.point ~ dunif(18,22)		#where the stick breaks
  
  Occ <- sum(o)/(nind+nz)
  N <- sum(w)		
  D <- N/(2*L*Bx)			#burrow density
  Nt <- N*Occ 			 
  Dt <- Nt/(2*L*Bx)			#tort density
  
}
"

nind <- nrow(Found)
nz <-  ifelse(nrow(Burrs) >= 400, 550, 275)

# 
# w.clust <- c(0.50, 0.1, 0.33, .03, .04)
# shape.clust <- c(10, 20, 30, 40, 66)
# rate.clust <- c(1,.7,1, 2, 2)
# y.clust <- rep(Found$size, 3)
# Z <- do.call(cbind, lapply(1:5, function(j)
#   w.clust[j]*dgamma(y.clust, shape.clust[j], rate = rate.clust[j])))
# Z <- apply(Z, 1, function(x) which(x==max(x))[1])
# res <- mixgam.rjmcmc(y = y.clust, nsweep = 80000, kmax = 10, k = 5, w = w.clust, 
#                      shape = shape.clust, rate = rate.clust, Z,
#                      delta=2, xi=NULL, kappa=NULL, alpha=NULL,
#                      beta=NULL, g=NULL, h=NULL, verbose=TRUE)
# 
# ksave <- res$k.save
# groups <- round(table(ksave[-(1:40000)])/40000,4)
# groups
# Nclust <- as.numeric(names(groups)[which(groups == max(groups))])+1
# if(length(Nclust) != 1){Nclust <- sample(Nclust,1)}
#Nclust <- 2

L = sim.Burrs[[i]]$total.len*10^-4
x = c(Found$dist,rep(NA,nz))		
Bx = max(x, na.rm = TRUE)
y = c(rep(1,nind), rep(0,nz))
z <- c(Found$size,rep(NA,nz))
w = c(rep(1,nind), rep(NA,nz))

clust = rep(NA,(nind +nz)) 	# no idea which cluster anything is in, so unknown	
clust[which.min(z)]=1 # smallest value assigned to cluster 1; generally represents juvis
clust[which.max(z)]=Nclust # highest value assigned to largest cluster value; generally represents large adults


### Run the model 
jd.Foo = list(nind= nind, Nclust= Nclust, nz = nz, z= z, L = L, Bx = Bx, x=x, y=y, onesRepNclust = rep(1,Nclust), clust = clust)
ji.Foo <- list(sigma.int = 1, sigma.beta = .025, w= y, pClust = rep(1/Nclust, Nclust), 
               shape = c(13, 22, 48, 33), betaClust =c(2,.6, 1, .5)) 
jp.Foo <- c("pClust","D", "N", "sigma.int", "sigma.beta", "p.online", "shape", "betaClust")
Foo <- autorun.jags(modelstring.Foo, data = jd.Foo, inits = list(ji.Foo, ji.Foo, ji.Foo), 
	monitor= jp.Foo, startburnin = 1000, startsample = 5000, adapt = 200, 
	method = 'parallel',max.time = '1m', n.chain = 3, psrf.target = 1.2, thin = 1, silent.jags = FALSE, summarise = TRUE)

Covmod <- summary(Foo)


results = list()
results$simulation <- paste("Dataset #", i)
results$true <- c(paste('Truth'), n/(1000*1000*10^-4))
results$DISTANCE <- c(paste('Distance'), ci.est.null, pt.est.null)
results$FoundN <- nind
results$Covsum <- Covmod
results

lapply(results, function(x) write.table(data.frame(x), 'Covariate_Simulation_Clusts.csv'  , append= T, sep=',' ))
}
