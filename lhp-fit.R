# here we demonstrate how to estimate the spatial scale of phenotypic variation
# for a single species from data simulated under our model. this procedure can
# be repeated for a second species to obtain upper and lower bounds on the
# spatial scale of coevolution (as defined in our paper).

# this script has been adapted from example script provided by 
# RandomFields documentation:
#                     https://rdrr.io/cran/RandomFields/man/RandomFields.html

library(RandomFields)

#
# section 1: define and simulate GRFs for each species
#

# dimensions of the total region under consideration
width = 10
height = 10

# global/colocated mean trait variances
V.h = 10
V.p = 10

# individual trait variance at a given point in space
v.h = 0.1
v.p = 0.1

# spatial scales of intraspecific variation
lambda.h = 2
lambda.p = 0.2

# this defines a GRF with Matern spatial autocovarince function
# with parameter nu=1 in agreement with our results for intraspecific
# spatial variation. to simplify this demonstration, we assume
# independent GRFs for each species
model.h = RMmatern(nu=1, var=V.h, scale=lambda.h)
model.p = RMmatern(nu=1, var=V.p, scale=lambda.p)

# defines the grid over which the random fields will be simulated
x = seq(0,width,.1)
y = seq(0,height,.1)

# this provides a single (n=1) realization for each simulated random field
z.h <- RFsimulate(model.h, x, y, n=1)
z.p <- RFsimulate(model.p, x, y, n=1)

#
# section 2: take survey sample from the simulated GRF
#

# number of individuals sampled for each species during initial survey
h.num = 100
p.num = 100

# take sample
xy.h <- coordinates(z.h)
xy.p <- coordinates(z.p)
pts.h <- sample(nrow(xy.h), min(h.num, nrow(xy.p) / 2))
pts.p <- sample(nrow(xy.p), min(p.num, nrow(xy.h) / 2))
dta.h <- matrix(nrow=nrow(xy), as.vector(z.h))[pts.h, ]
dta.p <- matrix(nrow=nrow(xy), as.vector(z.p))[pts.p, ]

# this step adds sampling noise due to variation among individuals
# at each location
dta.h = dta.h + rnorm(h.num,sd=sqrt(v.h))
dta.p = dta.p + rnorm(p.num,sd=sqrt(v.p))

dta.h <- cbind(xy[pts.h, ], dta.h)
dta.p <- cbind(xy[pts.p, ], dta.p)

# visualize the sampled points
plot(z.h, dta.h)
plot(z.p, dta.p)

#
# section 3: estimate model parameters (spatial scale in
#                             particular) from simulated data
#

# in the model being estimated, we fix nu, but let
# the colocated variance and spatial scale be fit 
# via maximum likelihood. we only define this one time
# since the form of the model is the same for each species
estmodel <- RMmatern(nu=1, var=NA, scale=NA)

# these lines do the actual fitting (may take a while)
(fit.h <- RFfit(estmodel, data=dta.h))
(fit.p <- RFfit(estmodel, data=dta.p))

# these are the estimated spatial scales for each species
lmda.h = fit.h@table$ml[1]
lmda.p = fit.p@table$ml[1]

#
# section 4: identify sample sites
#

# this demonstrates that samples taken 3 times the spatial scale
# will have only small statistical dependence
matern(3,nu=1)

# this demonstrates that samples taken within 1/4 the spatial scale
# will be highly statistically dependent
matern(1/4,nu=1)

# to approximate the covariance, we then toss down a grid of sample locations
# that are 3 times the upper bound of the spatial scale of coevolution 
# (i.e., 3*max(lmda.h,lmda.p)) apart. we then sample individuals at each site,
# each with a radius 1/8 of the lower bound of the spatial scale of coevolution
# (i.e., min(lmda.h,lmda.p)/8), where 1/8 is used since the diameter
# will then be min(lmda.h,lmda.p)/4. with these samples we compute local mean 
# traits for each species at each samples site. finally, we then compute
# the covariance of mean traits across sample sites.

# lower and upper bounds on the spatial scale of coevolution
lwr.bd = min(lmda.h,lmda.p)
upr.bd = max(lmda.h,lmda.p)

# note, we add lwr.bd/4 to make sure the edges of the sampling sites are
# 3*upr.bd apart.
xs = seq(0,width,by=(3*upr.bd+lwr.bd/4))
ys = seq(0,height,by=(3*upr.bd+lwr.bd/4))

# we then center these sample sites in the total region considered
xs = xs + max(xs)/2
ys = ys + max(ys)/2

#
# section 5: sample individuals in each species at the sample sites chosen
#

# these are the number of individuals sampled at each site
hs = 20
ps = 20

# loop through each site and draw individual sample locations at each site
h.locs = c()
p.locs = c()
for(x in xs){
  for(y in ys){
    xmn = x-lwr.bd/8
    xmx = x+lwr.bd/8
    ymn = y-lwr.bd/8
    ymx = y+lwr.bd/8
    new.h = cbind(runif(hs,xmn,xmx),runif(hs,ymn,ymx))
    new.p = cbind(runif(ps,xmn,xmx),runif(ps,ymn,ymx))
    h.locs = rbind(h.locs,new.h)
    p.locs = rbind(p.locs,new.p)
  }
}

# here we sample individual traits at each individual location
z.h = RFsimulate(model.h, h.locs[,1], h.locs[,2], n=1)
z.p = RFsimulate(model.p, p.locs[,1], p.locs[,2], n=1)

# total number of individuals sampled per species
h.tot = hs*length(xs)*length(ys)
p.tot = ps*length(xs)*length(ys)

# this adds sampling noise due to trait variation among individuals
zh = as.vector(z.h) + rnorm(h.tot,sd=v.h)
zp = as.vector(z.p) + rnorm(p.tot,sd=v.p)

#
# section 6: compute local mean traits per species at each sample site
#

zh.bar = c()
h.part = seq(1,h.tot,by=hs)
for(i in 1:length(h.part)){
  zhb = mean(zh[h.part[i]:(h.part[i]+hs-1)])
  zh.bar = c(zh.bar, zhb)
}

zp.bar = c()
p.part = seq(1,ps*length(xs)*length(ys),by=ps)
for(i in 1:length(p.part)){
  zpb = mean(zp[p.part[i]:(p.part[i]+ps-1)])
  zp.bar = c(zp.bar, zpb)
}

#
# section 7: estimating the interspecific colocated covariance C_HP(0)
#

CHP = cov(cbind(zh.bar,zp.bar))[2,1]

print(paste("Estimated interspecific colocated covariance:",CHP,sep=" "))

# note: under the assumptions made in this script, this covariance will be
# zero on average.
