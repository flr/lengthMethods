library(FLCore)
library(FLBRP)
library(FLasher)
library(ggplotFL)
library(FLife)
library(mydas)

library(popbio)
library(plyr)
library(reshape)

library(LBSPR)

setwd("/home/laurence-kell/Desktop/flr/lengthMethods/tests")

source("../R/oemLn.R")
source("../R/lbspr.R")

load("../data/fpp_om.Rdata")
load("../data/fpp_om.Rdata")

## needed for selection pattern a50 and ato95
par=lhPar(VBpars,units="NA")

## needed to estimate m/k 
priors=popdyn(par)

set.seed(789)
srDev=rlnoise(10,FLQuant(0,dimnames=list(year=2020:2040)),.5,b=0)
params(sr)=iter(params(sr),seq(10))

ak   =invAlk(par,cv=0.1,age=0:10)  

f    =FLQuant(0.2,dimnames=list(year=2020:2040,iter=seq(10)))
om   =fwd(iter(stk,seq(10)),fbar=f,sr=sr,residuals=srDev)

lfd =lenSample(catch.n(om)[,ac(2011:2030)],ak,nsample=250)  

ggplot(melt(lfd[,ac(seq(2011,2030,5))]))+ 
  geom_histogram(aes(len,weight=value),binwidth=1)+
  facet_grid(year~iter,scale="free")+
  xlab("Length (cm)")+ylab("Frequency")+
  coord_cartesian(xlim=c(0,mean(lh["linf"])))

lb=lbspr(lfd,priors)  

plot(lb)


## L-CC
source('~/Desktop/flr/lengthMethods/R/haupt.R')


z=haupt(lfd,priors,priors["sel1"]*0.9)

## MLZ
library(MLZ)

eq=lhEql(par)
om=as(eq,"FLStock")

source('~/Desktop/flr/mydas/R/omOut.R')

ts   =omSmry(om,eq,par)
mnLen=as.FLQuant(with(ts,data.frame(data=cln,year=year,iter=iter)))

source('~/Desktop/flr/mydas/R/mlz.R')

res=mlz(mnLen[,ac(91)],priors)
