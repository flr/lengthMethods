if (FALSE){
ms<-foreach(i=toRun, 
             .combine=list,
             .multicombine=TRUE,
             .export=c("srDev","lhs","design","haupt"),
             .packages=c("FLCore","ggplotFL","FLBRP","FLasher",
           "FLife","mydas","FLCandy","JABBA","mpb",
           "popbio","spatstat",
           "plyr","dplyr","reshape","GGally")) %dopar% {
             
      fl=paste("/home/laurence-kell/Desktop/papers/COM3/R/runs/om/om",i,"RData",sep=".") 
      load(fl)
      maxf=max(fbar(om))*1.000 
             
      om.lmean=om
      rm(om)
             
      set.seed(789)
             
      par        =lhs[[design[i,"Stock"]]]
      par["s"]   =design[i,"s"]
      par["sel2"]=design[i,"sel2"]
      par["sel3"]=design[i,"sel3"]
      ak         =invAlk(par,cv=0.1)  
             
      #r=c(popdyn(lhs[[design[i,"Stock"]]])["r"])          
             
      for (iYr in 106:130){ 
         lfd  =lenSample(catch.n(om.lmean)[,ac(iYr-(5:1))],ak,nsample=design[i,"nsample"])
         ind=transform(subset(as.data.frame(lfd,drop=TRUE),data>0),
         
         wt  =c(par["a"])*len^c(par["b"]),
         lopt=c(2/3*par["linf"]))
         
         ind=ddply(ind, .(year,iter), with, lenInd(len,data,wt,lopt))
         ind=cbind(ind,linf=c(par["linf"]),l50=c(par["l50"]))
         lmean=as.FLQuant(transmute(ind,year=year,iter=iter,data=1.5*lmean/par["linf",drop=T]))
         
         ctc=hcrSBTD(iYr,control=FLPar(c(k1=7.5,k2=7.5,gamma=1)),
                     index =lmean,
                     catch =catch(om.lmean)[,ac(iYr-(1))])
         om.lmean=fwd(om.lmean,catch=ctc,sr=eq,residuals=srDev[[design[i,"CV"]]],maxF=maxf)
     
      save(iYr,om.lmean,file=paste("/home/laurence-kell/Desktop/papers/COM3/R/runs/om/om.lmean",i,"RData",sep="."))
    }
           }
}
