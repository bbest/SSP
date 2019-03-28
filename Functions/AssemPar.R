## AssemPar: Script to estimate parameters of the assemblage
#Install and load the following libraries and code:
library(vegan)
library(ape)

AssemPar<-function(data, type, spatial.str = FALSE, coords, moranI.alfa=0.05, spatial.unsee.sp= "positive"){
  if (type == "P/A"){
    Y <- 1 * (data>0) #Ensure data is of presence/absence
    richness<-specpool(Y)
    Par<-NULL
    Par$type <-"P/A"
    Par$Sobs <- as.numeric(richness[1]) #Species observed
    Par$Sest <-round(as.numeric(richness[2])) # Extrapoled richness using Chao2
    fo<-apply(Y, 2, sum)
    fo<-sort(fo, decreasing = TRUE)
    sp<-fo[1:Par$Sobs]
    sp.names<-names(sp)
    Y<-subset(Y, select= sp.names)
    data<-subset(data, select= sp.names) #ensure exclude unregistered species in data 
    Par$n<-nrow(Y) # number of samples
    if (Par$Sest > Par$Sobs){
      Par$par<-data.frame(Species=c(rep(NA, Par$Sest)), fo=c(rep(NA, Par$Sest)), mu=c(rep(NA, Par$Sest)), Spatial.cor=c(rep(NA, Par$Sest)))
      Par$par[,1]<-c(variable.names(Y),c(rep("unseen species", Par$Sest-Par$Sobs)))
      Par$par[1:Par$Sobs,2]<-sp #frequency of occurrence of each species
      Par$par[(Par$Sobs+1):Par$Sest,2]<-rep(Par$par[Par$Sobs,2],Par$Sest-Par$Sobs)
      Par$par[,3]<- Par$par[,2]/Par$n
      
    } else {
      Par$par<-data.frame(Species=c(rep(NA, Par$Sobs)), fo=c(rep(NA, Par$Sobs)), mu=c(rep(NA, Par$Sest)), Spatial.cor=c(rep(NA, Par$Sobs)))
      Par$par[,1]<- variable.names(Y)
      Par$par[,2]<- fo# frequency of occurrence of each species
      Par$par[,3]<- Par$par[,2]/Par$n
    }
    Par$sumFo<-sum(Par$par[,2])
    
    if (spatial.str == TRUE){ #test for spatial autocorrelation
      d<-as.matrix(dist(coords)) 
      w<-1/d
      diag(w)<-0
      moranI.test<-NULL 
      for (i in 1:Par$Sobs){
        moranI.test[[i]]<-tryCatch(Moran.I(data[,i], w), error= function(e){
          z<-list(observed=0, expected=0, sd=0, p.value=1)
          return(z)
        })     
      }
      for (i in 1:Par$Sobs){
        if (Par$par$fo[i]>=3 & moranI.test[[i]][4] < moranI.alfa & as.numeric(moranI.test[[i]][1]) > as.numeric(moranI.test[[i]][2])){
          Par$par$Spatial.cor[i]<-"positive"
        }
        if (Par$par$fo[i]>=3 & moranI.test[[i]][4] < moranI.alfa & as.numeric(moranI.test[[i]][1]) < as.numeric(moranI.test[[i]][2])){
          Par$par$Spatial.cor[i]<-"negative"
        }
        if (is.na(Par$par$Spatial.cor[i]))
        {Par$par$Spatial.cor[i]<-"random"}
        Par$par[(Par$Sobs+1):Par$Sest,4]<-rep(spatial.unsee.sp, Par$Sest-Par$Sobs)
      }
      
    } else {
      for (i in 1:Par$Sobs){
        Par$par$Spatial.cor[i]<-"random"
      }
      Par$par[(Par$Sobs+1):Par$Sest,4]<-rep(spatial.unsee.sp, Par$Sest-Par$Sobs)
    }
    
  }        
  
  if (type == "counts"){
    Y <- data
    richness<-specpool(Y)
    Par<-NULL
    Par$type <-"counts"
    Par$Sobs <- as.numeric(richness[1]) #Species observed
    Par$Sest <-round(as.numeric(richness[2])) # Extrapoled richness using Chao2
    fo<-apply(Y, 2, sum)
    fo<-sort(fo, decreasing = TRUE)
    sp<-fo[1:Par$Sobs]
    sp.names<-names(sp)
    Y<-subset(Y, select= sp.names)
    data<-subset(data, select= sp.names)
    Par$n<-nrow(Y) # number of samples
    if (Par$Sest > Par$Sobs){
      Par$par<-data.frame(Species=c(rep(NA, Par$Sest)), fo=c(rep(NA, Par$Sest)), mean=c(rep(NA, Par$Sest)), 
                          var=c(rep(NA, Par$Sest)), D=c(rep(NA, Par$Sest)), prob=c(rep(NA, Par$Sest)), d=c(rep(NA, Par$Sest)), 
                          Spatial.cor=c(rep(NA, Par$Sest)))
      Ypa <- 1 * (Y>0) 
      M.0rm<-replace(Y, Y==0, NA)# remove zero by NA
      Par$par[,1]<-c(variable.names(Y),c(rep("unseen species", Par$Sest-Par$Sobs)))
      Par$par[1:Par$Sobs,2]<-(apply(Ypa, 2, sum))
      Par$par[(Par$Sobs+1):Par$Sest,2]<-rep(Par$par[Par$Sobs,2],Par$Sest-Par$Sobs)
      Par$par[1:Par$Sobs,3]<-(apply(M.0rm, 2, mean, na.rm=TRUE))
      Par$par[(Par$Sobs+1):Par$Sest,3]<- Par$par[(Par$Sobs+1):Par$Sest,4]<-rpois(Par$Sest-Par$Sobs,mean(Par$par$mean[1:Par$Sobs]))
      Par$par[1:Par$Sobs,4]<-(apply(M.0rm, 2, var, na.rm=TRUE))
      for (i in 1:Par$Sobs){
        if (is.na(Par$par$var[i])){
          Par$par$var[i]<-0
        }
      }
      DW<-dispweight(Y)
      Par$par[1:Par$Sobs,5]<-attr(DW,"D")
      Par$par[(Par$Sobs+1):Par$Sest,5]<-c(rep(1, Par$Sest-Par$Sobs))
      Par$par[1:Par$Sobs,6]<-attr(DW,"p")
      for (i in 1:Par$Sobs){
        if (Par$par$fo[i]>=3 & Par$par$D[i] > 1 & Par$par$prob[i]< 0.05 & Par$par$var[i] > Par$par$mean[i] & Par$par$mean[i]>0){
          Par$par$d[i]<-(Par$par$var[i] - Par$par$mean[i]) / Par$par$mean[i]^2
        } else {Par$par$d[i]<-0}
      }
      Par$par$d[(Par$Sobs+1):Par$Sest]<-0
      if (spatial.str == TRUE){ #test for spatial autocorrelation
        d<-as.matrix(dist(coords)) 
        w<-1/d
        diag(w)<-0
        moranI.test<-NULL 
        for (i in 1:Par$Sobs){
          moranI.test[[i]]<-Moran.I(data[,i], w)     
        }
        for (i in 1:Par$Sobs){
          if (Par$par$fo[i]>=3 & moranI.test[[i]][4] < moranI.alfa & as.numeric(moranI.test[[i]][1]) > as.numeric(moranI.test[[i]][2])){
            Par$par$Spatial.cor[i]<-"positive"
          }
          if (Par$par$fo[i]>=3 & moranI.test[[i]][4] < moranI.alfa & as.numeric(moranI.test[[i]][1]) < as.numeric(moranI.test[[i]][2])){
            Par$par$Spatial.cor[i]<-"negative"
          }
          if (is.na(Par$par$Spatial.cor[i]))
          {Par$par$Spatial.cor[i]<-"random"}
        }
        Par$par[(Par$Sobs+1):Par$Sest,8]<-rep(spatial.unsee.sp, Par$Sest-Par$Sobs)
      } else {
        for (i in 1:Par$Sobs){
          Par$par$Spatial.cor[i]<-"random"
        }
        Par$par[(Par$Sobs+1):Par$Sest,8]<-rep(spatial.unsee.sp, Par$Sest-Par$Sobs)
      }
      
    } else {
      Par$par<-data.frame(Species=c(rep(NA, Par$Sobs)), fo=c(rep(NA, Par$Sobs)), mean=c(rep(NA, Par$Sobs)), 
                          var=c(rep(NA, Par$Sobs)), D=c(rep(NA, Par$Sobs)), prob=c(rep(NA, Par$Sobs)),d=c(rep(NA, Par$Sobs)),
                          Spatial.cor=c(rep(NA, Par$Sobs)))
      Ypa <- 1 * (Y>0) 
      M.0rm<-replace(Y, Y==0, NA)# remove zero by NA
      Par$par[,1]<-variable.names(Y)
      Par$par[,2]<-(apply(Ypa, 2, sum))
      Par$par[,3]<-(apply(M.0rm, 2, mean, na.rm=TRUE))
      Par$par[,4]<-(apply(M.0rm, 2, var, na.rm=TRUE))
      for (i in 1:Par$Sobs){
        if (is.na(Par$par$var[i])){
          Par$par$var[i]<-0
        }
      }
      DW<-dispweight(Y)
      Par$par[,5]<-attr(DW,"D")
      Par$par[,6]<-attr(DW,"p")
      for (i in 1:nrow(Par$par)){## porqu? este es diferente?
        if (Par$par$fo[i]>=3 & Par$par$D[i] > 1 & Par$par$prob[i]< 0.05 & Par$par$var[i] > Par$par$mean[i] & Par$par$mean[i]>0){
          Par$par$d[i]<-(Par$par$var[i] - Par$par$mean[i]) / Par$par$mean[i]^2
        } else {Par$par$d[i]<-0}
      }
      if (spatial.str == TRUE){ #test for spatial autocorrelation
        d<-as.matrix(dist(coords)) 
        w<-1/d
        diag(w)<-0
        moranI.test<-NULL 
        for (i in 1:Par$Sobs){
          moranI.test[[i]]<-Moran.I(data[,i], w)     
        }
        for (i in 1:Par$Sobs){
          if (Par$par$fo[i]>=3 & moranI.test[[i]][4] < moranI.alfa & as.numeric(moranI.test[[i]][1]) > as.numeric(moranI.test[[i]][2])){
            Par$par$Spatial.cor[i]<-"positive"
          }
          if (Par$par$fo[i]>=3 & moranI.test[[i]][4] < moranI.alfa & as.numeric(moranI.test[[i]][1]) < as.numeric(moranI.test[[i]][2])){
            Par$par$Spatial.cor[i]<-"negative"
          }
          if (is.na(Par$par$Spatial.cor[i]))
          {Par$par$Spatial.cor[i]<-"random"}
        }
      } else {
        for (i in 1:Par$Sobs){
          Par$par$Spatial.cor[i]<-"random"
        }
      }
    }
  }
  if (type == "cover"){
    Y <- data
    richness<-specpool(Y)
    Par<-NULL
    Par$type <-"cover"
    Par$Sobs <- as.numeric(richness[1]) #Species observed
    Par$Sest <-round(as.numeric(richness[2])) # Extrapoled richness using Chao2
    fo<-apply(Y, 2, sum)
    fo<-sort(fo, decreasing = TRUE)
    sp<-fo[1:Par$Sobs]
    sp.names<-names(sp)
    Y<-subset(Y, select= sp.names)
    data<-subset(data, select= sp.names)
    Par$n<-nrow(Y) # number of samples
    if (Par$Sest > Par$Sobs){
      Par$par<-data.frame(Species=c(rep(NA, Par$Sest)), fo=c(rep(NA, Par$Sest)), 
                          logmean=c(rep(NA, Par$Sest)), logvar=c(rep(NA, Par$Sest)), Spatial.cor=c(rep(NA, Par$Sest)))
      Ypa <- 1 * (Y>0) 
      M.0rm<-replace(Y, Y==0, NA)# remove zero by NA
      Par$par[,1]<-c(variable.names(Y),c(rep("unseen species", Par$Sest-Par$Sobs)))
      Par$par[1:Par$Sobs,2]<-(apply(Ypa, 2, sum))
      Par$par[(Par$Sobs+1):Par$Sest,2]<-rep(Par$par[Par$Sobs,2],Par$Sest-Par$Sobs)
      Par$par[1:Par$Sobs,3]<-(apply(log(M.0rm), 2, mean, na.rm=TRUE))
      Par$par[(Par$Sobs+1):Par$Sest,3]<- Par$par[(Par$Sobs+1):Par$Sest,4]<-rpois(1,mean(Par$par$logmean[1:Par$Sobs]))
      Par$par[1:Par$Sobs,4]<-(apply(log(M.0rm), 2, var, na.rm=TRUE))
      for (i in 1:Par$Sobs){
        if (is.na(Par$par$logvar[i])){
          Par$par$logvar[i]<-0
        }
      }
      if (spatial.str == TRUE){ #test for spatial autocorrelation
        d<-as.matrix(dist(coords)) 
        w<-1/d
        diag(w)<-0
        moranI.test<-NULL 
        for (i in 1:Par$Sobs){
          moranI.test[[i]]<-Moran.I(data[,i], w)     
        }
        for (i in 1:Par$Sobs){
          if (Par$par$fo[i]>=3 & moranI.test[[i]][4] < moranI.alfa & as.numeric(moranI.test[[i]][1]) > as.numeric(moranI.test[[i]][2])){
            Par$par$Spatial.cor[i]<-"positive"
          }
          if (Par$par$fo[i]>=3 & moranI.test[[i]][4] < moranI.alfa & as.numeric(moranI.test[[i]][1]) < as.numeric(moranI.test[[i]][2])){
            Par$par$Spatial.cor[i]<-"negative"
          }
          if (is.na(Par$par$Spatial.cor[i]))
          {Par$par$Spatial.cor[i]<-"random"}
        }
        Par$par[(Par$Sobs+1):Par$Sest,5]<-rep(spatial.unsee.sp, Par$Sest-Par$Sobs)
      } else {
        for (i in 1:Par$Sobs){
          Par$par$Spatial.cor[i]<-"random"
        }
        Par$par[(Par$Sobs+1):Par$Sest,5]<-rep(spatial.unsee.sp, Par$Sest-Par$Sobs)
      } 
    } else {
      Par$par<-data.frame(Species=c(rep(NA, Par$Sobs)), fo=c(rep(NA, Par$Sobs)), 
                          logmean=c(rep(NA, Par$Sobs)), logvar=c(rep(NA, Par$Sobs)), Spatial.cor=c(rep(NA, Par$Sobs)))
      Ypa <- 1 * (Y>0) 
      M.0rm<-replace(Y, Y==0, NA)# remove zero by NA
      Par$par[,1]<-c(variable.names(Y))
      Par$par[,2]<-(apply(Ypa, 2, sum))
      Par$par[,3]<-(apply(log(M.0rm), 2, mean, na.rm=TRUE))
      Par$par[,4]<-(apply(log(M.0rm), 2, var, na.rm=TRUE))
      for (i in 1:Par$Sobs){
        if (is.na(Par$par$logvar[i])){
          Par$par$logvar[i]<-0
        }}
      if (spatial.str == TRUE){ #test for spatial autocorrelation
        d<-as.matrix(dist(coords)) 
        w<-1/d
        diag(w)<-0
        moranI.test<-NULL 
        for (i in 1:Par$Sobs){
          moranI.test[[i]]<-Moran.I(data[,i], w)     
        }
        for (i in 1:Par$Sobs){
          if (Par$par$fo[i]>=3 & moranI.test[[i]][4] < moranI.alfa & as.numeric(moranI.test[[i]][1]) > as.numeric(moranI.test[[i]][2])){
            Par$par$Spatial.cor[i]<-"positive"
          }
          if (Par$par$fo[i]>=3 & moranI.test[[i]][4] < moranI.alfa & as.numeric(moranI.test[[i]][1]) < as.numeric(moranI.test[[i]][2])){
            Par$par$Spatial.cor[i]<-"negative"
          }
          if (is.na(Par$par$Spatial.cor[i]))
          {Par$par$Spatial.cor[i]<-"random"}
        }
        Par$par[(Par$Sobs+1):Par$Sest,5]<-rep(spatial.unsee.sp, Par$Sest-Par$Sobs)
      } else {
        for (i in 1:Par$Sobs){
          Par$par$Spatial.cor[i]<-"random"
        }
      }
    }
  }
  return(Par)
}

##Description
# Assempar is a function that extracts the main parameters of the pilot data using basic R functions as well as functions like specpool{Vegan} and Moran.I{ape}.
# The function returns an object of class List, to be used by SimDat.

#Usage
AssemPar(data, type, spatial.str, coords, moranI.alfa, spatial.unsee.sp)

## Arguments
# data                  Dataframe giving P/A or abundance (counts or coverage) of species.

# type                  Nature of the data to be processed. It may be presence / absence ("P/A"),
#                       counts of individuals ("counts"), or coverage ("cover")

# spatial.str           Logical argument indicating whether spatial information (coordinates) of samples is available. If set TRUE, a matrix or dataframe with X and Y
#                       coordinates must be supplied. This file must have the same number of rows than data to be used by the function Moran.I{ape}.
#                       The analisys is only valid for individual counts (argument type = "counts"), and will gave an Error if two or more samples share the same 
#                       X and Y values.

# coords                Matrix or dataframe given Y and X coordinates of samples.

# moranI.alfa           Probability for rejecting random spatial distribution. Default is < 0.05.

#spatial.unsee.sp       Spatial distribution to be assumed for unseen species. It might be "random", "positive" or "negative".



