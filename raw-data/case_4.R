# Epibiont on mangrove roots

library(vegan)
library (ape)
library(sampling)
library(ggplot2)

#Import data, and split by sector
dat <- read.csv(choose.files(), row.names=1, sep=",")
dat.xy <- read.csv(choose.files(), row.names=1, sep=",")

dat1<-dat[dat$sector=="EXTERNAL",3:154]
dat2<-dat[dat$sector=="INTERMEDIATE",3:154]
dat3<-dat[dat$sector=="INTERNAL",3:154]

dat1.xy<-dat.xy[dat.xy$sector=="EXTERNAL",3:4]
dat2.xy<-dat.xy[dat.xy$sector=="INTERMEDIATE",3:4]
dat3.xy<-dat.xy[dat.xy$sector=="INTERNAL",3:4]

#####################################
#Estimation of parameters using pilot data and function AssemPar

dat.par1<-AssemPar (data = dat1, 
                    type= "P/A",
                    spatial.str = TRUE,
                    coords = dat1.xy,
                    moranI.alfa=0.05,
                    spatial.unsee.sp= "positive")

dat.par2<-AssemPar (data = dat2, 
                    type= "P/A",
                    spatial.str = TRUE,
                    coords = dat2.xy,
                    moranI.alfa=0.05,
                    spatial.unsee.sp= "positive")

dat.par3<-AssemPar (data = dat3, 
                    type= "P/A",
                    spatial.str = TRUE,
                    coords = dat3.xy,
                    moranI.alfa=0.05,
                    spatial.unsee.sp= "positive")


#Simulation of data
ptmSIM <- proc.time()
dat.sim1<-SimData(dat.par1, cases= 100, n= 30, sites = 20)
dat.sim2<-SimData(dat.par2, cases= 100, n= 30, sites = 20)
dat.sim3<-SimData(dat.par3, cases= 100, n= 30, sites = 20)
proc.time() - ptmSIM


#Sampling and estimation of MSE for each combination of sites and samples with function SampSD

ptmSAMP <- proc.time()
samples1<-SampSD(dat.sim1, 
               Sest = 107,
               transformation = "P/A",
               method = "jaccard",
               multi.site = TRUE,
               n=30, 
               p.n = 20,
               sites = 20,
               p.s = 20,
               k=20)

samples2<-SampSD(dat.sim2, 
                 Sest = 72,
                 transformation = "P/A",
                 method = "jaccard",
                 multi.site = TRUE,
                 n=30, 
                 p.n = 20,
                 sites = 20,
                 p.s = 20,
                 k=20)

samples3<-SampSD(dat.sim3, 
                 Sest = 56,
                 transformation = "P/A",
                 method = "jaccard",
                 multi.site = TRUE,
                 n=30, 
                 p.n = 20,
                 sites = 20,
                 p.s = 20,
                 k=2)
proc.time() - ptmSAMP

#Get the means MSE for each sampling design with funcion datplot and then plot with ggplot

mse.means1<-MSEplot(results = samples1, multi.site = TRUE)
mse.means2<-MSEplot(results = samples2, multi.site = TRUE)
mse.means3<-MSEplot(results = samples3, multi.site = TRUE)

mse.means1$sector<-rep("External", nrow(mse.means1))
mse.means2$sector<-rep("Intermediate", nrow(mse.means2))
mse.means3$sector<-rep("Internal", nrow(mse.means3))
mse.means<-rbind(mse.means1, mse.means2, mse.means3)


##Plot MultSE

fig5<-ggplot(mse.means, aes(x=levels, y=mean, colour=sector))+
        geom_point(position = position_dodge(0.5), size = 1)+
        geom_errorbar(aes(ymin=lower, ymax=upper),position = position_dodge(0.5), size=0.4, width=.4)+
        facet_grid(.~sv)+
        theme_bw(base_size=16) +
        ylab ("Multivariate pseudo SE")+ 
        xlab("Sampling effort")+
        scale_y_continuous(breaks=seq(0.0, 0.3, 0.05))+
        theme(axis.text.x = element_text(colour="black", size=rel(0.7)),
        axis.text.y = element_text(colour="black", size=rel(0.7)),
        axis.title.x = element_text(colour="black", size=rel(0.9)),
        axis.title.y = element_text(colour="black", size=rel(0.9)),
        panel.grid.major = element_line(size=0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=0.4),
        axis.ticks= element_line(size=0.2),
        legend.position=c(1,1),
        legend.justification=c(1,1),
        legend.background=element_rect(fill="white", size=0.2, colour="black"))+
        labs(colour="Sector")  
fig5

#quality of data

qua.1<-DatQuality(pilotdata = dat1,
                     simdata = dat.sim1,
                     assempar = dat.par1,
                     transformation = "P/A",
                     method = "jaccard")

qua.2<-DatQuality(pilotdata = dat2,
                  simdata = dat.sim2,
                  assempar = dat.par2,
                  transformation = "P/A",
                  method = "jaccard")

qua.3<-DatQuality(pilotdata = dat3,
                  simdata = dat.sim3,
                  assempar = dat.par3,
                  transformation = "P/A",
                  method = "jaccard")