# micromollusk from 'Cayo Nuevo', Campeche bank, Gulf of Mexico, Yucatán, Mexico

library(vegan)
library (ape)
library(sampling)
library(ggplot2)

#import micromollusk.csv
dat <- read.csv(choose.files(), row.names=1)

#import micromollusk_spatial.csv
dat.xy <- read.csv(choose.files(), row.names=1)

#####################################
#Estimation of parameters using pilot data and function AssemPar

dat.par<-AssemPar (data = dat, 
                    type= "P/A",
                    spatial.str = TRUE,
                    coords = dat.xy,
                    moranI.alfa=0.05,
                    spatial.unsee.sp= "random")

#Simulation of data
ptmSIM <- proc.time()
dat.sim<-SimData(dat.par, cases= 100, n= 1000, sites = 1)
proc.time() - ptmSIM


#Sampling and estimation of MSE for each combination of sites and samples with function SampSD

ptmSAMP <- proc.time()
samples<-SampSD(dat.sim, 
               Sest = 98,
               transformation = "P/A",
               method = "jaccard",
               multi.site = FALSE,
               n=1000, 
               p.n = 100,
               sites = 1,
               p.s = 1,
               k=100)

proc.time() - ptmSAMP

#Get the means MSE for each sampling design with funcion datplot and then plot with ggplot

mse.means<-MSEplot(results = samples,
                   multi.site = FALSE)

#Figure 2
mse.means50<-mse.means[1:48,]
fig2<-ggplot(mse.means50, aes(x=levels, y=mean))+
        geom_point(size = 0.5)+
        geom_errorbar(aes(ymin=lower, ymax=upper), size=0.1, width=.2)+
        theme_bw(base_size=16) +
        ylab ("Multivariate pseudo SE")+ 
        xlab("Sampling effort (n)")+
        scale_y_continuous(breaks=seq(0.0, 0.3, 0.025))+
        scale_x_discrete(breaks=seq(4, 50, 2))+
        theme(axis.text.x = element_text(colour="black", size=rel(0.7)),
              axis.text.y = element_text(colour="black", size=rel(0.7)),
              axis.title.x = element_text(colour="black", size=rel(0.9)),
              axis.title.y = element_text(colour="black", size=rel(0.9)),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(size=0.4),
              axis.ticks= element_line(size=0.2))
fig2 

#quality of data

ptmQ <- proc.time()
qua.mico<-DatQuality(pilotdata = dat,
                     simdata = dat.sim,
                     assempar = dat.par,
                     transformation = "P/A",
                     method = "jaccard")
proc.time() - ptmQ