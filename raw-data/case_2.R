##Sponges from Alacranes National Park, Gulf of Mexico

library(vegan)
library (ape)
library(sampling)
library(ggplot2)

#import files sponges_data.csv and sponges_coord.csv
dat <- sponges_data
dat.xy <- sponges_coord

#####################################
#Estimation of parameters using pilot data and function AssemPar

dat.par<-AssemPar (data = dat, 
                   type= "counts",
                   spatial.str = TRUE,
                   coords = dat.xy,
                   moranI.alfa=0.05,
                   spatial.unsee.sp= "random")

#Simulation of data
ptmSIM <- proc.time()
dat.sim<-SimData(dat.par, cases= 10, n= 20, sites = 20)
proc.time() - ptmSIM


#Sampling and estimation of MSE for each combination of sites and samples with function SampSD

ptmSAMP <- proc.time()
samples<-SampSD(dat.sim, 
                Sest = 50,
                transformation = "square root",
                method = "bray",
                multi.site = TRUE,
                n=20, 
                p.n = 20,
                sites = 20,
                p.s = 20,
                k=50)

proc.time() - ptmSAMP

#Get the means MSE for each sampling design with funcion datplot and then plot with ggplot

ptmMSE <- proc.time()

mse.means<-MSEplot(results = samples,
                   multi.site = TRUE)

proc.time() - ptmMSE

mse.means[20:37,5]<-"transects"

fig3<-ggplot(mse.means, aes(x=levels, y=mean))+
        geom_point()+
        geom_errorbar(aes(ymin=lower, ymax=upper), width=.1)+
        facet_grid(.~sv)+
        theme_bw(base_size=16) +
        ylab ("Multivariate pseudo SE")+ 
        xlab("Sampling effort")+
        scale_y_continuous(breaks=seq(0.0, 0.3, 0.05))+
        theme(axis.text.x = element_text(colour="black", size=rel(0.7)),
        axis.text.y = element_text(colour="black", size=rel(0.7)),
        axis.title.x = element_text(colour="black", size=rel(0.9)),
        axis.title.y = element_text(colour="black", size=rel(0.9)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=0.4),
        axis.ticks= element_line(size=0.2))
fig3


#quality of data


qua.sp<-DatQuality(pilotdata = dat,
                     simdata = dat.sim,
                     assempar = dat.par,
                     transformation = "square root",
                     method = "bray")
