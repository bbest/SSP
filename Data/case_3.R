# Intertidal rocky shores

library(vegan)
library (ape)
library(sampling)
library(ggplot2)

##Import file La_Pared.csv
dat1 <- read.csv(choose.files(), row.names=1)

##Import file Chirimena.csv
dat2 <- read.csv(choose.files(), row.names=1)

#####################################
##Estimation of parameters using pilot data and function AssemPar

#La Pared
dat.par1<-AssemPar (data = dat1,
                    type= "cover",
                    spatial.str = FALSE,
                    moranI.alfa=0.05,
                    spatial.unsee.sp= "random")

#Chirimena
dat.par2<-AssemPar (data = dat2,
                    type= "cover",
                    spatial.str = FALSE,
                    moranI.alfa=0.05,
                    spatial.unsee.sp= "random")

##Simulation of data
#La Pared
dat.sim1<-SimData(dat.par1, cases= 100, n= 100, sites = 1)

#Chirimena
dat.sim2<-SimData(dat.par2, cases= 100, n= 100, sites = 1)

##Sampling and estimation of MSE for each combination of sites and samples with function SampSD
#La Pared
samples1<-SampSD(dat.sim1,
                 Sest = 13,
                 transformation = "square root",
                 method = "bray",
                 multi.site = FALSE,
                 n=100,
                 p.n = 50,
                 sites = 1,
                 p.s = 1,
                 k=100)
#Chirimena
samples2<-SampSD(dat.sim2,
                 Sest = 20,
                 transformation = "square root",
                 method = "bray",
                 multi.site = FALSE,
                 n=100,
                 p.n = 50,
                 sites = 1,
                 p.s = 1,
                 k=100)

##Get the means MSE for each sampling design with funcion datplot and then plot with ggplot
#La Pared
mse.means1<-MSEplot(results = samples1, multi.site = FALSE)

#Chirimena
mse.means2<-MSEplot(results = samples2, multi.site = FALSE)

##Quality of data

qua.1<-DatQuality(pilotdata = dat1,
                simdata = dat.sim1,
                assempar = dat.par1,
                transformation = "square root",
                method = "bray")


qua.2<-DatQuality(pilotdata = dat2,
                simdata = dat.sim2,
                assempar = dat.par2,
                transformation = "square root",
                method = "bray")

##Plot MultSE

mse.means1$site<-rep("La Pared", nrow(mse.means1))
mse.means2$site<-rep("Chirimena", nrow(mse.means2))
mse.means<-rbind(mse.means1, mse.means2)
mse.means

fig4<-ggplot(mse.means, aes(x=levels, y=mean, colour=site))+
                geom_point(position = position_dodge(0.5), size = 0.5)+
                geom_errorbar(aes(ymin=lower, ymax=upper),position = position_dodge(0.5), size=0.1, width=.2)+
                theme_bw(base_size=16) +
                ylab ("Multivariate pseudo SE")+
                xlab("Sampling effort (n)")+
                scale_y_continuous(breaks=seq(0.0, 0.3, 0.025))+
                scale_x_discrete(breaks=seq(4, 50, 2))+
                theme(axis.text.x = element_text(colour="black", size=rel(0.7)),
                axis.text.y = element_text(colour="black", size=rel(0.7)),
                axis.title.x = element_text(colour="black", size=rel(0.9)),
                axis.title.y = element_text(colour="black", size=rel(0.9)),
                panel.grid.major = element_line(size=0.2),
                panel.grid.minor = element_blank(),
                panel.border = element_rect(size=0.2),
                axis.ticks= element_line(size=0.2),
                legend.position=c(1,1),
                legend.justification=c(1,1),
                legend.background=element_rect(fill="white", size=0.2, colour="black"))+
                labs(colour="Sites")  
fig4

