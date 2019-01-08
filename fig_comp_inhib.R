###############################################################################
###############################################################################
#A figure produce with a dummy dataset of V.inaequalis type
###############################################################################
###############################################################################

#loading the libraries
library(drc)
library(plotrix)
library(gdata)

#load the global dataset
creadat<-read.table(file="data/reg_inhib.txt",header=T,sep="\t")
creadat$species<-factor(creadat$species,levels=rev(levels(creadat$species)))
collist<-c("forestgreen","black")


###############################################################################
#This is the figure for the boscalid
###############################################################################

#subsetting the global dataset
bosc.dat<-creadat[creadat$active_substance=="boscalid",]

#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
bosc_rez<-as.character(bosc.dat[bosc.dat$dose=="30" & bosc.dat$perc_croiss>50,
                                "sample_ID"])
REZbos<-data.frame("strain_ID"=as.character(),"ED50"=as.character())
bosc.dat<-drop.levels(bosc.dat,reorder=FALSE)

i<-1
datatemp<-bosc.dat[bosc.dat$strain_ID==names(table(bosc.dat$strain_ID))[i],]
couleur<-collist[as.numeric(datatemp$strain_type)]
typeline<-as.numeric(datatemp$species)[1]
if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
  tempx<-data.frame("strain_ID"=names(table(bosc.dat$strain_ID))[i],
                    "ED50"=c("NA"))
  REZbos<-rbind(REZbos,tempx)
  plot(0,1,ylim=c(-10,120),xlim=c(0,30),col="red",
       main=names(table(bosc.dat$strain_ID))[i])
} else {
  temp.m1<-drm(perc_croiss~dose,
               data=datatemp,
               fct=LL.4())
  plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col="red",lwd=3,
       lty=typeline,main=names(table(bosc.dat$strain_ID))[i])
  plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col="red",
       lty=typeline,type="confidence",add=TRUE)
  temp<-ED(temp.m1,50,type="absolute")
  tempx<-data.frame("strain_ID"=names(table(bosc.dat$strain_ID))[i],
                    "ED50"=temp[1])
  REZbos<-rbind(REZbos,tempx)
}

i<-2
datatemp<-bosc.dat[bosc.dat$strain_ID==names(table(bosc.dat$strain_ID))[i],]
couleur<-collist[as.numeric(datatemp$strain_type)]
typeline<-as.numeric(datatemp$species)[1]
if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
  tempx<-data.frame("strain_ID"=names(table(bosc.dat$strain_ID))[i],
                    "ED50"=c("NA"))
  REZbos<-rbind(REZbos,tempx)
  plot(0,1,ylim=c(-10,120),xlim=c(0,30),col="darkgreen",
       main=names(table(bosc.dat$strain_ID))[i],add=TRUE)
} else {
  temp.m1<-drm(perc_croiss~dose,
               data=datatemp,
               fct=LL.4())
  plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col="darkgreen",lwd=3,
       lty=typeline,main=names(table(bosc.dat$strain_ID))[i],add=TRUE)
  plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col="darkgreen",
       lty=typeline,type="confidence",add=TRUE)
  temp<-ED(temp.m1,50,type="absolute")
  tempx<-data.frame("strain_ID"=names(table(bosc.dat$strain_ID))[i],
                    "ED50"=temp[1])
  REZbos<-rbind(REZbos,tempx)
}

#export to pdf 6 x 5 inches


###############################################################################
#END
###############################################################################