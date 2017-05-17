###############################################################################
###############################################################################
#code for the PBO tests ...
###############################################################################
###############################################################################

#loading the libraries
library(drc)
library(plotrix)

#set the working directory
setwd("~/work/Rfichiers/Githuber/probit_data")


###############################################################################
#comparison between repetitions
###############################################################################

#load the dataset
myzpbo<-read.table("pbo_280317.txt",header=T,sep="\t")

#in order to obtain the same results than with priprobit, the Finney equivalent
#method, we have to remove the constrain on lowerlimit and chose a log-normal 
#model instead of a log-logistic model
myzus.mpbo<-drm(dead/total~dose,weights=total,data=myzpbo,fct=LN.3u(),
                type="binomial")
plot(myzus.mpbo)
plot(myzus.mpbo,type="obs")
plot(myzus.mpbo,type="confidence")
plot(myzus.mpbo,type="obs",add=TRUE)
#the ED50 obtained is identical to the one obtained with priprobit, the SD 
#still differ a little bit
ed50val_mpbo<-ED(myzus.mpbo,50,interval="delta",reference="control")


op<-par(mfrow=c(2,1))
#let's evaluate the different model for each clone with and without pbo####
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                               myzpbo$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(myzpbo[,2]))[1])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                               myzpbo$pesticide=="thiaclopride_PBO",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                               myzpbo$pesticide=="thiaclopride" & 
                               myzpbo$date==levels(myzpbo$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo5<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[5],],
                  fct=LN.2(),
                  type="binomial")
plot(myzus.sspbo,type="obs",main=names(summary(myzpbo[,2]))[1])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue")
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2")
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="obs",add=TRUE,col="blue3")

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo5<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[5],],
                  fct=LN.2(),
                  type="binomial")
plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red")
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2")
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo4,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo4,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo5,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo5,type="obs",add=TRUE,col="red3")

#plot by date
op<-par(mfrow=c(5,1))
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[1],
                                       (levels(myzpbo$date)[1])))
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue")
plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red")
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[1],
                                       (levels(myzpbo$date)[2])))
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2")
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2")
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[1],
                                       (levels(myzpbo$date)[3])))
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3")
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3")
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[1],
                                       (levels(myzpbo$date)[4])))
plot(myzus.sspbo4,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="obs",add=TRUE,col="blue3")
plot(myzus.avpbo4,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo4,type="obs",add=TRUE,col="red3")
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[1],
                                       (levels(myzpbo$date)[5])))
plot(myzus.sspbo5,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="obs",add=TRUE,col="blue3")
plot(myzus.avpbo5,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo5,type="obs",add=TRUE,col="red3")
par(op)




#let's evaluate the different model for each clone with and without pbo####
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                               myzpbo$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(myzpbo[,2]))[2])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                               myzpbo$pesticide=="thiaclopride_PBO",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo5<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[5],],
                  fct=LN.2(),
                  type="binomial")
plot(myzus.sspbo,type="obs",main=names(summary(myzpbo[,2]))[2])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue")
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2")
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="obs",add=TRUE,col="blue3")

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo5<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[5],],
                  fct=LN.2(),
                  type="binomial")
plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red")
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2")
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo4,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo4,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo5,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo5,type="obs",add=TRUE,col="red3")

#plot by date
op<-par(mfrow=c(5,1))
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[2],
                                       (levels(myzpbo$date)[1])))
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue")
plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red")
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[2],
                                       (levels(myzpbo$date)[2])))
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2")
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2")
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[2],
                                       (levels(myzpbo$date)[3])))
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3")
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3")
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[2],
                                       (levels(myzpbo$date)[4])))
plot(myzus.sspbo4,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="obs",add=TRUE,col="blue3")
plot(myzus.avpbo4,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo4,type="obs",add=TRUE,col="red3")
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[2],
                                       (levels(myzpbo$date)[5])))
plot(myzus.sspbo5,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="obs",add=TRUE,col="blue3")
plot(myzus.avpbo5,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo5,type="obs",add=TRUE,col="red3")
par(op)




#let's evaluate the different model for each clone with and without pbo####
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                               myzpbo$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(myzpbo[,2]))[3])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                               myzpbo$pesticide=="thiaclopride_PBO",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo5<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[5],],
                  fct=LN.2(),
                  type="binomial")
plot(myzus.sspbo,type="obs",main=names(summary(myzpbo[,2]))[3])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue")
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2")
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="obs",add=TRUE,col="blue3")

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo5<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[5],],
                  fct=LN.2(),
                  type="binomial")
plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red")
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2")
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo4,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo4,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo5,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo5,type="obs",add=TRUE,col="red3")

#plot by date
op<-par(mfrow=c(5,1))
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[3],
                                       (levels(myzpbo$date)[1])))
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue")
plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red")
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[3],
                                       (levels(myzpbo$date)[2])))
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2")
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2")
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[3],
                                       (levels(myzpbo$date)[3])))
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3")
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3")
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[3],
                                       (levels(myzpbo$date)[4])))
plot(myzus.sspbo4,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="obs",add=TRUE,col="blue3")
plot(myzus.avpbo4,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo4,type="obs",add=TRUE,col="red3")
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[3],
                                       (levels(myzpbo$date)[5])))
plot(myzus.sspbo5,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="obs",add=TRUE,col="blue3")
plot(myzus.avpbo5,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo5,type="obs",add=TRUE,col="red3")
par(op)




#let's evaluate the different model for each clone with and without pbo####
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                               myzpbo$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(myzpbo[,2]))[4])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                               myzpbo$pesticide=="thiaclopride_PBO",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo5<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[5],],
                  fct=LN.2(),
                  type="binomial")
plot(myzus.sspbo,type="obs",main=names(summary(myzpbo[,2]))[4])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue")
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2")
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="obs",add=TRUE,col="blue3")

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo5<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[5],],
                  fct=LN.2(),
                  type="binomial")
plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red")
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2")
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo4,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo4,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo5,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo5,type="obs",add=TRUE,col="red3")

par(op)

#plot by date
op<-par(mfrow=c(5,1))
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[4],
                                       (levels(myzpbo$date)[1])))
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue")
plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red")
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[4],
                                       (levels(myzpbo$date)[2])))
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2")
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2")
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[4],
                                       (levels(myzpbo$date)[3])))
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3")
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3")
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[4],
                                       (levels(myzpbo$date)[4])))
plot(myzus.sspbo4,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="obs",add=TRUE,col="blue3")
plot(myzus.avpbo4,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo4,type="obs",add=TRUE,col="red3")
plot(myzus.sspbo,type="obs",main=paste(names(summary(myzpbo[,2]))[4],
                                       (levels(myzpbo$date)[5])))
plot(myzus.sspbo5,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="obs",add=TRUE,col="blue3")
plot(myzus.avpbo5,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo5,type="obs",add=TRUE,col="red3")
par(op)



###############################################################################
#sum up the results of the different tests
###############################################################################


#load the dataset
myzpbo<-read.table("pbo_280317.txt",header=T,sep="\t")
myzpbo<-cbind(myzpbo,"perc"=myzpbo$dead*100/myzpbo$total)
myzpbo<-cbind(myzpbo,"testID"=paste(myzpbo$clone,myzpbo$date))

#building a table summarising the results of the different repetition
recaprepet<-aggregate(dose~date+clone,
                      data=myzpbo[myzpbo$pesticide=="thiaclopride",],"mean")
recaprepet<-cbind("testID"=paste(recaprepet$clone,recaprepet$date),
                  recaprepet[,2:1])
#data of test with only thiaclopride
recaprepet<-cbind(recaprepet,"mort rate T- thia"=
  (aggregate(dead~date+clone,data=myzpbo[myzpbo$pesticide=="thiaclopride" 
                                      & myzpbo$dose==0,],"sum")$dead *100 / 
   aggregate(total~date+clone,data=myzpbo[myzpbo$pesticide=="thiaclopride" 
                                      & myzpbo$dose==0,],"sum")$total)
)
recaprepet<-cbind(recaprepet,"number T- thia"= aggregate(total~date+clone,
                  data=myzpbo[myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$dose==0,],"sum")$total)
recaprepet<-cbind(recaprepet,"number tested thia"= aggregate(total~date+clone,
                  data=myzpbo[myzpbo$pesticide=="thiaclopride",],"sum")$total)
#data of test with pbo added to thiaclopride
recaprepet<-cbind(recaprepet,"mort rate T- pbo"=
  (aggregate(dead~date+clone,data=myzpbo[myzpbo$pesticide!="thiaclopride" 
                                      & myzpbo$dose==0,],"sum")$dead *100 / 
   aggregate(total~date+clone,data=myzpbo[myzpbo$pesticide!="thiaclopride" 
                                      & myzpbo$dose==0,],"sum")$total)
)
recaprepet<-cbind(recaprepet,"number T- pbo"= aggregate(total~date+clone,
                  data=myzpbo[myzpbo$pesticide!="thiaclopride" & 
                                myzpbo$dose==0,],"sum")$total)
recaprepet<-cbind(recaprepet,"number tested pbo"= aggregate(total~date+clone,
                  data=myzpbo[myzpbo$pesticide!="thiaclopride",],"sum")$total)



#let's evaluate the different model for each clone with and without pbo
DL50<-c()
for (i in 1:19){
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$testID==recaprepet$testID[i] & 
                               myzpbo$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
ed50val<-ED(myzus.sspbo,50,interval="delta",reference="control")
DL50<-c(DL50,ed50val[1])
}

DL50pbo<-c()
for (i in 1:5){
  myzus.sspbo<-drm(dead/total~dose,weights=total,
                   data=myzpbo[myzpbo$testID==recaprepet$testID[i] & 
                                 myzpbo$pesticide!="thiaclopride",],
                   fct=LN.2(),
                   type="binomial")
  ed50val<-ED(myzus.sspbo,50,interval="delta",reference="control")
  DL50pbo<-c(DL50pbo,ed50val[1])
}
DL50pbo<-c(DL50pbo,NA)
for (i in 7:19){
  myzus.sspbo<-drm(dead/total~dose,weights=total,
                   data=myzpbo[myzpbo$testID==recaprepet$testID[i] & 
                                 myzpbo$pesticide!="thiaclopride",],
                   fct=LN.2(),
                   type="binomial")
  ed50val<-ED(myzus.sspbo,50,interval="delta",reference="control")
  DL50pbo<-c(DL50pbo,ed50val[1])
}

recaprepet<-cbind(recaprepet,"DL50"=DL50,"DL50pbo"=DL50pbo)
recaprepet<-cbind(recaprepet,"Synergy"=recaprepet$DL50/recaprepet$DL50pbo)
recaprepet[order(recaprepet$clone,as.Date(recaprepet$date,"%d-%b-%y")),]
write.table(recaprepet[order(recaprepet$clone,
                             as.Date(recaprepet$date,"%d-%b-%y")),],
            file="tablerecap2.txt",sep="\t",quote=FALSE,row.names=FALSE)


###############################################################################
#building the tables to make an analyze that merges the different repetition
###############################################################################

byclone<-aggregate(cbind(dead,total)~dose+pesticide+clone,data=myzpbo,"sum")

#without pbo
DL50clone<-c()
modelclone<-vector("list",4)
for (i in 1:4){
  myzus.sspbo<-drm(dead/total~dose,weights=total,
                   data=byclone[byclone$clone==levels(byclone$clone)[i] & 
                                  byclone$pesticide=="thiaclopride",],
                   fct=LN.2(),
                   type="binomial")
  modelclone[[i]]<-myzus.sspbo
  ed50val<-ED(myzus.sspbo,50,interval="delta",reference="control")
  DL50clone<-c(DL50clone,ed50val[1])
}
plot(modelclone[[1]],main=levels(byclone$clone)[1],type="bars",
     col="blue",ylim=c(-0.05,1.05))
plot(modelclone[[2]],type="average",col="green3",add=TRUE)
plot(modelclone[[3]],type="confidence",col="grey30",add=TRUE)
plot(modelclone[[4]],type="all",col="orange",add=TRUE)

#with pbo
DL50clonepbo<-c()
modelclonepbo<-vector("list",4)
for (i in 1:4){
  myzus.sspbo<-drm(dead/total~dose,weights=total,
                   data=byclone[byclone$clone==levels(byclone$clone)[i] & 
                                  byclone$pesticide!="thiaclopride",],
                   fct=LN.2(),
                   type="binomial")
  modelclonepbo[[i]]<-myzus.sspbo
  ed50val<-ED(myzus.sspbo,50,interval="delta",reference="control")
  DL50clonepbo<-c(DL50clonepbo,ed50val[1])
}
plot(modelclonepbo[[1]],main=levels(byclone$clone)[1],type="bars",
     col="blue",ylim=c(-0.05,1.05))
plot(modelclonepbo[[2]],type="average",col="green3",add=TRUE)
plot(modelclonepbo[[3]],type="confidence",col="grey30",add=TRUE)
plot(modelclonepbo[[4]],type="all",col="orange",add=TRUE)


#because there was important "natural" mortality in some tests, we remove 
#repetitions with a mortality rate greater than 25 % in either modality 
#(with or without pbo)

listrepcor<-recaprepet[c(which(recaprepet$`mort rate T- thia`>25),
                         which(recaprepet$`mort rate T- pbo` >25)),"testID"]
byclonecor<-aggregate(cbind(dead,total)~dose+pesticide+clone,
                      data=myzpbo[!myzpbo$testID %in% listrepcor,],"sum")

#we export the table
write.table(byclonecor,file="effectifpartest.txt",sep="\t",quote=FALSE,
            row.names=FALSE)

#another way which is simplier to do the model for every different clone
myzus.sspbocor<-drm(dead/total~dose,weights=total,clone,
                    data=byclonecor[byclonecor$pesticide=="thiaclopride",],
                    fct=LN.2(),type="binomial")
plot(myzus.sspbocor)
DL50clonecor<-ED(myzus.sspbocor,50,interval="delta",reference="control")

myzus.avpbocor<-drm(dead/total~dose,weights=total,clone,
                    data=byclonecor[byclonecor$pesticide!="thiaclopride",],
                    fct=LN.2(),type="binomial")
plot(myzus.avpbocor,add=TRUE,col="red")
DL50clonepbocor<-ED(myzus.avpbocor,50,interval="delta",reference="control")



###############################################################################
#New set of clones: comparison between repetitions
###############################################################################

#load the dataset
myzpbo<-read.table("pbo_120517.txt",header=T,sep="\t")


myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                               myzpbo$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(myzpbo[,2]))[1])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                               myzpbo$pesticide=="thiaclopride_PBO",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")


#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo5<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[5],],
                  fct=LN.2(),
                  type="binomial")
plot(myzus.sspbo,type="obs",main=names(summary(myzpbo[,2]))[1])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue")
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2")
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="obs",add=TRUE,col="blue3")


#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo5<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[1] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[5],],
                  fct=LN.2(),
                  type="binomial")
plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red")
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2")
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo4,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo4,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo5,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo5,type="obs",add=TRUE,col="red3")



#for the second clone
#let's evaluate the different model for each clone with and without pbo####
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                               myzpbo$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(myzpbo[,2]))[2])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                               myzpbo$pesticide=="thiaclopride_PBO",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")


#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo5<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[5],],
                  fct=LN.2(),
                  type="binomial")
plot(myzus.sspbo,type="obs",main=names(summary(myzpbo[,2]))[2])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue")
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2")
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="obs",add=TRUE,col="blue3")

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo5<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[2] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[5],],
                  fct=LN.2(),
                  type="binomial")
plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red")
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2")
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo4,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo4,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo5,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo5,type="obs",add=TRUE,col="red3")



#let's analysis the third clone
#let's evaluate the different model for each clone with and without pbo####
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                               myzpbo$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(myzpbo[,2]))[3])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                               myzpbo$pesticide=="thiaclopride_PBO",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo5<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[5],],
                  fct=LN.2(),
                  type="binomial")
plot(myzus.sspbo,type="obs",main=names(summary(myzpbo[,2]))[3])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue")
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2")
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="obs",add=TRUE,col="blue3")

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo5<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[3] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[5],],
                  fct=LN.2(),
                  type="binomial")
plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red")
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2")
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo4,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo4,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo5,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo5,type="obs",add=TRUE,col="red3")


#for the fourth clone
#let's evaluate the different model for each clone with and without pbo####
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                               myzpbo$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(myzpbo[,2]))[4])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                               myzpbo$pesticide=="thiaclopride_PBO",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo5<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride" & 
                                myzpbo$date==levels(myzpbo$date)[5],],
                  fct=LN.2(),
                  type="binomial")
plot(myzus.sspbo,type="obs",main=names(summary(myzpbo[,2]))[4])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue")
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2")
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="obs",add=TRUE,col="blue3")

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo5<-drm(dead/total~dose,weights=total,
                  data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                                myzpbo$pesticide=="thiaclopride_PBO" & 
                                myzpbo$date==levels(myzpbo$date)[5],],
                  fct=LN.2(),
                  type="binomial")
plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red")
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2")
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo4,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo4,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo5,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo5,type="obs",add=TRUE,col="red3")

#the results for the "fourth" date (in fact the third repetition of the test) 
#are not very good. So we remove them for the final analysis
levels(myzpbo$date)[4]
myzpbofinal<-myzpbo[myzpbo$date!=levels(myzpbo$date)[4],]


myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                               myzpbo$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(myzpbo[,2]))[4])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[4] & 
                               myzpbo$pesticide=="thiaclopride_PBO",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")


#another way which is simplier to do the model for every different clone
myzus.sspbofin<-drm(dead/total~dose,weights=total,clone,
                    data=myzpbofinal[myzpbofinal$pesticide=="thiaclopride",],
                    fct=LN.3u(),type="binomial")
plot(myzus.sspbofin)
DL50clonefin<-ED(myzus.sspbofin,50,interval="delta",reference="control")

myzus.avpbofin<-drm(dead/total~dose,weights=total,clone,
                    data=myzpbofinal[myzpbofinal$pesticide!="thiaclopride",],
                    fct=LN.3u(),type="binomial")
plot(myzus.avpbofin,add=TRUE,col="red")
DL50clonepbofin<-ED(myzus.avpbofin,50,interval="delta",reference="control")

