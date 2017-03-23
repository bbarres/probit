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
myzpbo<-read.table("pbo_140317.txt",header=T,sep="\t")

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
#let's evaluate the different model for each clone with and without pbo
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
plot(myzus.sspbo,type="obs",main=names(summary(myzpbo[,2]))[1])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue")
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2")
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="obs",add=TRUE,col="blue3")

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
plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red")
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2")
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo4,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo4,type="obs",add=TRUE,col="red3")


#let's evaluate the different model for each clone with and without pbo
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
plot(myzus.sspbo,type="obs",main=names(summary(myzpbo[,2]))[2])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue")
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2")
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="obs",add=TRUE,col="blue3")

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
plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red")
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2")
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo4,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo4,type="obs",add=TRUE,col="red3")

#let's evaluate the different model for each clone with and without pbo
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
plot(myzus.sspbo,type="obs",main=names(summary(myzpbo[,2]))[3])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue")
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2")
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="obs",add=TRUE,col="blue3")

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
plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red")
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2")
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo4,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo4,type="obs",add=TRUE,col="red3")

#let's evaluate the different model for each clone with and without pbo
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
plot(myzus.sspbo,type="obs",main=names(summary(myzpbo[,2]))[4])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue")
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2")
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="obs",add=TRUE,col="blue3")

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
plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red")
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2")
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3")
plot(myzus.avpbo4,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo4,type="obs",add=TRUE,col="red3")


myzpbo<-cbind(myzpbo,"perc"=myzpbo$dead*100/myzpbo$total)
myzpbo<-cbind(myzpbo,"testID"=paste(myzpbo$clone,myzpbo$date))
#concatenation of the data
agraa<-aggregate(cbind(dead,total)~dose+pesticide+clone,data=myzpbo,"sum")
agraa<-cbind(agraa,"perc"=agraa$dead*100/agraa$total)
barplot(as.matrix(agraa[1:11,c("dead","total","perc")]),beside=TRUE)
plot(agraa$dead/agraa$total)

recapall<-myzpbo[myzpbo$dose==0 & myzpbo$pesticide=="thiaclopride",
                 c("testID","perc")]
recapall<-cbind(recapall,"total"=aggregate(cbind(dead,total)~pesticide+clone+date,
          data=myzpbo[myzpbo$pesticide=="thiaclopride",],"sum")$total)
recapall<-cbind(recapall,"percpbo"=myzpbo[myzpbo$dose==0 & 
                                  myzpbo$pesticide!="thiaclopride",
                                c("perc")])
recapall<-cbind(recapall,"totalpbo"=aggregate(cbind(dead,total)~pesticide+clone+date,
          data=myzpbo[myzpbo$pesticide!="thiaclopride",],"sum")$total)


#let's evaluate the different model for each clone with and without pbo
DL50<-c()
for (i in 1:16){
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=myzpbo[myzpbo$testID==recapall$testID[i] & 
                               myzpbo$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
ed50val<-ED(myzus.sspbo,50,interval="delta",reference="control")
DL50<-c(DL50,ed50val[1])
}

DL50pbo<-c()
for (i in 1:13){
  myzus.sspbo<-drm(dead/total~dose,weights=total,
                   data=myzpbo[myzpbo$testID==recapall$testID[i] & 
                                 myzpbo$pesticide!="thiaclopride",],
                   fct=LN.2(),
                   type="binomial")
  ed50val<-ED(myzus.sspbo,50,interval="delta",reference="control")
  DL50pbo<-c(DL50pbo,ed50val[1])
}
DL50pbo<-c(DL50pbo,NA)
for (i in 15:16){
  myzus.sspbo<-drm(dead/total~dose,weights=total,
                   data=myzpbo[myzpbo$testID==recapall$testID[i] & 
                                 myzpbo$pesticide!="thiaclopride",],
                   fct=LN.2(),
                   type="binomial")
  ed50val<-ED(myzus.sspbo,50,interval="delta",reference="control")
  DL50pbo<-c(DL50pbo,ed50val[1])
}

recapall<-cbind(recapall,"DL50"=DL50,"DL50pbo"=DL50pbo)
recapall<-cbind(recapall,"FI"=recapall$DL50/recapall$DL50pbo)
recapall[order(recapall$testID),]
write.table(recapall[order(recapall$testID),], file="tablerecap.txt",
            sep="\t",quote=FALSE,row.names=FALSE)

par(op)
