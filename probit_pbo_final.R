###############################################################################
###############################################################################
#code for the PBO tests ...
###############################################################################
###############################################################################

#loading the libraries
library(drc)
library(plotrix)
library(gdata)

#set the working directory
setwd("~/work/Rfichiers/Githuber/probit_data")


###############################################################################
#comparison between repetitions
###############################################################################

#load the dataset
myzpbo<-read.table("final_av_ss_pbo.txt",header=T,sep="\t")
#because some concentration were only used for adapting the pesticide dose
#scale
myzpbo<-myzpbo[myzpbo$include!="n",]

#let's sum the number of individual tested per dose for each treatment 
#to check the representativeness of the obtained results
checkdat<-aggregate(cbind(dead,total)~dose+pesticide+clone,data=myzpbo,"sum")
write.table(checkdat,file="checkdat.txt",sep="\t",row.names=FALSE)



#for the first clone of the series in 2015, 2016####
clone1<-myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[8],]
clone1<-drop.levels(clone1)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=clone1[clone1$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(myzpbo[,2]))[8])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=clone1[clone1$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide=="thiaclopride" & 
                                clone1$date==levels(clone1$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide=="thiaclopride" & 
                                clone1$date==levels(clone1$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide=="thiaclopride" & 
                                clone1$date==levels(clone1$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(clone1[,2]))[1])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue",pch=19)
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2",pch=19)
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3",pch=19)

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide!="thiaclopride" & 
                                clone1$date==levels(clone1$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide!="thiaclopride" & 
                                clone1$date==levels(clone1$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide!="thiaclopride" & 
                                clone1$date==levels(clone1$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red",pch=19)
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2",pch=19)
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3",pch=19)






#for the second clone of the series in 2015, 2016####
clone1<-myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[15],]
clone1<-drop.levels(clone1)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=clone1[clone1$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(myzpbo[,2]))[15])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=clone1[clone1$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")


#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide=="thiaclopride" & 
                                clone1$date==levels(clone1$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide=="thiaclopride" & 
                                clone1$date==levels(clone1$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide=="thiaclopride" & 
                                clone1$date==levels(clone1$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide=="thiaclopride" & 
                                clone1$date==levels(clone1$date)[4],],
                  fct=LN.2(),
                  type="binomial")


plot(myzus.sspbo,type="obs",main=names(summary(clone1[,2]))[1])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue",pch=19)
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2",pch=19)
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3",pch=19)
plot(myzus.sspbo4,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="obs",add=TRUE,col="blue3",pch=19)


#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide!="thiaclopride" & 
                                clone1$date==levels(clone1$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide!="thiaclopride" & 
                                clone1$date==levels(clone1$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide!="thiaclopride" & 
                                clone1$date==levels(clone1$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide!="thiaclopride" & 
                                clone1$date==levels(clone1$date)[4],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red",pch=19)
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2",pch=19)
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3",pch=19)
plot(myzus.avpbo4,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo4,type="obs",add=TRUE,col="red3",pch=19)






#for the second clone of the series in 2015, 2016####
clone1<-myzpbo[myzpbo$clone==names(summary(myzpbo[,2]))[17],]
clone1<-drop.levels(clone1)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=clone1[clone1$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(myzpbo[,2]))[17])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=clone1[clone1$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide=="thiaclopride" & 
                                clone1$date==levels(clone1$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide=="thiaclopride" & 
                                clone1$date==levels(clone1$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide=="thiaclopride" & 
                                clone1$date==levels(clone1$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide=="thiaclopride" & 
                                clone1$date==levels(clone1$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo5<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide=="thiaclopride" & 
                                clone1$date==levels(clone1$date)[5],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo6<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide=="thiaclopride" & 
                                clone1$date==levels(clone1$date)[6],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo7<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide=="thiaclopride" & 
                                clone1$date==levels(clone1$date)[7],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo8<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide=="thiaclopride" & 
                                clone1$date==levels(clone1$date)[8],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo9<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide=="thiaclopride" & 
                                clone1$date==levels(clone1$date)[9],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo10<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide=="thiaclopride" & 
                                clone1$date==levels(clone1$date)[10],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(clone1[,2]))[1])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue",pch=19)
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2",pch=19)
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3",pch=19)
plot(myzus.sspbo4,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo4,type="obs",add=TRUE,col="blue3",pch=19)
plot(myzus.sspbo5,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo5,type="obs",add=TRUE,col="blue3",pch=19)
plot(myzus.sspbo6,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo6,type="obs",add=TRUE,col="blue3",pch=19)
plot(myzus.sspbo7,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo7,type="obs",add=TRUE,col="blue3",pch=19)
plot(myzus.sspbo8,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo8,type="obs",add=TRUE,col="blue3",pch=19)
plot(myzus.sspbo9,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo9,type="obs",add=TRUE,col="blue3",pch=19)
plot(myzus.sspbo10,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo10,type="obs",add=TRUE,col="blue3",pch=19)

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide!="thiaclopride" & 
                                clone1$date==levels(clone1$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide!="thiaclopride" & 
                                clone1$date==levels(clone1$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide!="thiaclopride" & 
                                clone1$date==levels(clone1$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide!="thiaclopride" & 
                                clone1$date==levels(clone1$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo5<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide!="thiaclopride" & 
                                clone1$date==levels(clone1$date)[5],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo6<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide!="thiaclopride" & 
                                clone1$date==levels(clone1$date)[6],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo7<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide!="thiaclopride" & 
                                clone1$date==levels(clone1$date)[7],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo8<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide!="thiaclopride" & 
                                clone1$date==levels(clone1$date)[8],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo9<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide!="thiaclopride" & 
                                clone1$date==levels(clone1$date)[9],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo10<-drm(dead/total~dose,weights=total,
                  data=clone1[clone1$pesticide!="thiaclopride" & 
                                clone1$date==levels(clone1$date)[10],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red",pch=19)
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2",pch=19)
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3",pch=19)
plot(myzus.avpbo4,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo4,type="obs",add=TRUE,col="red3",pch=19)
plot(myzus.avpbo5,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo5,type="obs",add=TRUE,col="red3",pch=19)
plot(myzus.avpbo6,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo6,type="obs",add=TRUE,col="red3",pch=19)
plot(myzus.avpbo7,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo7,type="obs",add=TRUE,col="red3",pch=19)
plot(myzus.avpbo8,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo8,type="obs",add=TRUE,col="red3",pch=19)
plot(myzus.avpbo9,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo9,type="obs",add=TRUE,col="red3",pch=19)
plot(myzus.avpbo10,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo10,type="obs",add=TRUE,col="red3",pch=19)



