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
#scale. It also include repetition that were flawed because of insufficient
#number of individual for the entire repetition or because there was a problem
#during the lab experiment
myzpbo<-myzpbo[myzpbo$include!="n",]

#let's sum the number of individual tested per dose for each treatment 
#to check the representativeness of the obtained results
checkdat<-aggregate(cbind(dead,total)~dose+pesticide+clone,data=myzpbo,"sum")
write.table(checkdat,file="checkdat.txt",sep="\t",row.names=FALSE)


barplot(cloneX$dead[cloneX$pesticide!="thiaclopride"]/
          cloneX$total[cloneX$pesticide!="thiaclopride"])

###############################################################################
#2015serie
###############################################################################

serie2015<-myzpbo[myzpbo$serie=="2015serie",]
serie2015<-drop.levels(serie2015)

#for the first clone of the series in 2015 11-037-001####
cloneX<-serie2015[serie2015$clone==names(summary(serie2015[,2]))[1],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2015[,2]))[1])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue",pch=19)
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2",pch=19)
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3",pch=19)

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red",pch=19)
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2",pch=19)
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3",pch=19)

#maybe removing 19-Nov-15


#for the second clone of the series in 2015 12-069-001####
cloneX<-serie2015[serie2015$clone==names(summary(serie2015[,2]))[2],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2015[,2]))[2])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue",pch=19)
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2",pch=19)
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3",pch=19)

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red",pch=19)
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2",pch=19)
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3",pch=19)

#maybe removing 19-Nov-15


#for the third clone of the series in 2015 4106A####
cloneX<-serie2015[serie2015$clone==names(summary(serie2015[,2]))[3],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2015[,2]))[3])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
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
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
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

#maybe removing 29-Oct-15


#for the fourth clone of the series in 2015 5191A####
cloneX<-serie2015[serie2015$clone==names(summary(serie2015[,2]))[4],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2015[,2]))[4])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo5<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[5],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo6<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[6],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
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

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo5<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[5],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo6<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[6],],
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






###############################################################################
#2016serie
###############################################################################

serie2016<-myzpbo[myzpbo$serie=="2016serie",]
serie2016<-drop.levels(serie2016)


#for the first clone of the series in 2016 13-001-048####
cloneX<-serie2016[serie2016$clone==names(summary(serie2016[,2]))[1],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2016[,2]))[1])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue",pch=19)
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2",pch=19)
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3",pch=19)

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red",pch=19)
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2",pch=19)
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3",pch=19)


#enlever le 24-Mar-16 ?


#for the second clone of the series in 2016 13-001-032####
cloneX<-serie2016[serie2016$clone==names(summary(serie2016[,2]))[2],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2016[,2]))[2])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
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
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
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





#for the third clone of the series in 2016 384C####
cloneX<-serie2016[serie2016$clone==names(summary(serie2016[,2]))[3],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2016[,2]))[3])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue",pch=19)
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2",pch=19)
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3",pch=19)

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red",pch=19)
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2",pch=19)
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3",pch=19)






#for the fourth clone of the series in 2016 4916A####
cloneX<-serie2016[serie2016$clone==names(summary(serie2016[,2]))[4],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2016[,2]))[4])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
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
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
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


#enlever 21-Jan-16 ?






#for the fifth clone of the series in 2016 5191A####
cloneX<-serie2016[serie2016$clone==names(summary(serie2016[,2]))[5],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2016[,2]))[5])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
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
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
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



#enlever le 21-Jan-16 ?




###############################################################################
#2017serie1
###############################################################################

serie2017s1<-myzpbo[myzpbo$serie=="2017serie1",]
serie2017s1<-drop.levels(serie2017s1)


#for the first clone of the series in 2017 serie 1 11-037-018####
cloneX<-serie2017s1[serie2017s1$clone==names(summary(serie2017s1[,2]))[1],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2017s1[,2]))[1])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo5<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[5],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
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

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo5<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[5],],
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




#for the second clone of the series in 2017 serie 1 11-037-001####
cloneX<-serie2017s1[serie2017s1$clone==names(summary(serie2017s1[,2]))[2],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2017s1[,2]))[2])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo5<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[5],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
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


#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo5<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[5],],
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








#for the third clone of the series in 2017 serie 1 13-001-050####
cloneX<-serie2017s1[serie2017s1$clone==names(summary(serie2017s1[,2]))[3],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2017s1[,2]))[3])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
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
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
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









#for the fourth clone of the series in 2017 serie 1 13-001-032####
cloneX<-serie2017s1[serie2017s1$clone==names(summary(serie2017s1[,2]))[4],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2017s1[,2]))[4])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo5<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[5],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo6<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[6],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo7<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[7],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
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


#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo5<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[5],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo6<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[6],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo7<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[7],],
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



#remove the 28-Mar-17 for 13-001-032
#remove the entire 14-Mar-17 repetition




###############################################################################
#2017serie2
###############################################################################

serie2017s2<-myzpbo[myzpbo$serie=="2017serie2",]
serie2017s2<-drop.levels(serie2017s2)


#for the first clone of the series in 2017 serie 2 11-039-003####
cloneX<-serie2017s2[serie2017s2$clone==names(summary(serie2017s2[,2]))[1],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2017s2[,2]))[1])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue",pch=19)
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2",pch=19)
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3",pch=19)

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red",pch=19)
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2",pch=19)
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3",pch=19)







#for the second clone of the series in 2017 serie 2 11-037-001####
cloneX<-serie2017s2[serie2017s2$clone==names(summary(serie2017s2[,2]))[2],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2017s2[,2]))[2])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue",pch=19)
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2",pch=19)
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3",pch=19)

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red",pch=19)
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2",pch=19)
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3",pch=19)







#for the third clone of the series in 2017 serie 2 13-001-045####
cloneX<-serie2017s2[serie2017s2$clone==names(summary(serie2017s2[,2]))[3],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2017s2[,2]))[3])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue",pch=19)
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2",pch=19)
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3",pch=19)

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red",pch=19)
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2",pch=19)
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3",pch=19)


















#for the fourth clone of the series in 2017 serie 2 13-001-041####
cloneX<-serie2017s2[serie2017s2$clone==names(summary(serie2017s2[,2]))[4],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2017s2[,2]))[4])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue",pch=19)
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2",pch=19)
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3",pch=19)

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red",pch=19)
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2",pch=19)
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3",pch=19)




###############################################################################
#2017serie3
###############################################################################

serie2017s3<-myzpbo[myzpbo$serie=="2017serie3",]
serie2017s3<-drop.levels(serie2017s3)


#for the first clone of the series in 2017 serie 3 11-060-006####
cloneX<-serie2017s3[serie2017s3$clone==names(summary(serie2017s3[,2]))[1],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2017s3[,2]))[1])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue",pch=19)
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2",pch=19)
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3",pch=19)

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red",pch=19)
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2",pch=19)
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3",pch=19)






#for the second clone of the series in 2017 serie 3 11-037-012####
cloneX<-serie2017s3[serie2017s3$clone==names(summary(serie2017s3[,2]))[2],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2017s3[,2]))[2])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
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
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
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




#for the third clone of the series in 2017 serie 3 11-037-001####
cloneX<-serie2017s3[serie2017s3$clone==names(summary(serie2017s3[,2]))[3],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2017s3[,2]))[3])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
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
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo4<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[4],],
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









#for the fourth clone of the series in 2017 serie 3 12-068-006####
cloneX<-serie2017s3[serie2017s3$clone==names(summary(serie2017s3[,2]))[4],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2017s3[,2]))[4])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue",pch=19)
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2",pch=19)
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3",pch=19)

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red",pch=19)
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2",pch=19)
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3",pch=19)

















#for the fifth clone of the series in 2017 serie 3 12-067-023####
cloneX<-serie2017s3[serie2017s3$clone==names(summary(serie2017s3[,2]))[5],]
cloneX<-drop.levels(cloneX)
myzus.sspbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide=="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.sspbo,type="confidence",main=names(summary(serie2017s3[,2]))[5])
plot(myzus.sspbo,type="obs",add=TRUE)

myzus.avpbo<-drm(dead/total~dose,weights=total,
                 data=cloneX[cloneX$pesticide!="thiaclopride",],
                 fct=LN.2(),
                 type="binomial")
plot(myzus.avpbo,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo,type="obs",add=TRUE,col="red")

#curves for each repetition without PBO
myzus.sspbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.sspbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide=="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.sspbo,type="obs",main=names(summary(cloneX[,2]))[1])
plot(myzus.sspbo1,type="confidence",add=TRUE,col="blue")
plot(myzus.sspbo1,type="obs",add=TRUE,col="blue",pch=19)
plot(myzus.sspbo2,type="confidence",add=TRUE,col="blue2")
plot(myzus.sspbo2,type="obs",add=TRUE,col="blue2",pch=19)
plot(myzus.sspbo3,type="confidence",add=TRUE,col="blue3")
plot(myzus.sspbo3,type="obs",add=TRUE,col="blue3",pch=19)

#curves for each repetition with PBO
myzus.avpbo1<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[1],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo2<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[2],],
                  fct=LN.2(),
                  type="binomial")
myzus.avpbo3<-drm(dead/total~dose,weights=total,
                  data=cloneX[cloneX$pesticide!="thiaclopride" & 
                                cloneX$date==levels(cloneX$date)[3],],
                  fct=LN.2(),
                  type="binomial")

plot(myzus.avpbo1,type="confidence",add=TRUE,col="red")
plot(myzus.avpbo1,type="obs",add=TRUE,col="red",pch=19)
plot(myzus.avpbo2,type="confidence",add=TRUE,col="red2")
plot(myzus.avpbo2,type="obs",add=TRUE,col="red2",pch=19)
plot(myzus.avpbo3,type="confidence",add=TRUE,col="red3")
plot(myzus.avpbo3,type="obs",add=TRUE,col="red3",pch=19)
