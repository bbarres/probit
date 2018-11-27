###############################################################################
###############################################################################
#Venturia inaequalis Monitoring data analysis
###############################################################################
###############################################################################

#loading the libraries
library(drc)
library(plotrix)
library(gdata)

#load the global dataset
ventu18<-read.table(file="data/18-venturia.txt",header=T,sep=";")
# ventu18$bioagr_ID<-factor(ventu18$bioagr_ID,
#                           levels=rev(levels(ventu18$bioagr_ID)))
collist<-c("forestgreen","black")
levels(ventu18$pest_sa_id)


###############################################################################
#Analysis for the boscalid
###############################################################################

#subsetting the global dataset
databysa<-ventu18[ventu18$pest_sa_id=="BOSCALID" & 
                    ventu18$lect_echec!=1,]

#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
rez<-as.character(databysa[databysa$dose=="30" & databysa$rslt_03>50,
                           "ech_id"])
REZ<-data.frame("ech_id"=as.character(),"ED50"=as.character())
#we limit the dataset to the sample that reach somehow a IC of 50%
#databysa<-databysa[!(databysa$ech_id %in% rez),]
databysa<-drop.levels(databysa,reorder=FALSE)
pdf(file="output/boscalid.pdf",width=20,height=12)
op<-par(mfrow=c(3,6))
for (i in 1: dim(table(databysa$ech_id))[1]) {
  datatemp<-databysa[databysa$ech_id==names(table(databysa$ech_id))[i],]
  couleur<-collist[as.numeric(datatemp$bioagr_id)]
  typeline<-as.numeric(datatemp$bioagr_id)[1]
  if (is.na(datatemp[1,"rslt_03"])==TRUE) {
    tempx<-data.frame("ech_id"=names(table(databysa$ech_id))[i],
                      "ED50"=c("NA"))
    REZ<-rbind(REZ,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(databysa$ech_id))[i])
  } else {
    temp.m1<-drm(rslt_03~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(databysa$ech_id))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
    temp<-ED(temp.m1,50,type="absolute")
    tempx<-data.frame("ech_id"=names(table(databysa$ech_id))[i],
                      "ED50"=temp[1])
    REZ<-rbind(REZ,tempx)
  }
}
par(op)
dev.off()

REZbos<-REZ


###############################################################################
#Analysis for the captane
###############################################################################

#subsetting the global dataset
databysa<-ventu18[ventu18$pest_sa_id=="CAPTANE" & 
                    ventu18$lect_echec!=1,]

#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
rez<-as.character(databysa[databysa$dose=="30" & databysa$rslt_03>50,
                           "ech_id"])
REZ<-data.frame("ech_id"=as.character(),"ED50"=as.character())
#we limit the dataset to the sample that reach somehow a IC of 50%
#databysa<-databysa[!(databysa$ech_id %in% rez),]
databysa<-drop.levels(databysa,reorder=FALSE)
pdf(file="output/captane.pdf",width=20,height=12)
op<-par(mfrow=c(3,6))
for (i in 1: dim(table(databysa$ech_id))[1]) {
  datatemp<-databysa[databysa$ech_id==names(table(databysa$ech_id))[i],]
  couleur<-collist[as.numeric(datatemp$bioagr_id)]
  typeline<-as.numeric(datatemp$bioagr_id)[1]
  if (is.na(datatemp[1,"rslt_03"])==TRUE) {
    tempx<-data.frame("ech_id"=names(table(databysa$ech_id))[i],
                      "ED50"=c("NA"))
    REZ<-rbind(REZ,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(databysa$ech_id))[i])
  } else {
    temp.m1<-drm(rslt_03~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(databysa$ech_id))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
    temp<-ED(temp.m1,50,type="absolute")
    tempx<-data.frame("ech_id"=names(table(databysa$ech_id))[i],
                      "ED50"=temp[1])
    REZ<-rbind(REZ,tempx)
  }
}
par(op)
dev.off()

REZcap<-REZ


###############################################################################
#Analysis for the cyprodinil
###############################################################################

#subsetting the global dataset
databysa<-ventu18[ventu18$pest_sa_id=="CYPRODINIL" & 
                    ventu18$lect_echec!=1,]

#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
rez<-as.character(databysa[databysa$dose=="30" & databysa$rslt_03>50,
                           "ech_id"])
REZ<-data.frame("ech_id"=as.character(),"ED50"=as.character())
#we limit the dataset to the sample that reach somehow a IC of 50%
#databysa<-databysa[!(databysa$ech_id %in% rez),]
databysa<-drop.levels(databysa,reorder=FALSE)
pdf(file="output/cyprodinil.pdf",width=20,height=12)
op<-par(mfrow=c(3,6))
for (i in 1: dim(table(databysa$ech_id))[1]) {
  datatemp<-databysa[databysa$ech_id==names(table(databysa$ech_id))[i],]
  couleur<-collist[as.numeric(datatemp$bioagr_id)]
  typeline<-as.numeric(datatemp$bioagr_id)[1]
  if (is.na(datatemp[1,"rslt_03"])==TRUE) {
    tempx<-data.frame("ech_id"=names(table(databysa$ech_id))[i],
                      "ED50"=c("NA"))
    REZ<-rbind(REZ,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(databysa$ech_id))[i])
  } else {
    temp.m1<-drm(rslt_03~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(databysa$ech_id))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
    temp<-ED(temp.m1,50,type="absolute")
    tempx<-data.frame("ech_id"=names(table(databysa$ech_id))[i],
                      "ED50"=temp[1])
    REZ<-rbind(REZ,tempx)
  }
}
par(op)
dev.off()

REZcyp<-REZ


###############################################################################
#Analysis for the dithianon
###############################################################################

#subsetting the global dataset
databysa<-ventu18[ventu18$pest_sa_id=="DITHIANON" & 
                    ventu18$lect_echec!=1,]

#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
rez<-as.character(databysa[databysa$dose=="30" & databysa$rslt_03>50,
                           "ech_id"])
REZ<-data.frame("ech_id"=as.character(),"ED50"=as.character())
#we limit the dataset to the sample that reach somehow a IC of 50%
#databysa<-databysa[!(databysa$ech_id %in% rez),]
databysa<-drop.levels(databysa,reorder=FALSE)
pdf(file="output/dithianon.pdf",width=20,height=12)
op<-par(mfrow=c(3,6))
for (i in 1: dim(table(databysa$ech_id))[1]) {
  datatemp<-databysa[databysa$ech_id==names(table(databysa$ech_id))[i],]
  couleur<-collist[as.numeric(datatemp$bioagr_id)]
  typeline<-as.numeric(datatemp$bioagr_id)[1]
  if (is.na(datatemp[1,"rslt_03"])==TRUE) {
    tempx<-data.frame("ech_id"=names(table(databysa$ech_id))[i],
                      "ED50"=c("NA"))
    REZ<-rbind(REZ,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(databysa$ech_id))[i])
  } else {
    temp.m1<-drm(rslt_03~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(databysa$ech_id))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
    temp<-ED(temp.m1,50,type="absolute")
    tempx<-data.frame("ech_id"=names(table(databysa$ech_id))[i],
                      "ED50"=temp[1])
    REZ<-rbind(REZ,tempx)
  }
}
par(op)
dev.off()

REZdit<-REZ


###############################################################################
#Analysis for the dodine
###############################################################################

#subsetting the global dataset
databysa<-ventu18[ventu18$pest_sa_id=="DODINE" & 
                    ventu18$lect_echec!=1,]

#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
rez<-as.character(databysa[databysa$dose=="30" & databysa$rslt_03>50,
                           "ech_id"])
REZ<-data.frame("ech_id"=as.character(),"ED50"=as.character())
#we limit the dataset to the sample that reach somehow a IC of 50%
#databysa<-databysa[!(databysa$ech_id %in% rez),]
databysa<-drop.levels(databysa,reorder=FALSE)
pdf(file="output/dodine.pdf",width=20,height=12)
op<-par(mfrow=c(3,6))
for (i in 1: dim(table(databysa$ech_id))[1]) {
  datatemp<-databysa[databysa$ech_id==names(table(databysa$ech_id))[i],]
  couleur<-collist[as.numeric(datatemp$bioagr_id)]
  typeline<-as.numeric(datatemp$bioagr_id)[1]
  if (is.na(datatemp[1,"rslt_03"])==TRUE) {
    tempx<-data.frame("ech_id"=names(table(databysa$ech_id))[i],
                      "ED50"=c("NA"))
    REZ<-rbind(REZ,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(databysa$ech_id))[i])
  } else {
    temp.m1<-drm(rslt_03~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(databysa$ech_id))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
    temp<-ED(temp.m1,50,type="absolute")
    tempx<-data.frame("ech_id"=names(table(databysa$ech_id))[i],
                      "ED50"=temp[1])
    REZ<-rbind(REZ,tempx)
  }
}
par(op)
dev.off()

REZdod<-REZ


###############################################################################
#combined results and export
###############################################################################

REZ<-rbind(REZbos,REZcap,REZcyp,REZdit,REZdod)
REZ<-cbind(REZ,"active_substance"=c(rep("boscalid",18),
                                    rep("captane",16),
                                    rep("cyprodinile",3),
                                    rep("dithianon",16),
                                    rep("dodine",17)))

write.table(REZ,file="output/result_CI50.txt",quote=FALSE,col.names=TRUE, 
            row.names=FALSE,sep="\t")


###############################################################################
#END
###############################################################################
