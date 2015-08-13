###############################################################################
###############################################################################
#determining IC50, ED50, DL50 ...
###############################################################################
###############################################################################

#loading the libraries
library(drc)
library(plotrix)

#set the working directory
setwd("~/work/Rfichiers/Githuber/probit_data")


###############################################################################
#example for Venturia inaequalis resistance to Difenoconazol
###############################################################################

#load the dataset
testCI50<-read.table("testCI50.txt",header=T,sep="\t")

#building of the model
tavelure.m1<-drm((mycel-7.5)~dose,data=testCI50,
                 fct=LL.4(names=c("Slope","Lower Limit","Upper Limit","ED50")))
summary(tavelure.m1)
plot(tavelure.m1)
tavelure.m2<-drm((mycel-7.5)~dose,data=testCI50,fct=LL.3())
summary(tavelure.m2)
plot(tavelure.m2)

#comparison of different model
mselect(tavelure.m1,list(LL.3(),LN.3(),LN.4()))

#plot with 95% confidence interval
plot(tavelure.m1,broken=TRUE,type="confidence")
plot(tavelure.m1,broken=TRUE,add=TRUE)
#a simplier plot
plot(tavelure.m1,broken=TRUE)

#evaluating the ED50
ed50val<-ED(tavelure.m1,50,interval="delta")
#this is also working for a list of value for ED10, ED90...
ed_val<-ED(tavelure.m1,c(10,50,90),interval="delta")

#predict the value of the measured trait for the ED50
predict(tavelure.m1,data.frame(dose=ed50val[1]),se.fit=TRUE)
#values predicted for the set of dose used in the test
predict(tavelure.m1, interval = "confidence")

#different possible plot for the results
plot(tavelure.m1,broken=TRUE)
abline(v=ed50val[1],col="red",lty=2)
abline(v=ed50val[3],col="red",lty=3)
abline(v=ed50val[4],col="red",lty=3)

#still doesn't work
plot(tavelure.m1,broken=TRUE)
segments(0,50,ed50val[1],50,lty=2,col="red")
#the not very good solution (to say the least)
segments(0.0001,predict(tavelure.m1,data.frame(dose=ed50val[1])),
         0.001,predict(tavelure.m1,data.frame(dose=ed50val[1])),
         lty=2,col="red")
segments(0.001,predict(tavelure.m1,data.frame(dose=ed50val[1])),
         0.01,predict(tavelure.m1,data.frame(dose=ed50val[1])),
         lty=2,col="red")
segments(0.01,predict(tavelure.m1,data.frame(dose=ed50val[1])),
         0.1,predict(tavelure.m1,data.frame(dose=ed50val[1])),
         lty=2,col="red")
segments(0.1,predict(tavelure.m1,data.frame(dose=ed50val[1])),
         ed50val[1],predict(tavelure.m1,data.frame(dose=ed50val[1])),
         lty=2,col="red")
segments(ed50val[1],0,
         ed50val[1],predict(tavelure.m1,data.frame(dose=ed50val[1])),
         lty=2,col="red")


###############################################################################
#example for Myzus persicae resistance to Difenoconazol
###############################################################################

#load the dataset
testMyz<-read.table("imida_cum_11037_18ter.txt",header=T,sep="\t")

#here is the same model described in the help of the R package, both
#the min and max are constrained with 0 and 1, respectively
myzus.m1<-drm(dead/total~dose,weights=total,data=testMyz,fct=LL.2(),
              type="binomial")
plot(myzus.m1)
ed50val_myz1<-ED(myzus.m1,50,interval="delta")

#this model constrains only the max (upperlimit) with 1. Probably 
#the most elegant model for THIS particular dataset
myzus.m2 <- drm(dead/total~dose,weights=total,data=testMyz,fct=LL.3u(),
                type="binomial")
plot(myzus.m2,type="confidence")
plot(myzus.m2)
ed50val_myz2<-ED(myzus.m2,50,interval="delta",reference="control")

#in order to obtain the same results than with priprobit, the Finney equivalent
#method, we have to remove the constrain on lowerlimit and chose a log-normal 
#model instead of a log-logistic model
myzus.m3<-drm(dead/total~dose,weights=total,data=testMyz,fct=LN.3u(),
              type="binomial")
plot(myzus.m3)
#the ED50 obtained is identical to the one obtained with priprobit, the SD 
#still differ a little bit
ed50val_myz3<-ED(myzus.m3,50,interval="delta",reference="control")

#another model could the model with all the 4 parameters unconstrained. 
#Probably the most ubiquituous model that will adapt to most datasets
myzus.m4<-drm(dead/total~dose,weights=total,data=testMyz,fct=LL.4(),
              type="binomial")
plot(myzus.m4)
ed50val_myz4<-ED(myzus.m4,50,interval="delta")

#graphical comparison of the four different models
plot(myzus.m1,col="green",lty=2)
plot(myzus.m2,add=TRUE,col="red",lty=2)
plot(myzus.m3,add=TRUE,col="blue",lty=2)
plot(myzus.m4,add=TRUE)

#comparison and selection of the model, the smaller the IC, the better, 
#the higher the lack of fit test, the better
mselect(myzus.m1,list(LL.3u(),LN.3u(),LL.4(),LN.4(),W1.3u(),W2.3u()))
#here the best model is LL.3u


###############################################################################
#Test on the leafhopper (Scaphoïdeus titanus) data of Jérémy
###############################################################################

#load the first dataset
ScaTitc<-read.table("0c.txt",header=T,sep="\t")

#here is the same model described in the help of the R package, both
#the min and max are constrained with 0 and 1, respectively
ScaTitc.m1<-drm(dead/total~dose,weights=total,data=ScaTitc,fct=LL.2(),
              type="binomial")
plot(ScaTitc.m1,type="confidence")
plot(ScaTitc.m1, main="constraint 0/1")
ed50val_ScaTitc1<-ED(ScaTitc.m1,50,interval="delta")

#this model constrains only the max (upperlimit) with 1. 
ScaTitc.m2 <- drm(dead/total~dose,weights=total,data=ScaTitc,fct=LL.3u(),
                type="binomial")
plot(ScaTitc.m2,type="confidence")
plot(ScaTitc.m2)
ed50val_ScaTitc2<-ED(ScaTitc.m2,50,interval="delta",reference="control")

#in order to obtain the same results than with priprobit, the Finney equivalent
#method, we have to remove the constrain on lowerlimit and chose a log-normal 
#model instead of a log-logistic model
ScaTitc.m3<-drm(dead/total~dose,weights=total,data=ScaTitc,fct=LN.3u(),
              type="binomial")
plot(ScaTitc.m3,type="confidence")
plot(ScaTitc.m3)
#the ED50 obtained is identical to the one obtained with priprobit, the SD 
#still differ a little bit
ed50val_ScaTitc3<-ED(ScaTitc.m3,50,interval="delta",reference="control")

#another model could the model with all the 4 parameters unconstrained. 
#Probably the most ubiquituous model that will adapt to most datasets
ScaTitc.m4<-drm(dead/total~dose,weights=total,data=ScaTitc,fct=LL.4(),
              type="binomial")
plot(ScaTitc.m4,type="confidence")
plot(ScaTitc.m4)
ed50val_ScaTitc4<-ED(ScaTitc.m4,50,interval="delta")

#graphical comparison of the four different models
plot(ScaTitc.m1,col="green",lty=2)
plot(ScaTitc.m2,add=TRUE,col="red",lty=2)
plot(ScaTitc.m3,add=TRUE,col="blue",lty=2)
plot(ScaTitc.m4,add=TRUE)

#comparison and selection of the model, the smaller the IC, the better, 
#the higher the lack of fit test, the better
mselect(ScaTitc.m1,list(LL.3u(),LN.3u(),LL.4(),LN.4(),W1.3u(),W2.3u()))
#here the best model is LL.3u



#load the second dataset
ScaTitc<-read.table("0d.txt",header=T,sep="\t")

#here is the same model described in the help of the R package, both
#the min and max are constrained with 0 and 1, respectively
ScaTitc.m1<-drm(dead/total~dose,weights=total,data=ScaTitc,fct=LL.2(),
                type="binomial")
plot(ScaTitc.m1,type="confidence")
plot(ScaTitc.m1, main="constraint 0/1")
ed50val_ScaTitc1<-ED(ScaTitc.m1,50,interval="delta")

#this model constrains only the max (upperlimit) with 1. 
ScaTitc.m2 <- drm(dead/total~dose,weights=total,data=ScaTitc,fct=LL.3u(),
                  type="binomial")
plot(ScaTitc.m2,type="confidence")
plot(ScaTitc.m2)
ed50val_ScaTitc2<-ED(ScaTitc.m2,50,interval="delta",reference="control")

#in order to obtain the same results than with priprobit, the Finney equivalent
#method, we have to remove the constrain on lowerlimit and chose a log-normal 
#model instead of a log-logistic model
ScaTitc.m3<-drm(dead/total~dose,weights=total,data=ScaTitc,fct=LN.3u(),
                type="binomial")
plot(ScaTitc.m3,type="confidence")
plot(ScaTitc.m3)
#the ED50 obtained is identical to the one obtained with priprobit, the SD 
#still differ a little bit
ed50val_ScaTitc3<-ED(ScaTitc.m3,50,interval="delta",reference="control")

#another model could the model with all the 4 parameters unconstrained. 
#Probably the most ubiquituous model that will adapt to most datasets
ScaTitc.m4<-drm(dead/total~dose,weights=total,data=ScaTitc,fct=LL.4(),
                type="binomial")
plot(ScaTitc.m4,type="confidence")
plot(ScaTitc.m4)
ed50val_ScaTitc4<-ED(ScaTitc.m4,50,interval="delta")

#graphical comparison of the four different models
plot(ScaTitc.m1,col="green",lty=2)
plot(ScaTitc.m2,add=TRUE,col="red",lty=2)
plot(ScaTitc.m3,add=TRUE,col="blue",lty=2)
plot(ScaTitc.m4,add=TRUE)

#comparison and selection of the model, the smaller the IC, the better, 
#the higher the lack of fit test, the better
mselect(ScaTitc.m1,list(LL.3u(),LN.3u(),LL.4(),LN.4(),W1.3u(),W2.3u()))
#here the best model is LL.3u


###############################################################################
#code for the figure of CLaire's article
###############################################################################

#load the dataset
ImidaMyz<-read.table("imidata.txt",header=T,sep="\t")

#let's isolate the name of the clone and their genotype for the R81T resistance
clone_gen<-ImidaMyz[ImidaMyz$dose==0, 1:2]

#we change the coding for the genotype in order to have the color 
#green, orange and red associated to the RR, RT and TT genotypes
levels(ImidaMyz$Rgeno)<-c("green3","darkorange","firebrick3")

#in order to be consistent with the content of the paper, the default model 
#used for all the different clones will be 'LN.3u()'. This is the equivalent
#to the Finney method. There is no constrain on the lowerlimit and chose a 
#log-normal link is chosen. 
temp<-ImidaMyz[ImidaMyz$ind_ID==clone_gen[1,1] & ImidaMyz$tota!=0,]
temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.3u(),
              type="binomial")

plot(temp.mod,xlim=c(0,max(ImidaMyz$dose)+30000),type="obs",broken=TRUE)
plot(temp.mod,xlim=c(0,100000),type="obs",broken=TRUE)
plot(temp.mod,xlim=c(0,100000),add=TRUE,lwd=2.5,
     type="none",col=as.character(temp$Rgeno[1]))

for (i in 2:dim(clone_gen)[1]) {
  temp<-ImidaMyz[ImidaMyz$ind_ID==clone_gen[i,1] & ImidaMyz$tota!=0,]
  temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.3u(),
                type="binomial")
  plot(temp.mod,xlim=c(0,100000),type="obs",add=TRUE)
  plot(temp.mod,xlim=c(0,100000),add=TRUE,lwd=2.5,
       type="none",col=as.character(temp$Rgeno[1]))
}

#export pdf 15 x 8 inches

###############################################################################
#END
###############################################################################