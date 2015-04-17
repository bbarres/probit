###############################################################################
###############################################################################
#determining IC50, ED50, DL50 ...
###############################################################################
###############################################################################

#loading the library
library(drc)

#set the working directory
setwd("~/work/Rfichiers/Githuber/probit_data")


###############################################################################
#example for Venturia inaequalis resistance to Difenoconazol
###############################################################################

#load the dataset
testCI50<-read.table("testCI50.txt",header=T,sep="\t")

#building of the model
tavelure.m1<-drm(mycel~dose,data=testCI50,
                 fct=LL.4(names=c("Slope","Lower Limit","Upper Limit","ED50")))
summary(tavelure.m1)

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

library(plotrix)
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

#this model bound min and max with 0 and 1 respectively
myzus.m1 <- drm(dead/total~dose,weights=total,data=testMyz,fct=LL.2(),
                type="binomial")
plot(myzus.m1,type="confidence")
plot(myzus.m1)

ed50val_myz<-ED(myzus.m1,50,interval="delta",reference="control")

#in order to obtain the same results than with priprobit, the Finney equivalent
#method, we have to remove the constrain the lower bound and chose a log-normal 
#model instead of a log-logistic model
myzus.m1 <- drm(dead/total~dose,weights=total,data=testMyz,fct=LN.3u(),
                type="binomial")
plot(myzus.m1,type="confidence")
plot(myzus.m1)

#the ED50 obtained is identical to the one obtained with priprobit, the SD 
#still differ a little bit
ed50val_myz<-ED(myzus.m1,50,interval="delta",reference="control")


