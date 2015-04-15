###############################################################################
#determining IC50, ED50, DL50 ...
###############################################################################

#loading the library
library(drc)

#set the working directory
setwd("~/work/Rfichiers/Githuber/probit_data")

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
segments(0,50,ed50val[1],50,lty=2,col="red")



