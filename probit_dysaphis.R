##############################################################################/
##############################################################################/
#code for the dysaphis tests ...
##############################################################################/
##############################################################################/

#loading the libraries
library(drc)
library(plotrix)


#load the dataset
dysaphis<-read.table("data/dysaphis.txt",header=T,sep="\t")

#reformating the data by "merging" the different wells of the same repetition
dysaphis<-cbind(aggregate(dead~dose+date+clone+pesticide,data=dysaphis,"sum"), 
                "total"=aggregate(total~dose+date+clone+pesticide,
                                  data=dysaphis,"sum")$total)[,c(2,3,1,5,6,4)]


##############################################################################/
#First attempt to estimate the LD50 for Dysaphis plantaginae vs flonicamid####
##############################################################################/

dysa_mod1<-drm(dead/total~dose,weights=total,
               data=dysaphis[dysaphis$clone==names(summary(dysaphis[,2]))[1] & 
                             dysaphis$date==levels(dysaphis$date)[1],],
                  fct=LN.3u(),
                  type="binomial")

dysa_mod2<-drm(dead/total~dose,weights=total,
               data=dysaphis[dysaphis$clone==names(summary(dysaphis[,2]))[1] & 
                               dysaphis$date==levels(dysaphis$date)[2],],
               fct=LN.3u(),
               type="binomial")
dysa_mod3<-drm(dead/total~dose,weights=total,
               data=dysaphis[dysaphis$clone==names(summary(dysaphis[,2]))[1] & 
                               dysaphis$date==levels(dysaphis$date)[3],],
               fct=LN.3u(),
               type="binomial")

plot(dysa_mod1,type="confidence",main=names(summary(dysaphis[,2]))[1])
plot(dysa_mod1,type="obs",add=TRUE)
plot(dysa_mod2,type="confidence",col="red",add=TRUE)
plot(dysa_mod3,type="confidence",col="blue",add=TRUE)


dysa_mod1<-drm(dead/total~dose,weights=total,
               data=dysaphis[dysaphis$clone==names(summary(dysaphis[,2]))[2] & 
                               dysaphis$date==levels(dysaphis$date)[1],],
               fct=LN.3u(),
               type="binomial")

dysa_mod2<-drm(dead/total~dose,weights=total,
               data=dysaphis[dysaphis$clone==names(summary(dysaphis[,2]))[2] & 
                               dysaphis$date==levels(dysaphis$date)[2],],
               fct=LN.3u(),
               type="binomial")
dysa_mod3<-drm(dead/total~dose,weights=total,
               data=dysaphis[dysaphis$clone==names(summary(dysaphis[,2]))[2] & 
                               dysaphis$date==levels(dysaphis$date)[3],],
               fct=LN.3u(),
               type="binomial")

plot(dysa_mod1,type="confidence",main=names(summary(dysaphis[,2]))[2])
plot(dysa_mod1,type="obs",add=TRUE)
plot(dysa_mod2,type="confidence",col="red",add=TRUE)
plot(dysa_mod3,type="confidence",col="blue",add=TRUE)


##############################################################################/
#additionnal code for pooling individuals from different wells ...####
##############################################################################/

#load the dataset
dysa<-read.table("data/dysa_080917.txt",header=T,sep="\t")

#reformating the data by "merging" the different wells of the same repetition
dysa2<-cbind(aggregate(dead~dose+date+clone+pesticide,data=dysa,"sum"), 
                "total"=aggregate(total~dose+date+clone+pesticide,
                                  data=dysa,"sum")$total)[,c(2,3,1,5,6,4)]


##############################################################################/
#END
##############################################################################/