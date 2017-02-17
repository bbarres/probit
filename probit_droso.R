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
#example for Drosophila suzukii resistance to Difenoconazol
###############################################################################

#load the dataset
droso<-read.table("Droso2.txt",header=T,sep="\t")

#in order to obtain the same results than with priprobit, the Finney equivalent
#method, we have to remove the constrain on lowerlimit and chose a log-normal 
#model instead of a log-logistic model
droso.m1<-drm(mort/total~dose,weights=total,data=droso,fct=LN.3u(),
              type="binomial")
plot(droso.m1)
plot(droso.m1,type="confidence")

#the ED50 obtained is identical to the one obtained with priprobit, the SD 
#still differ a little bit
ed50val_droso1<-ED(droso.m1,50,interval="delta",reference="control")

droso.m2<-drm(mort/total~dose+age,weights=total,data=droso,fct=LN.3u(),
              type="binomial")
plot(droso.m2)
plot(droso.m2,type="confidence")

#to test for the effect of age and sexe, we use a logistic regression model

library(dplyr)
drosocomb<-droso %>%
           group_by(age,sexe,dose) %>% 
           summarise(mort = sum(mort), total=sum(total))


test<-glm(cbind((drosocomb$total-drosocomb$mort),drosocomb$mort)~dose*age*sexe,
          family=binomial(link=probit),data=drosocomb)
summary(test)

test<-glm(cbind((drosocomb$total-drosocomb$mort),drosocomb$mort)~dose+age*sexe,
          family=binomial(link=probit),data=drosocomb)
summary(test)

#without combining the test by date
test<-glm(cbind((droso$total-droso$mort),droso$mort)~dose+age*sexe*date,
          family=binomial(link=probit),data=droso)
summary(test)

