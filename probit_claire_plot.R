###############################################################################
###############################################################################
#the R code for the Figure 1 of Claire's paper
###############################################################################
###############################################################################

#loading the libraries
library(drc)
library(plotrix)

#set the working directory
setwd("~/work/Rfichiers/Githuber/probit_data")


###############################################################################
#First, the plot for the Imidaclopride
###############################################################################

#load the dataset
imida<-read.table("imidata.txt",header=T,sep="\t")

#let's isolate the name of the clone and their genotype for the R81T resistance
clone_gen<-imida[imida$dose==0, 1:2]

#we change the coding for the genotype in order to have the color 
#green, orange and red associated to the RR, RT and TT genotypes
levels(imida$Rgeno)<-c(1,2,3)
colist<-c("green3","darkorange","firebrick3")

#in order to be consistent with the content of the paper, the default model 
#used for all the different clones will be 'LN.3u()'. This is the equivalent
#to the Finney method. There is no constrain on the lowerlimit and chose a 
#log-normal link is chosen. 
temp<-imida[imida$ind_ID==clone_gen[1,1] & imida$total!=0,]
temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.3u(),
              type="binomial")

plot(temp.mod,xlim=c(0,max(imida$dose)+30000),type="obs",broken=TRUE)
plot(temp.mod,xlim=c(0,100000),type="obs",broken=TRUE)
plot(temp.mod,xlim=c(0,100000),add=TRUE,lwd=5,
     type="none",col=colist[as.numeric(temp$Rgeno[1])],
     lty=as.numeric(temp$Rgeno[1])+4)

for (i in 2:dim(clone_gen)[1]) {
  temp<-imida[imida$ind_ID==clone_gen[i,1] & imida$tota!=0,]
  temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.3u(),
                type="binomial")
  plot(temp.mod,xlim=c(0,100000),type="obs",add=TRUE)
  plot(temp.mod,xlim=c(0,100000),add=TRUE,lwd=5,
       type="none",col=colist[as.numeric(temp$Rgeno[1])],
       lty=as.numeric(temp$Rgeno[1])+4)
}

#export pdf 15 x 8 inches

temp<-imida[imida$ind_ID==clone_gen[1,1] & imida$total!=0,]
temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.3u(),
              type="binomial")
#apply the Abbott correction to the death rate
Abbottemp<- 100*((temp[,5]/temp[,4])-(temp[1,5]/temp[1,4]))/
  (100-(temp[1,5]/temp[1,4]))
#when Abbott's correction lead to negative results, we turn them to 0
Abbottemp[Abbottemp<0]<-0
plot(temp.mod,xlim=c(0,max(imida$dose)+30000),type="obs",broken=TRUE)
plot(temp.mod,xlim=c(0,100000),type="obs",broken=TRUE)
plot(temp.mod,xlim=c(0,100000),add=TRUE,lwd=5,
     type="none",col=colist[as.numeric(temp$Rgeno[1])],
     lty=as.numeric(temp$Rgeno[1])+4)

for (i in 2:dim(clone_gen)[1]) {
  temp<-imida[imida$ind_ID==clone_gen[i,1] & imida$tota!=0,]
  temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.3u(),
                type="binomial")
  plot(temp.mod,xlim=c(0,100000),type="obs",add=TRUE)
  plot(temp.mod,xlim=c(0,100000),add=TRUE,lwd=5,
       type="none",col=colist[as.numeric(temp$Rgeno[1])],
       lty=as.numeric(temp$Rgeno[1])+4)
}


#in order to fix the problem of negative rate of dead aphid, we set 
#conditional model selection

op<-par(mar=c(5.1,5.5,2.1,2.1))
temp<-imida[imida$ind_ID==clone_gen[1,1] & imida$total!=0,]
if (temp[temp$dose==0,]$dead!=0) {
  temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.3u(),
                type="binomial")
} else {
  temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.2(),
                type="binomial")
}
plot(temp.mod,xlim=c(0,100000),type="obs",broken=TRUE,
     xlab=expression(paste("Log10(imidacloprid concentration) ",µg.liter^-1)),
     ylab="Percentage of mortality",cex.lab=1.5)
plot(temp.mod,xlim=c(0,100000),add=TRUE,lwd=4,
     type="none",col=colist[as.numeric(temp$Rgeno[1])])
for (i in 2:dim(clone_gen)[1]) {
  temp<-imida[imida$ind_ID==clone_gen[i,1] & imida$tota!=0,]
  if (temp[temp$dose==0,]$dead!=0) {
    temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.3u(),
                  type="binomial")
  } else {
    temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.2(),
                  type="binomial")
  }
  plot(temp.mod,xlim=c(0,100000),type="obs",add=TRUE)
  plot(temp.mod,xlim=c(0,100000),add=TRUE,lwd=4,
       type="none",col=colist[as.numeric(temp$Rgeno[1])])
}
par(op)


#same plot, but with Abbott correction applied
op<-par(mar=c(5.1,5.5,2.1,2.1))
temp<-imida[imida$ind_ID==clone_gen[1,1] & imida$total!=0,]
#apply the Abbott correction to the death rate
Abbottemp<- 100*((temp[,5]/temp[,4])-(temp[1,5]/temp[1,4]))/
  (100-(temp[1,5]/temp[1,4]))
#when Abbott's correction lead to negative results, we turn them to 0
Abbottemp[Abbottemp<0]<-0
if (temp[temp$dose==0,]$dead!=0) {
  temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.3u(),
                type="binomial")
} else {
  temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.2(),
                type="binomial")
}
plot(temp.mod,xlim=c(0,100000),type="obs",broken=TRUE,pch="",
     xlab=expression(paste("Log10(imidacloprid concentration) ",µg.liter^-1)),
     ylab="Abbott corrected \npercentage of mortality",cex.lab=1.5)
points(cbind(temp$dose,Abbottemp))
plot(temp.mod,xlim=c(0,100000),add=TRUE,lwd=4,
     type="none",col=colist[as.numeric(temp$Rgeno[1])])
for (i in 2:dim(clone_gen)[1]) {
  temp<-imida[imida$ind_ID==clone_gen[i,1] & imida$tota!=0,]
  #apply the Abbott correction to the death rate
  Abbottemp<- 100*((temp[,5]/temp[,4])-(temp[1,5]/temp[1,4]))/
    (100-(temp[1,5]/temp[1,4]))
  #when Abbott's correction lead to negative results, we turn them to 0
  Abbottemp[Abbottemp<0]<-0
  if (temp[temp$dose==0,]$dead!=0) {
    temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.3u(),
                  type="binomial")
  } else {
    temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.2(),
                  type="binomial")
  }
  plot(temp.mod,xlim=c(0,100000),type="obs",add=TRUE,pch="")
  points(cbind(temp$dose,Abbottemp))
  plot(temp.mod,xlim=c(0,100000),add=TRUE,lwd=4,
       type="none",col=colist[as.numeric(temp$Rgeno[1])])
}
par(op)

#export pdf 15 x 8 inches


###############################################################################
#Second, the plot for the thiaclopride
###############################################################################

#load the dataset
thia<-read.table("thiadata.txt",header=T,sep="\t")

#let's isolate the name of the clone and their genotype for the R81T resistance
clone_gen<-thia[thia$dose==0, 1:2]

#we change the coding for the genotype in order to have the color 
#green, orange and red associated to the RR, RT and TT genotypes
levels(thia$Rgeno)<-c(1,2,3)
colist<-c("green3","darkorange","firebrick3")

#in order to be consistent with the content of the paper, the default model 
#used for all the different clones will be 'LN.3u()'. This is the equivalent
#to the Finney method. There is no constrain on the lowerlimit and chose a 
#log-normal link is chosen. 
temp<-thia[thia$ind_ID==clone_gen[26,1] & thia$total!=0,]
temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.3u(),
              type="binomial")

plot(temp.mod,xlim=c(0,max(thia$dose)+30000),type="obs",broken=TRUE)
plot(temp.mod,xlim=c(0,70000),ylim=c(0,1),type="obs",broken=TRUE)
plot(temp.mod,xlim=c(0,70000),add=TRUE,lwd=4,
     type="none",col=colist[as.numeric(temp$Rgeno[1])])

for (i in 2:dim(clone_gen)[1]) {
  temp<-thia[thia$ind_ID==clone_gen[i,1] & thia$total!=0,]
  temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.3u(),
                type="binomial")
  plot(temp.mod,xlim=c(0,70000),type="obs",add=TRUE)
  plot(temp.mod,xlim=c(0,70000),add=TRUE,lwd=4,
       type="none",col=colist[as.numeric(temp$Rgeno[1])])
}

#export pdf 15 x 8 inches

#in order to fix the problem of negative rate of dead aphid, we set 
#conditional model selection

op<-par(mar=c(5.1,5.5,2.1,2.1))
temp<-thia[thia$ind_ID==clone_gen[26,1] & thia$total!=0,]
if (temp[temp$dose==0,]$dead!=0) {
  temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.3u(),
                type="binomial")
} else {
  temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.2(),
                type="binomial")
}

plot(temp.mod,xlim=c(0,70000),ylim=c(0,1),type="obs",broken=TRUE,
     xlab=expression(paste("Log10(thiachloprid concentration) ",µg.liter^-1)),
     ylab="Percentage of mortality",cex.lab=1.5)
plot(temp.mod,xlim=c(0,70000),add=TRUE,lwd=4,
     type="none",col=colist[as.numeric(temp$Rgeno[1])])

for (i in 2:dim(clone_gen)[1]) {
  temp<-thia[thia$ind_ID==clone_gen[i,1] & thia$total!=0,]
  if (temp[temp$dose==0,]$dead!=0) {
    temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.3u(),
                  type="binomial")
  } else {
    temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.2(),
                  type="binomial")
  }
  plot(temp.mod,xlim=c(0,70000),type="obs",add=TRUE)
  plot(temp.mod,xlim=c(0,70000),add=TRUE,lwd=4,
       type="none",col=colist[as.numeric(temp$Rgeno[1])])
}

par(op)

#export pdf 15 x 8 inches


#same graph but with Abbott correction

op<-par(mar=c(5.1,5.5,2.1,2.1))
temp<-thia[thia$ind_ID==clone_gen[26,1] & thia$total!=0,]
#apply the Abbott correction to the death rate
Abbottemp<- 100*((temp[,5]/temp[,4])-(temp[1,5]/temp[1,4]))/
  (100-(temp[1,5]/temp[1,4]))
#when Abbott's correction lead to negative results, we turn them to 0
Abbottemp[Abbottemp<0]<-0
if (temp[temp$dose==0,]$dead!=0) {
  temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.3u(),
                type="binomial")
  } else {
    temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.2(),
                  type="binomial")
}

plot(temp.mod,xlim=c(0,70000),ylim=c(0,1),type="obs",broken=TRUE,
xlab=expression(paste("Log10(thiacloprid concentration) ",µg.liter^-1)),
ylab="Abbott corrected \npercentage of mortality",cex.lab=1.5,pch="")
points(cbind(temp$dose,Abbottemp))
plot(temp.mod,xlim=c(0,70000),add=TRUE,lwd=4,
type="none",col=colist[as.numeric(temp$Rgeno[1])])

for (i in 2:dim(clone_gen)[1]) {
  temp<-thia[thia$ind_ID==clone_gen[i,1] & thia$total!=0,]
  #apply the Abbott correction to the death rate
  Abbottemp<- 100*((temp[,5]/temp[,4])-(temp[1,5]/temp[1,4]))/
    (100-(temp[1,5]/temp[1,4]))
  #when Abbott's correction lead to negative results, we turn them to 0
  Abbottemp[Abbottemp<0]<-0
  if (temp[temp$dose==0,]$dead!=0) {
    temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.3u(),
                  type="binomial")
    } else {
      temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.2(),
                    type="binomial")
      }
  plot(temp.mod,xlim=c(0,70000),type="obs",add=TRUE,pch="")
  points(cbind(temp$dose,Abbottemp))
  plot(temp.mod,xlim=c(0,70000),add=TRUE,lwd=4,
       type="none",col=colist[as.numeric(temp$Rgeno[1])])
}

par(op)

#export pdf 15 x 8 inches


###############################################################################
#plot with the probit on the y-axis
###############################################################################

#imida plot, with Abbott correction applied and probit on the y-axis
op<-par(mar=c(5.1,5.5,2.1,2.1))
temp<-imida[imida$ind_ID==clone_gen[1,1] & imida$total!=0,]
#apply the Abbott correction to the death rate
Abbottemp<- 100*((temp[,5]/temp[,4])-(temp[1,5]/temp[1,4]))/
  (100-(temp[1,5]/temp[1,4]))
#when Abbott correction lead to negative results, we turn them to 0
Abbottemp[Abbottemp<0]<-0
qAbbott<-qnorm(Abbottemp)
temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.2(),
                type="binomial")
plot(temp.mod,xlim=c(0,100000),ylim=c(-3,3),type="obs",broken=TRUE,pch="",
     xlab=expression(paste("Log10(imidacloprid concentration) ",µg.liter^-1)),
     ylab="probit of Abbott corrected \npercentage of mortality",cex.lab=1.5)
points(cbind(temp$dose,qAbbott))
plotemp<-plot(temp.mod,new=FALSE,xlim=c(0,100000),add=TRUE,lwd=0,
     type="none",col=colist[as.numeric(temp$Rgeno[1])],lty=0)
lines(plotemp[,1],qnorm(plotemp[,2]),col=colist[as.numeric(temp$Rgeno[1])],
      lwd=4)

for (i in 2:dim(clone_gen)[1]) {
  temp<-imida[imida$ind_ID==clone_gen[i,1] & imida$tota!=0,]
  #apply the Abbott correction to the death rate
  Abbottemp<- 100*((temp[,5]/temp[,4])-(temp[1,5]/temp[1,4]))/
    (100-(temp[1,5]/temp[1,4]))
  #when Abbott correction lead to negative results, we turn them to 0
  Abbottemp[Abbottemp<0]<-0
  qAbbott<-qnorm(Abbottemp)
  temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.2(),
                  type="binomial")
  plot(temp.mod,xlim=c(0,100000),type="obs",add=TRUE,pch="")
  points(cbind(temp$dose,qAbbott))
  plotemp<-plot(temp.mod,new=FALSE,xlim=c(0,100000),add=TRUE,lwd=0,
                type="none",col=colist[as.numeric(temp$Rgeno[1])],lty=0)
  plotemp[plotemp[,2]<0,2]<-0.0001
  lines(plotemp[,1],qnorm(plotemp[,2]),col=colist[as.numeric(temp$Rgeno[1])],
        lwd=4)
  }
par(op)

#export pdf 15 x 8 inches


#imida plot, with Abbott correction applied and probit on the y-axis
op<-par(mar=c(5.1,5.5,2.1,2.1))
temp<-thia[thia$ind_ID==clone_gen[26,1] & thia$total!=0,]
#apply the Abbott correction to the death rate
Abbottemp<- 100*((temp[,5]/temp[,4])-(temp[1,5]/temp[1,4]))/
  (100-(temp[1,5]/temp[1,4]))
#when Abbott's correction lead to negative results, we turn them to 0
Abbottemp[Abbottemp<0]<-0
qAbbott<-qnorm(Abbottemp)
temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.2(),
                type="binomial")
plot(temp.mod,xlim=c(0,70000),ylim=c(-3,3),type="obs",broken=TRUE,
     xlab=expression(paste("Log10(thiachloprid concentration) ",µg.liter^-1)),
     ylab="probit of Abbott corrected \npercentage of mortality",cex.lab=1.5,pch="")
points(cbind(temp$dose,qAbbott))
plotemp<-plot(temp.mod,new=FALSE,xlim=c(0,100000),add=TRUE,lwd=0,
              type="none",col=colist[as.numeric(temp$Rgeno[1])],lty=0)
lines(plotemp[,1],qnorm(plotemp[,2]),col=colist[as.numeric(temp$Rgeno[1])],
      lwd=4)

for (i in 2:dim(clone_gen)[1]) {
  temp<-thia[thia$ind_ID==clone_gen[i,1] & thia$total!=0,]
  #apply the Abbott correction to the death rate
  Abbottemp<- 100*((temp[,5]/temp[,4])-(temp[1,5]/temp[1,4]))/
    (100-(temp[1,5]/temp[1,4]))
  #when Abbott's correction lead to negative results, we turn them to 0
  Abbottemp[Abbottemp<0]<-0
  qAbbott<-qnorm(Abbottemp)
  temp.mod<-drm(dead/total~dose,weights=total,data=temp,fct=LN.2(),
                  type="binomial")
  plot(temp.mod,xlim=c(0,70000),type="obs",add=TRUE,pch="")
  points(cbind(temp$dose,qAbbott))
  plotemp<-plot(temp.mod,xlim=c(0,70000),add=TRUE,lwd=0,
                type="none",col=colist[as.numeric(temp$Rgeno[1])],lty=0)
  plotemp[plotemp[,2]<0,2]<-0.0001
  lines(plotemp[,1],qnorm(plotemp[,2]),col=colist[as.numeric(temp$Rgeno[1])],
        lwd=4)
  }
par(op)


###############################################################################
#END
###############################################################################