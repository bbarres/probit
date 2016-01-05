###############################################################################
###############################################################################
#the R code for the Figure of Claire's paper
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