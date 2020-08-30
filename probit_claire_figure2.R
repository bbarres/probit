##############################################################################/
##############################################################################/
#the R code for the Figure 2 of Claire's paper
##############################################################################/
##############################################################################/


##############################################################################/
#plotting the LC50 with imdaclopride against LC50 with thiaclopride####
##############################################################################/

op<-par(mar=c(5.1,5.1,2.1,2.1))
ITdata<-read.table("data/imida-thia.txt",header=TRUE,sep="\t")
plot(ITdata,log="xy",pch=19,bg="black",cex=1.5,cex.lab=1.5,las=1,
     xlab=expression(paste("imidacloprid LC50 ",µg.liter^-1)),
     ylab=expression(paste("thiacloprid LC50 ",µg.liter^-1)))
rezlm<-lm(log10(ITdata$thiacloprid.LC50)~log10(ITdata$imidacloprid.LC50))
abline(rezlm,col="grey60",lwd=2)
par(op)

#export to pdf file 8 x 6 inches


##############################################################################/
#END
##############################################################################/
