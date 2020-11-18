##############################################################################/
##############################################################################/
#code for the drosophila spinosad tests
##############################################################################/
##############################################################################/

#loading the libraries
library(drc)
library(plotrix)


#load the dataset
drosoSpino<-read.table("data/spin201117.txt",header=T,sep="\t")

#combining moribund and dead into dead
drosoSpino$dead<-drosoSpino$dead+drosoSpino$moribund

#reformating the data by "merging" the different wells of the same repetition
dysaphis<-cbind(aggregate(dead~dose+date+clone+pesticide,data=dysaphis,"sum"), 
                "total"=aggregate(total~dose+date+clone+pesticide,
                                  data=dysaphis,"sum")$total)[,c(2,3,1,5,6,4)]


##############################################################################/
#Drosophila suzukii vs spinosad first test####
##############################################################################/

droSpi_fem<-drm(dead/total~dose,weights=total,
               data=drosoSpino[drosoSpino$sex=="female",],
               fct=LN.3u(),
               type="binomial")
droSpi_mal<-drm(dead/total~dose,weights=total,
                data=drosoSpino[drosoSpino$sex=="male",],
                fct=LN.3u(),
                type="binomial")

op<-par(mfrow=c(2,1))
plot(droSpi_fem,type="confidence",main="spinosad - female")
plot(droSpi_fem,type="obs",add=TRUE)
plot(droSpi_fem,type="average",col="darkblue",pch=19,add=TRUE)
segments(ED(droSpi_fem,50,interval="delta")[1],-0.1,
         ED(droSpi_fem,50,interval="delta")[1],0.5,
         col="red",lwd=2)
segments(0.001,0.5,
         ED(droSpi_fem,50,interval="delta")[1],0.5,
         col="red",lwd=2)
text(50,0.4,
     labels=round(ED(droSpi_fem,50,interval="delta")[1],2),cex=2)

plot(droSpi_mal,type="confidence",main="spinosad - male")
plot(droSpi_mal,type="obs",add=TRUE)
plot(droSpi_mal,type="average",col="darkblue",pch=19,add=TRUE)
segments(ED(droSpi_mal,50,interval="delta")[1],-0.1,
         ED(droSpi_mal,50,interval="delta")[1],0.5,
         col="red",lwd=2)
segments(0.001,0.5,
         ED(droSpi_mal,50,interval="delta")[1],0.5,
         col="red",lwd=2)
text(50,0.4,
     labels=round(ED(droSpi_mal,50,interval="delta")[1],2),cex=2)
par(op)

#export to .pdf 10 x 5 inches


##############################################################################/
#END
##############################################################################/