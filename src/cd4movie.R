## this is to make a movie of the cd4 distribution through time

setwd("data")
CD4 <- read.table(file="results4.dat",header=FALSE,sep="\t")
CD4 <- CD4[,-ncol(CD4)]

n <- 20
nmz <- paste("[",seq(from=0,length=n,by=50),",",seq(from=50,length=n,by=50),")",sep="")

#names(tmp) <- nmz
system("if [ ! -d 'tmp' ]; then mkdir tmp; fi") #make dir if doesn't exist
setwd("tmp")

#make graphs
for(i in 1:nrow(CD4)){
  png(sprintf("%03d.png",i))
  tmp <- as.numeric(CD4[i,6:ncol(CD4)])
  par(las=3,mar=c(7, 4, 4, 2)+.1)
  barplot(tmp,names.arg=nmz,main=paste("t=",format(CD4[i,1],digits=5,nsmall=1),sep=""))
  par(las=0,mar=c(5, 4, 4, 2)+.1)
  mtext("CD4 count per 10^6L",side=1,line=3)
  dev.off()
}

## make movie
mvcmd <- "ffmpeg -qscale 5 -r 12 -b 9600 -i %03d.png cd4movie.mp4"
system(mvcmd,ignore.stdout=T,ignore.stderr=T)
system("mv cd4movie.mp4 ../../plots/")

## ## tidy
system("rm *")
