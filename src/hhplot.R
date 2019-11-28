dr<-getwd()
loc<-substr(dr,1,nchar(dr))#NB this is the directory the script is called from! not where it is!!
########
#FIRST, BY GENDER
## nmcname<-paste(loc,"/data/nmcounts.dat",sep="")
## NMC<-read.delim(file=nmcname,sep=" ",header=FALSE)
## nmc<-as.matrix(NMC)
## nmc<-nmc/sum(nmc)
## nmcvec<-0
## nvec<-rep(0,10)
## #less than 10 phh
## for( i in 1:10 ){
##   for(j in 1:10){
##     if((i+j)<12){
##       nmcvec<-c(nmcvec,nmc[i,j])	
##       nvec[i+j-1]<-nvec[i+j-1]+nmc[i,j]
##     }
##   }	
## }
## nmcvec<-nmcvec[-1]

## NMC now not used

fname<-paste(loc,"/data/dss1.dat",sep="")
first<-read.csv(file=fname, sep=",", comment.char=";",header=FALSE)
F<-as.matrix(first)


lname<-paste(loc,"/data/dss2.dat",sep="")
last<-read.csv(file=lname, sep=",", comment.char=";",header=FALSE)
L<-as.matrix(last)

H<-matrix(0,nrow=2,ncol=55)
## H[1,]<-nmcvec
#less than 10phh
H[1,]<-t(F[((F[,1]+F[,2])<10),3])
H[2,]<-t(L[((L[,1]+L[,2])<10),3])

Iz<-F[((F[,1]+F[,2])<10),1]
Jz<-F[((F[,1]+F[,2])<10),2]

nmz<-""
for( i in 1:55 ){
  nmz<-c(nmz,paste("(",Iz[i],",",Jz[i],")",sep=""))	
}
nmz<-nmz[-1]

##########
#NOW NOT BY GENDER
H0<-matrix(0,nrow=2,ncol=10)
## H0[1,]<-nvec
for(i in 1:length(Iz)){
  H0[1,(Iz[i]+Jz[i]+1)]<-H0[1,(Iz[i]+Jz[i]+1)] + H[1,i]
  H0[2,(Iz[i]+Jz[i]+1)]<-H0[2,(Iz[i]+Jz[i]+1)] + H[2,i]
}
H1<-H0[,-1]


###########
#ACTUALLY PLOTTING
plot1<-paste(loc,"/plots/hhgendersize.pdf",sep="")
pdf(file=plot1,width=11.7,height=8.3)
par(las=2)
barplot(H,beside=TRUE,main="Household size distribution\n (by gender)",names.arg=nmz,legend.text=c("start","end"),xlab="number of (men,women)")
dev.off()

#H1[1,]<-H1[1,]/sum(H1[1,])
#H1[2,]<-H1[2,]/sum(H1[2,])
#H1[3,]<-H1[3,]/sum(H1[3,])
plot2<-paste(loc,"/plots/hhsize.pdf",sep="")
pdf(file=plot2)
par(las=1)
barplot(H1,beside=TRUE,main="Household size distribution",names.arg=1:9,legend.text=c("start","end"),xlab="household size")
dev.off()
#print(sum(H1[1,]))
#print(sum(H1[2,]))
#print(sum(H1[3,]))
