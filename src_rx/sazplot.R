dr<-getwd()
loc<-substr(dr,1,nchar(dr))#NB this is the directory the script is called from! not where it is!!

fname<-paste(loc,"/data/dss1.dat",sep="")



## ----------------- plotting the time-series
##setwd(basedr)

P <- read.csv(file='shared/TBtargetsSA.csv')
load('shared/Stargs.Rdata')
## HIV prevalence graph
zth <- read.csv(file='shared/saps.csv')
zsat <- read.csv(file='shared/ZAMSTARtargetsSA.csv',header=F)
rownames(zsat) <- zsat[,1]
fs <- Stargs$zshiv / zth[21,'hivp']


R5 <- read.table(file='data/results5.dat',sep='\t',header=FALSE)
R4 <- read.table(file='data/results4.dat',sep='\t',header=FALSE)
R3 <- read.table(file='data/results3.dat',sep='\t',header=FALSE)
R2 <- read.table(file='data/results2.dat',sep='\t',header=FALSE)
R1 <- read.table(file='data/results1.dat',sep='\t',header=FALSE)
tz <- R4[,1]


setwd('plots')
## data
n <- dim(P)[1]
nn <- length(tz)
ff <- R4[nn,5]/P$inc100k[n]

pdf(file='ZS1.pdf')
par(xaxs='i',yaxs='i',bty='l',las=1,mar=c(5.1,4.1,4.1,4.1))
plot(x=c(1980.0,2010.5),y=c(0,1.15*max(R4[,5],ff*P$inchi)),type='n',xlab='Year',
     ylab='Incidence rate (per 100k per year)',
     main="", axes=F )
lines(tz,R4[,5],col='blue')
points(P$year,P$inc100k * ff,col='blue')
lines(P$year,P$inchi * ff, lty=2,col='blue')
lines(P$year,P$inclo * ff, lty=2,col='blue')
## HIV
lines(tz,R2[,7]*R4[,2],col='red')
points(P$year,P$tbhivinc100k * ff,col='red')
lines(P$year,P$tbhinchi * ff, lty=2,col='red')
lines(P$year,P$tbhinclo * ff, lty=2,col='red')
axis(1);axis(2)
axis(4,labels=c(0,P$inc100k[n]),at=c(0,R4[nn,5]),lty=2);mtext('WHO national incidence (per 100k per year)',side=4,las=0,line=2.5)
legend('topleft',col=c('blue','red','black','black'),
       lty=c(-1,-1,1,2),pch=c(22,22,-1,1),pt.bg=c('blue','red',NA,NA),pt.cex=c(2,2,1.5,1),
       bty='n',legend=c('Total','HIV-TB','model','WHO'),lwd=1.5)
dev.off()


## data prevalence
pdf(file='ZS2.pdf')
par(xaxs='i',yaxs='i',bty='l',las=1)
par(xaxs='i',yaxs='i',bty='l',las=1)
plot(x=c(1980.0,2010.5),y=c(0,4005),type='n',xlab='Year',
     ylab='TB prevalence (per 100k)',
     main="", axes=F )
n <- dim(P)[1]
lines(tz,R4[,2],col='blue')                        #prevalence
lines(tz,R4[,4],col='green')                        #treatment
## TS2&6
axis(1);axis(2)
points(P$year,P$prev100k * Stargs$zstbp/P$prev100k[n])
points(2010,Stargs$zstbp,lwd=2,col='blue',cex=1.5);
lines(c(2010,2010),(1-0.5*c(-Stargs$zstbpw,+Stargs$zstbpw))*Stargs$zstbp,lwd=2,col='blue') #todo
eps <- .1
points(2010-eps,Stargs$zsrx,lwd=2,col='green',cex=1.5);
lines(c(2010,2010)-eps,(1-0.5*c(-Stargs$zsrxw,+Stargs$zsrxw))*Stargs$zsrx,lwd=2,col='green') #todo
lines(P$year,P$prevhi * Stargs$zstbp/P$prev100k[n], lty=2)
lines(P$year,P$prevlo * Stargs$zstbp/P$prev100k[n], lty=2)
legend('topleft',fill=c('blue','green'),bty='n',legend=c('Undiagnosed TB','On TB treatment'))
dev.off()



## HIV prevalence for ZAMSTAR SA
pdf('ZS3.pdf')
ytop <- 30
plot(c(1980,2010),c(0,ytop),xlab='year',ylab = 'HIV prevalence 18+ (%)',type='n',main='SA ZAMSTAR')
points(zth$year,fs*zth$hivp)
for(i in 1:dim(zth)[1]){lines(c(zth$year[i],zth$year[i]),c(fs*zth$hivplo[i],fs*zth$hivphi[i]))}
lines(R1[,1],100*R1[,8])
points(2010,zsat['hiv',2],lwd=2,col='red'); lines(c(2010,2010),c(zsat['hivlo',2],zsat['hivhi',2]),lwd=2,col='red')
lines(R5[,1],100*R5[,7],col='blue')
eps <- .2
points(2010+eps,zsat['art',2],lwd=2,col='blue'); lines(c(2010,2010)+eps,c(zsat['artlo',2],zsat['arthi',2]),lwd=2,col='blue')
legend('topleft',col=c('black','black','blue','blue','red'),pch=c(1,-1,1,-1,1),lty=1,lwd=c(1,1,2,1,2),legend=c('HIV prevalence (adj. UNAIDS)','HIV prevalence (model)','ART coverage (Z\'STAR)','ART coverage (model)','HIV prevalence (Z\'STAR)'),bty='n')
dev.off()


pdf('ZS4.pdf')
par(xaxs='i',yaxs='i',bty='l',las=1)
plot(x=c(1980.0,2010.5),y=c(0,25),type='n',xlab='Year',
     ylab='ARI (% per year)',
     main="", axes=F )
n <- dim(P)[1]
lines(tz,100*R3[,3],col='blue')             #ARI
## 3,3
axis(1);axis(2)
points(2010,Stargs$zsari,col='blue',cex=1.5,lwd=2)
lines(c(2010,2010),(1+c(Stargs$zsariw,-Stargs$zsariw))*Stargs$zsari, lwd=2,col='blue')
dev.off()

 
