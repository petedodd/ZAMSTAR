## this is for analysing the population snapshot generated at the end of the
## IBM run

setwd("data")


## and some time-series outputs
D3 <- read.table(file='results3.dat',header=F,sep="\t")
hcop <- D3[,11]
tbcop <- D3[,10]*1e5
D4 <- read.table(file='results4.dat',header=F,sep="\t")
D5 <- read.table(file='results5.dat',header=F,sep="\t")
tbp <- D4[,2]
D1 <- read.table(file='results1.dat',header=F,sep="\t")
D2 <- read.table(file='results2.dat',header=F,sep="\t")
hp <- D1[,4]                            #HIV 15-49, 8 is 18+

tz <- D1[,1]
tbr <- tbcop/tbp
hivr <- hcop/hp
hivr <- hivr[hp>1e-5]

## SA  HIV prevalence graph
zsat <- read.csv(file='../shared/ZAMSTARtargetsSA.csv',header=F)
rownames(zsat) <- zsat[,1]
H <- read.csv(file='../shared/AIDSinfoHIV.csv') #SA national
zth <- H[H$country=='ZA',c('year','hivp','hivplo','hivphi')]

print('plotting SA national HIV timeseries...')
pdf('../plots/SAnationalHIV.pdf')
ytop <- 30
plot(c(1980,2010.5),c(0,ytop),xlab='year',ylab = 'HIV prevalence 15-49 (%)',type='n',main='SA')
points(zth$year,zth$hivp)
for(i in 1:dim(zth)[1]){lines(c(zth$year[i],zth$year[i]),c(zth$hivplo[i],zth$hivphi[i]))}
lines(D1[,1],100*D1[,4])
#points(2010,zsat['hiv',2],lwd=2,col='red'); lines(c(2010,2010),c(zsat['hivlo',2],zsat['hivhi',2]),lwd=2,col='red')
lines(D5[,1],100*D5[,7],col='blue')
eps <- .2
points(2010+eps,zsat['art',2],lwd=2,col='blue'); lines(c(2010,2010)+eps,c(zsat['artlo',2],zsat['arthi',2]),lwd=2,col='blue')
legend('topleft',col=c('black','black','blue','blue'),pch=c(1,-1,1,-1),lty=1,lwd=c(1,1,2,1),legend=c('HIV prevalence','HIV prevalence (model)','ART coverage (Z\'STAR)','ART coverage (model)'),bty='n')
dev.off()

## SA incidence plot
SI <- read.table(file='../shared/GDsainc.dat',sep = '\t');
colnames(SI) <- c('year','inc','inclo','inchi')
SIH <- read.table(file='../shared/GDsahinc.dat',sep = '\t');
colnames(SIH) <- c('year','inc','inclo','inchi')
pdf('../plots/SAnationalTBincidence.pdf')
par(xaxs='i',yaxs='i',bty='l',las=1)
plot(x=c(1980.0,2010.5),y=c(0,1.15*max(D4[,5],SI$inchi)),type='n',xlab='Year',
     ylab='Incidence rate (per 100k per year)',
     main="", axes=F )
lines(tz,D4[,5],col='blue')
points(SI$year,SI$inc,col='blue')
lines(SI$year,SI$inchi, lty=2,col='blue')
lines(SI$year,SI$inclo, lty=2,col='blue')
## HIV
lines(tz,D2[,7]*D4[,2],col='red')
points(SIH$year,SIH$inc,col='red')
lines(SIH$year,SIH$inchi, lty=2,col='red')
lines(SIH$year,SIH$inclo, lty=2,col='red')
axis(1);axis(2)
legend('topleft',fill=c('blue','red'),
       bty='n',legend=c('Total','HIV-TB'))
dev.off()

SP <- read.table(file='../shared/GDsainc.dat',sep = '\t');
colnames(SP) <- c('year','prev','prevlo','prevhi')

## SA prevalence plot
pdf('../plots/SAnationalTBprevalence.pdf')
par(xaxs='i',yaxs='i',bty='l',las=1)
plot(x=c(1980.0,2010.5),y=c(0,1.15*max(D4[,2],SP$prevhi)),type='n',xlab='Year',
     ylab='Prevalece (per 100k)',
     main="", axes=F )
lines(tz,D4[,2],col='blue')
points(SP$year,SP$prev,col='blue')
lines(SP$year,SP$prevhi, lty=2,col='blue')
lines(SP$year,SP$prevlo, lty=2,col='blue')
axis(1);axis(2)
dev.off()



library(ggplot2)
theme_set(theme_bw())
df <- data.frame(time= D1[,1], RR = tbr )
tit <- paste("TB household clustering: RR=",format(mean(tbr,na.rm=T),nsmall=2,digits=2))
qplot(data=df, x=time,y=RR)+geom_smooth() + opts(title=tit)
print(tit)
ggsave("../plots/hhTB.pdf")
df <- data.frame(time= D1[hp>1e-5,1], RR = hivr )
tit <- paste("HIV household clustering: RR=",format(mean(hivr,na.rm=T),nsmall=2,digits=2))
qplot(data=df, x=time,y=RR)+geom_smooth() + ylim(0,8) + opts(title=tit)
print(tit)
ggsave("../plots/hhHIV.pdf")

## and now the popualtion snapshot stuff...

PSSraw <- read.csv(file="PopSnapshot.csv",header=TRUE)
PSS <- PSSraw[PSSraw$isX==0,]           #alive!
## hist(as.numeric(PSS[PSS[,"age"]<80,"age"]))
## hist(as.numeric(PSS[PSS[,"TBu"]==1,"age"]))
## need to extract for given age groups and work out ratio...

atop <- 70
restr <- (as.numeric(PSS[,"age"])<atop)
PSSr <- PSS[restr,]
agebrz <- seq(from=0,to=atop,by=5)
TBu <- (PSSr[,"TBu"]==1)
HIVp <- (PSSr[,"HIV"]==1)
SMRp <- (PSSr[,"smrp"]==1)
male <- (PSSr[,"gender"]==1)
ltbi <- (PSSr[,"LTBI"]==1)

denomh <- hist(as.numeric(PSSr[HIVp,"age"]),breaks=agebrz,plot=FALSE)
## denoms <- hist(as.numeric(PSSr[,"age"]),breaks=agebrz,plot=FALSE)
denom <- hist(as.numeric(PSSr[,"age"]),breaks=agebrz,plot=FALSE)
denomm <- hist(as.numeric(PSSr[male,"age"]),breaks=agebrz,plot=FALSE)
denomf <- hist(as.numeric(PSSr[!male,"age"]),breaks=agebrz,plot=FALSE)
denommh <- hist(as.numeric(PSSr[male & HIVp,"age"]),breaks=agebrz,plot=FALSE)
denomfh <- hist(as.numeric(PSSr[!male & HIVp,"age"]),breaks=agebrz,plot=FALSE)

TB <- hist(as.numeric(PSSr[TBu,"age"]),breaks=agebrz,plot=FALSE) 
TBh <- hist(as.numeric(PSSr[TBu & HIVp,"age"]),breaks=agebrz,plot=FALSE)
TBs <- hist(as.numeric(PSSr[TBu & SMRp,"age"]),breaks=agebrz,plot=FALSE) 
TBm <- hist(as.numeric(PSSr[TBu & male,"age"]),breaks=agebrz,plot=FALSE)
TBf <- hist(as.numeric(PSSr[TBu & !male,"age"]),breaks=agebrz,plot=FALSE) 
TBmh <- hist(as.numeric(PSSr[TBu & male & HIVp,"age"]),breaks=agebrz,plot=FALSE)
TBfh <- hist(as.numeric(PSSr[TBu & !male & HIVp,"age"]),breaks=agebrz,plot=FALSE)

ltbi <- hist(as.numeric(PSSr[ltbi,"age"]),breaks=agebrz,plot=FALSE)

## ## NB kids are large numbers but small proportions!
## ## basic plot of proportions
## ## plot(agebrz[-1],TBc$counts/denom$counts)
## bpdata <- barplot(1000*(TB$counts/denom$counts),
##                   xlab="age group",ylab="prevalence undiagnosed TB/100,000",
##                   main="")
## bp <- (bpdata + 0.5*diff(bpdata)[1])
## bp <- c(bp[1]-diff(bp)[1],bp)
## axis(side=1,at=bp,labels=agebrz)


## ;;;;;;;;;;;; PLOTTING
setwd("../plots")
pdf("snapshot_TBprevalence.pdf",width=10,height=10)
par(mfrow=c(2,2))
## HIV+/-
TBH <- rbind((TB$counts-TBh$counts)/denom$counts, TBh$counts/denom$counts)
bpdata <- barplot(100000*(TBH),xlab="age group",
                  ylab="prevalence undiagnosed TB/100,000",beside=FALSE,
                  legend=c("HIV-","HIV+"),main="Total, by HIV status")
bp <- (bpdata + 0.5*diff(bpdata)[1])
bp <- c(bp[1]-diff(bp)[1],bp)
axis(side=1,at=bp,labels=agebrz)


## smr+/-
TBS <- rbind((TB$counts-TBs$counts)/denom$counts, TBs$counts/denom$counts)
bpdata <- barplot(100000*(TBS),xlab="age group",
                  ylab="prevalence undiagnosed TB/100,000",beside=FALSE,
                  legend=c("Smr-","Smr+"),main="Total, by smear status")
bp <- (bpdata + 0.5*diff(bpdata)[1])
bp <- c(bp[1]-diff(bp)[1],bp)
axis(side=1,at=bp,labels=agebrz)

## ## men, women
## ## men
## bpdata <- barplot(1000*(TBm$counts/denomm$counts),xlab="age group",ylab="prevalence undiagnosed TB/100,000")
## bp <- (bpdata + 0.5*diff(bpdata)[1])
## bp <- c(bp[1]-diff(bp)[1],bp)
## axis(side=1,at=bp,labels=agebrz)
## ## women
## bpdata <- barplot(1000*(TBf$counts/denomf$counts),xlab="age group",ylab="prevalence undiagnosed TB/100,000")
## bp <- (bpdata + 0.5*diff(bpdata)[1])
## bp <- c(bp[1]-diff(bp)[1],bp)
## axis(side=1,at=bp,labels=agebrz)

## men, women, HIV+/-
## women
TBFH <- rbind((TBf$counts-TBfh$counts)/denomf$counts, TBfh$counts/denomf$counts)
bpdata <- barplot(100000*(TBFH),xlab="age group",
                  ylab="prevalence undiagnosed TB/100,000",beside=FALSE,
                  legend=c("HIV-","HIV+"),main="Women, by HIV status")
bp <- (bpdata + 0.5*diff(bpdata)[1])
bp <- c(bp[1]-diff(bp)[1],bp)
axis(side=1,at=bp,labels=agebrz)


## men
TBMH <- rbind((TBm$counts-TBmh$counts)/denomm$counts, TBmh$counts/denomm$counts)
bpdata <- barplot(100000*(TBMH),xlab="age group",
                  ylab="prevalence undiagnosed TB/100,000",beside=FALSE,
                  legend=c("HIV-","HIV+"),main="Men, by HIV status")
bp <- (bpdata + 0.5*diff(bpdata)[1])
bp <- c(bp[1]-diff(bp)[1],bp)
axis(side=1,at=bp,labels=agebrz)

dev.off()
par(mfrow=c(1,1))

## ;;;;;;;;;;;;;;HIV plots
pdf("snapshot_HIVprevalence.pdf")
par(mfrow=c(2,1))
## women
HIVFH <- denomfh$counts/denomf$counts
bpdata <- barplot(100*(HIVFH),xlab="age group",
                  ylab="prevalence HIV (%)",
                  main="HIV prevalence, women")
bp <- (bpdata + 0.5*diff(bpdata)[1])
bp <- c(bp[1]-diff(bp)[1],bp)
axis(side=1,at=bp,labels=agebrz)
## men
HIVMH <- denommh$counts/denomm$counts
bpdata <- barplot(100*(HIVMH),xlab="age group",
                  ylab="prevalence HIV (%)",
                  main="HIV prevalence, men")
bp <- (bpdata + 0.5*diff(bpdata)[1])
bp <- c(bp[1]-diff(bp)[1],bp)
axis(side=1,at=bp,labels=agebrz)

par(mfrow=c(1,1))
dev.off()


## ;;;;;;;;;;;;;;LTBI plots
pdf("snapshot_LTBIprevalence.pdf")
#par(mfrow=c(2,1))
## women
LTBI <- ltbi$counts/denom$counts
bpdata <- barplot(100*(LTBI),xlab="age group",
                  ylab="LTBI (%)",
                  main="LTBI prevalence")
bp <- (bpdata + 0.5*diff(bpdata)[1])
bp <- c(bp[1]-diff(bp)[1],bp)
axis(side=1,at=bp,labels=agebrz)
## men
## HIVMH <- denommh$counts/denomm$counts
## bpdata <- barplot(100*(HIVMH),xlab="age group",
##                   ylab="prevalence HIV (%)",
##                   main="HIV prevalence, men")
## bp <- (bpdata + 0.5*diff(bpdata)[1])
## bp <- c(bp[1]-diff(bp)[1],bp)
## axis(side=1,at=bp,labels=agebrz)

#par(mfrow=c(1,1))
dev.off()

## ;;;;;;;;;;;;;;;HOUSEHOLD PLOTS
hhdat <- table(PSS[,"TBu"],PSS[,"hhsize"])
hhprop <- hhdat[2,] / colSums(hhdat)
pdf("snapshot_TBprevalenceByHHsize.pdf")
plot(100*hhprop,xlab="household size",ylab="TB prevalence (%)",main="TB prevalence by household size")
dev.off()

## ;;;;;;;;;;;PYRAMID PLOTS
library("plotrix")

## absolute burden
## pyramid.plot(100*TBm$counts/sum(TB$counts),100*TBf$counts/sum(TB$counts), labels=lbz,gap=1.5,lxcol="darkgray",rxcol="darkgray",main="Absolute burden of prevalent TB")
pdf("snapshot_absTB.pdf")
bpdata <- barplot(100*(TB$counts/sum(TB$counts)),xlab="age group",
                  ylab="Proportion of prevalent TB (%)",
                  main="Absolute burden pulmonary of TB")
bp <- (bpdata + 0.5*diff(bpdata)[1])
bp <- c(bp[1]-diff(bp)[1],bp)
axis(side=1,at=bp,labels=agebrz)
dev.off()

pdf("snapshot_absLTBI.pdf")
bpdata <- barplot(100*(ltbi$counts/sum(ltbi$counts)),xlab="age group",
                  ylab="Proportion of LTBI (%)",
                  main="Absolute burden of LTBI")
bp <- (bpdata + 0.5*diff(bpdata)[1])
bp <- c(bp[1]-diff(bp)[1],bp)
axis(side=1,at=bp,labels=agebrz)
dev.off()


## demography
pdf("snapshot_demography.pdf")
stuff <- c(100*denomm$counts/sum(denomm$counts),100*denomf$counts/sum(denomf$counts))
lbz <- paste("(",agebrz[-length(agebrz)],"-", agebrz[-1],"]",sep="")
pyramid.plot(100*denomm$counts/sum(denomm$counts),100*denomf$counts/sum(denomf$counts), labels=lbz,gap=.15*max(stuff),main="Demography")
dev.off()
## lxcol="black",rxcol="black",

## ;;;;;;;;;;;;;;;TEXT OUTPUT
## HH RR for HIV
PSSr2 <- PSS[ PSS[,"age"]>=15 & PSS[,"age"]<50,]
prevh <- mean(PSSr2[,"HIV"])
prevhhh <- mean(PSSr2[PSSr2[,"HHHIV"]==1,"HIV"])
print(paste("HIV prevalece=",format(prevh*100,digits=2),"%",sep=""))
print(paste("HIV prev in HIV hhs=",format(100*prevhhh,digits=2),"%",sep=""))
print(paste("crude HHHIV RR=",format(prevhhh/prevh,digits=2),sep=""))


## OR for prevalent TB from HIV
hp <- PSSr2[,"HIV"]==1
hn <- PSSr2[,"HIV"]==0
artp <- PSSr2[,"ART"]==1
tbu <- PSSr2[,"TBu"]==1
prevhp <- mean(tbu[hp])
prevhn <- mean(tbu[hn])
cOR <- ( prevhp/(1-prevhp) ) / ( prevhn/(1-prevhn) ) 

## adjusted OR
newdf <- cbind(PSSr2[,c("TBu","HIV","gender","age")],as.numeric(PSSr2[,"age"])^2)
names(newdf) <- c(names(newdf)[-5],"age2")
aamod <- glm(TBu ~ HIV + gender + age + age2, family = binomial, data=newdf)
aOR <- exp(coef(aamod)["HIV"])

print(paste("Prevalent TB odds-ratio given HIV: cOR=",
            format(cOR,digits=2),"; aOR=",
            format(aOR,digits=2),sep=""))

## ART in HIV...
artprev <- sum(hp*artp)/sum(hp)
print(paste("ART among HIV+ = ",format(100*artprev,digits=2),"%",sep=""))

## mortality among ART defaulters
## artdefdead <- mean(PSS[PSS[,'artdefr']==1,'isX'])
## print(paste("mortality among ART defaulters = ",format(100*artdefdead,digits=2),"%",sep=""))

