## THIS VERSION IS JUST TO DUMP THE HISTORY AND DO (IN CHUNKS) AN I.S.
## called like R --vanilla --slave --args <cAMIS.R $iR --vanilla --slave --args <cAMIS.R $i
## looking at new set of parms and flat priors
## this is going to just use sequential importance sampling

## ## TESTING
## setwd("/Users/pjd/Documents/code/TBibm/ZSclean/codata/fit/working")

test <- commandArgs(trailingOnly=T)
if(length(test)!=1){stop("require 1 argument,k!")}

k <- as.numeric(test)
library(lhs)
source('readnwrite.R')
##library(lattice)

## ;-------------------utilities

## change to use the geterrors function and capture time-series
writeandgo <- function(country,bdr,j){
  ## print(j)
  setwd(paste(bdr,j,sep=''))
  E1 <- file.exists(file='results1.dat')
  E2 <- file.exists(file='results2.dat')
  E3 <- file.exists(file='results3.dat')
  E4 <- file.exists(file='results4.dat')
  E5 <- file.exists(file='results5.dat')

  if(E1){
    R1 <- read.table(file='results1.dat',sep='\t',header=FALSE) #has THIV is 6
    hmf <- R1[,8]
  }
  if(E2){
    R2 <- read.table(file='results2.dat',sep='\t',header=FALSE) #has HIV in Inc(7)
    hivinITB <- mean(R2[R2[,1]>=2009,7],na.rm=T)
  }
  if(E3){
    R3 <- read.table(file='results3.dat',sep='\t',header=FALSE) #has THIV is 6
    hcopm <- R3[,11]
    tbcopm <- R3[,10]*1e5
  }
  if(E4) {
    R4 <- read.table(file='results4.dat',sep='\t',header=FALSE)  #has incidence
    ## 2 is prevalence; 5 is incidence
    prev <- mean(R4[R4[,1]>=2009,2],na.rm=T)
  }
  if(E5) {
    R5 <- read.table(file='results5.dat',sep='\t',header=FALSE)
  }

  if( E1&E2&E3&E4&E5){                     #all ok
    if(country=='SA')                      #choose country
      targs <- Stargs
    else
      targs <- Ztargs
    SSE <- rep(0,9)
    yr <- 1990:2010
    ht <- targs$whohiv                 #HIV prevalence in over 18s
    hm <- ht
    for(i in 1:length(hm)){
      hm[i] <- 100*hmf[which(R1[,1]==yr[i])]
    }
    tn <- length(ht)
    names(SSE) <- c('hiv','UNhiv','art','hhhiv','tbp','ari','hhtb','rx','hivrx')
    SSE[1] <- ((hm[tn]-ht[tn])/ht[tn]) / targs$zshivw # the ZAMSTAR HIV bit
    SSE[2] <- mean((hm-ht)^2/ht^2) / (targs$whohivw^2) #WHO HIV bit - corrected for length too
    SSE[3] <- ((100*R5[dim(R5)[1],7] - targs$zsart)/targs$zsart) / targs$zsartw #ZS ART coverage
    SSE[4] <- ((hcopm[length(hcopm)]/hmf[length(hcopm)] - targs$hhhivOR)/targs$hhhivOR) / targs$hhhivORw #hhhiv
    SSE[5] <- ((R4[dim(R4)[1],2]  - targs$zstbp)/targs$zstbp) / targs$zstbpw #TB prevalence
    SSE[6] <- ((100*R3[dim(R3)[1],3]  - targs$zsari)/targs$zsari) / targs$zsariw #ARi
    ## skipping tbhiv as this is not currently outputted?
    tbr <- tbcopm/R4[,2]                               #ratio
    tbor <- mean(tbr[(dim(R4)[1]-100):dim(R4)[1]])     #last 10 years
    SSE[7] <- ((tbor - targs$zshhtbOR)/tbor) / targs$zshhtbORw                 #hh/TB OR
    SSE[8] <- ((R4[dim(R4)[1],4] - targs$zsrx)/targs$zsrx) / targs$zsrxw #Rx
    hivinrx <- R1[,6]
    Lstyr <- (length(hmf)-10):length(hmf)
    SSE[9] <- mean( (hivinrx/(1-hivinrx))[Lstyr], na.rm=TRUE ) / mean( (hmf/(1-hmf))[Lstyr], na.rm=TRUE )
    hrxor <- ( hivinrx/(1-hivinrx) ) / ( hmf/(1-hmf) )
    hrxor[is.nan(hrxor)] <- NA
    ## SSE[9] <- ( (mean(hrxor[(length(hrxor)-10):length(hrxor)],na.rm=TRUE)  - targs$zsrxhivOR)/targs$zsrxhivOR)^2 / targs$zsrxhivORw^2
    inc <- R4[,5]                      #incidence
    prev <- R4[,2]                      #prevalence
    hii <- R2[,7]                      #HIV in Inc
    hiv <- hmf                         #HIV
    art <- R5[,7]                      #ART
    rx <- R4[,4]                       #Rx
    ari <- R3[,3]                      #ARI
    SSE[is.nan(SSE)] <- NA            #if HIV zero  etc...
    SSE[abs(SSE)==Inf] <- NA            #if HIV zero  etc...
    SSE[abs(SSE)>1e6] <- NA            #may as well be
    ## interventions
    interv <- R5[,1:5]; colnames(interv) <- c('time','PHH','HHART','HHIPT','ECF')
  } else {
    print(  paste(bdr,j,' missing results',sep='')) 
    SSE <- rep(NA,9)
    inc <- prev <- hii <- hiv <- art <- rx <- ari <- hivinrx <- rep(NA,300) ;
    tbor <- NA; hrxor <- rep(0,300); hcopm <- rep(0,300) ; hmf <- rep(1,300) ; interv <- matrix(NA,ncol=5,nrow=300)
  }
  setwd('..')
  ans <- list(SSE=SSE,inc=inc,prev=prev,hii=hii,hiv=hiv,art=art,rx=rx,ari=ari,hivinrx=hivinrx,
              c(hcopm[length(hcopm)]/hmf[length(hcopm)],
                tbor,mean(hrxor[(length(hrxor)-10):length(hrxor)])), interv=interv)
  return(ans)
}
                                        #taken from NM stuff


## doesn't appear to need changing
tidy <- function(i){
  if( i > 0 ) {
    system('mkdir tmp')
    system('mv */*.pdf tmp')
    setwd('tmp')
    system(paste('pdftk *.pdf cat output ALL',i,'.pdf',sep=''))
    system(paste('mv ALL',i,'.pdf ../',sep=''))
    setwd('..')
    system('rm -r tmp')
  }
  system('rm -r work*')
}

## 
writeoutINT <- function( parz, wkdir, IFname, country ){
  ## if(FALSE){ pnmz <- colnames(parz); parz <- cbind(5,parz);colnames(parz) <- c('noruns',pnmz);} #multiple runs!!!ls
  setwd( paste0(wkdir,'/',country) )
  olddirs <- list.files(pattern='work*',include.dirs=T)
  intdirs <- list.files(pattern='INT*',include.dirs=T)
  if(length(olddirs)>0){
    for(dr in olddirs){system(paste('rm -r',dr))}
    for(dr in intdirs){system(paste('rm -r',dr))}
  }
  bsnmz <- colnames(parz)
  extras <- c('HHIflag', 'hhdOR', 'hhdF','HHIcov','HIPflag','HARflag') #for interventions
  parzA <- parzB <- parzC <- parzD <- cbind(parz,matrix(1,ncol=length(extras),nrow=nrow(parz)))
  colnames(parzA) <- colnames(parzB) <- colnames(parzC) <- colnames(parzD) <- c(bsnmz,extras)
  extras <- c(extras,'tECFflag','tecfdOR','tecfdF','tECFhaz') #including also diffusion
  parzE <- parzF <- cbind(parz,matrix(1,ncol=length(extras),nrow=nrow(parz)))
  colnames(parzE) <- colnames(parzF) <- c(bsnmz,extras)

  ## A - HH no BC
  ## B - HH 12
  ## C - HH 21
  ## D - HH 1010
  parzA[,'HHIcov'] <- parzB[,'HHIcov'] <- parzC[,'HHIcov'] <- parzD[,'HHIcov'] <- parzE[,'HHIcov'] <- parzF[,'HHIcov'] <- 1
  parzA[,'HIPflag'] <- parzB[,'HIPflag'] <- parzC[,'HIPflag'] <- parzD[,'HIPflag'] <- parzE[,'HIPflag'] <- parzF[,'HIPflag'] <- 1
  parzA[,'HARflag'] <- parzB[,'HARflag'] <- parzC[,'HARflag'] <- parzD[,'HARflag'] <- parzE[,'HARflag'] <- parzF[,'HARflag'] <- 1
  parzA[,'HHIflag'] <- parzB[,'HHIflag'] <- parzC[,'HHIflag'] <- parzD[,'HHIflag'] <- parzE[,'HHIflag'] <- parzF[,'HHIflag'] <- 1
  parzA[,'hhdOR'] <- 1; parzB[,'hhdOR'] <- 1; parzC[,'hhdOR'] <- 2; parzD[,'hhdOR'] <- 10
  parzA[,'hhdF'] <- 1;  parzB[,'hhdF'] <- .5; parzC[,'hhdF'] <- 1;  parzD[,'hhdF'] <- .1
  ## diffusion
  ## E - HH22 and beta 5
  ## F - HH22 and beat 10
  parzE[,'tECFflag'] <- parzF[,'tECFflag'] <- 1  #on
  parzE[,'hhdOR'] <- parzF[,'hhdOR'] <- 1.5  #sure
  parzE[,'hhdF'] <- parzF[,'hhdF'] <- 0.67   #swift
  parzE[,'tecfdOR'] <- parzF[,'tecfdOR'] <- 1.5  #sure
  parzE[,'tecfdF'] <- parzF[,'tecfdF'] <- 0.67   #swift
  parzE[,'tECFhaz'] <- 2; parzF[,'tECFhaz'] <- 5 #beta
  for(i in 1:L){
    system(paste('mkdir work',i,sep=''))
    system(paste('mkdir INTA',i,sep=''))
    system(paste('mkdir INTB',i,sep=''))
    system(paste('mkdir INTC',i,sep=''))
    system(paste('mkdir INTD',i,sep=''))
    system(paste('mkdir INTE',i,sep=''))
    system(paste('mkdir INTF',i,sep=''))
  }
  nruns <- 1
  for( i in 1:L ){
    nm <- paste('work',i,'/parz.dat',sep="")
    changeparz(IFname,nm,changes=parz[i,])
    nm <- paste('INTA',i,'/parz.dat',sep="")
    changeparz(IFname,nm,changes=parzA[i,])
    nm <- paste('INTB',i,'/parz.dat',sep="")
    changeparz(IFname,nm,changes=parzB[i,])
    nm <- paste('INTC',i,'/parz.dat',sep="")
    changeparz(IFname,nm,changes=parzC[i,])
    nm <- paste('INTD',i,'/parz.dat',sep="")
    changeparz(IFname,nm,changes=parzD[i,])
    nm <- paste('INTE',i,'/parz.dat',sep="")
    changeparz(IFname,nm,changes=parzE[i,])
    nm <- paste('INTF',i,'/parz.dat',sep="")
    changeparz(IFname,nm,changes=parzF[i,])

  }
  setwd( wkdir )
}


## plotting prior and LHSing data
realpar <- function(parname,paramsi,parz){
  parno <- which(paramsi$names==parname)
  top <- paramsi$top[parno]; bot <- paramsi$bot[parno];
  s1 <- paramsi$s1[parno]; s2 <- paramsi$s2[parno];
  xz <- seq(from=0.75*bot,to=1.25*top,length=500);
  yz <- dbeta((xz-bot)/(top-bot),shape1=s1,shape2=s2)
  ans <- bot+(top-bot)*qbeta(parz[,parname],shape1=s1,shape2=s2)
  return(ans)
}


## changed to use wkdir and dadir
firstround <- function(BC,BS,wkdir,country){
  setwd(wkdir)
  ## parameter nos
  nbc <- length(BC$names)
  nbs <- length(BS$names)

  ## to get round country prob
  if(!file.exists('bgLHS.Rdata')){
    bgsample <- randomLHS(n = L * nosteps, k=nbc+nbs ) #make sample
    colnames(bgsample) <- c(BC$names,BS$names)
    save(file = 'bgsample.Rdata',bgsample)   #NB this is above the country level & raw
  } else {
    load('bgsample.Rdata')
  }

  for(i in nosteps:1){
    tmp <- bgsample[1:L + (i-1)*L,]
    ## make real
    for( nm in BC$names ){    #make real
      tmp[,nm] <- realpar(parname=nm,paramsi=BC,parz=tmp)
    }
    for( nm in BS$names ){    #make real
      tmp[,nm] <- realpar(parname=nm,paramsi=BS,parz=tmp)
    }
    parz <- tmp
    save(file=paste0(country,'/parz',i,'.Rdata'),parz)
  }

  ## -----------prep
  writeoutINT( parz=parz , wkdir=wkdir, IFname='SInfClean.txt', country=country )
  ## run scripts
  dnames <- c('work',paste0('INT',c('A','B','C','D','E','F')))
  ## dnames <- c('work',paste0('INT',c('E','F')))
  ## dnames <- c('work',paste0('INT',c('D','E','F')))
  for(wdn in dnames){
    sas <- paste0("#!/bin/bash\n#$ -cwd -V\n#$ -l h_rt=00:02:00\n#$  -l h_rss=1G,h_vmem=1.2G\n")
    sas <- c(sas,paste('#$ -t 1-',L,'\n',sep=''))
    binstr <- paste('cd ',wkdir,'/',country,'\n./TB ',wdn,'${SGE_TASK_ID}/',sep='')
    sas <- c(sas,binstr)
    con <- file(paste0(country,wdn,'runs.sh'))
    writeLines(sas,con)
    close(con)
    system(paste0('chmod +x ',country,wdn,'runs.sh'))
  }
  return(head(parz))
}



## contents:
## calculate the weights of the current generation
## update the weights of the previous generation
## tidy the working file structure
## estimate the new proposal parms
## write out the work data



getanalysepropose <- function(wkdir,k){
  ## load the parameters
  load('BC.Rdata'); load('BSSA.Rdata'); load('BSZM.Rdata')

  ## country loop to read in the results
  for(country in c('SA','ZM')){
    setwd(paste0(wkdir,'/',country))
    res <- resA <- resB <- resC <- resD <- resE <- resF <- list()
    errz <- matrix(NA,nrow=L,ncol=9)
    for (i in 1:L){
      res[[i]] <- writeandgo(country,'work',i)         #change to correct errors and folder
      errz[i,] <- res[[i]]$SSE^2
      errz[i,2] <- sqrt(errz[i,2])
      resA[[i]] <- writeandgo(country,'INTA',i)         #change to correct errors and folder
      resB[[i]] <- writeandgo(country,'INTB',i)         #change to correct errors and folder
      resC[[i]] <- writeandgo(country,'INTC',i)         #change to correct errors and folder
      resD[[i]] <- writeandgo(country,'INTD',i)         #change to correct errors and folder
      resE[[i]] <- writeandgo(country,'INTE',i)         #change to correct errors and folder
      resF[[i]] <- writeandgo(country,'INTF',i)         #change to correct errors and folder
    }
    save(file=paste0('results',k,'.Rdata'),res)    #this is going to need some post-processing...
    save(file=paste0('resultsA',k,'.Rdata'),resA)    #this is going to need some post-processing...
    save(file=paste0('resultsB',k,'.Rdata'),resB)    #this is going to need some post-processing...
    save(file=paste0('resultsC',k,'.Rdata'),resC)    #this is going to need some post-processing...
    save(file=paste0('resultsD',k,'.Rdata'),resD)    #this is going to need some post-processing...
    save(file=paste0('resultsE',k,'.Rdata'),resE)    #this is going to need some post-processing...
    save(file=paste0('resultsF',k,'.Rdata'),resF)    #this is going to need some post-processing...
    colnames(errz) <- names(res[[1]]$SSE)

    if(country=='SA'){                  #record
      SAerrz <- errz
    } else {
      ZMerrz <- errz
    }
    setwd('..')
  }

  
  ## tidying
  system('rm -r SA/work*')
  system('rm -r ZM/work*')
  system('rm -r SA/INT*')
  system('rm -r ZM/INT*')

  load(paste0('SA/parz',k,'.Rdata')); SAparz <- parz
  load(paste0('ZM/parz',k,'.Rdata')); ZMparz <- parz
  ## ;;;;;;;;;;;;;;; preparation ---------------
  writeoutINT( parz=SAparz , wkdir=wkdir, IFname='SInfClean.txt', country='SA' )
  writeoutINT( parz=ZMparz , wkdir=wkdir, IFname='SInfClean.txt', country='ZM' )
  return(head(parz))
}


## -----------testing work--------------
## ----------MAIN---------------------------
## ## cluster change
## wdir <- "/users/eidepdod/Documents/SArf1"
## ## wdir <- "~/Documents/code/TBibm/aannecy/R4h/fitting/fit3"
##  ddir <- '/users/eidepdod/Documents/SArf1'
## ## ddir <- "~/Documents/code/TBibm/aannecy/R4h/fitting/fit3"
## setwd(ddir)
## #N <- read.csv(file='SAtotalnotes.csv')
## P <- read.csv(file='TBtargetsSA.csv')

## wdir <- '/Users/pjd/Documents/code/TBibm/ZSclean/codata/fit/testing'
wdir <- '/home/cm1pjd/ZAMSTAR/intsmear'
## wdir <- "/Users/pjd/Documents/code/TBibm/ZSclean/codata/fit/testing"
setwd(wdir)

## how many samples
##L <- 2000
## L <- 10
## L <- 6
L <- 250
nosteps <- 10


load('shared/Ztargs.Rdata');
load('shared/Stargs.Rdata');


## new
## SA - reorganise into all and specific! check meanings and whether necessary
## chec f, g. deltapSmr
## HIV mortality to lock in above!
## f smr reduction if Lhiv - HR; smrfac -overall bodge in smear positivity
## fscale scales the fast risk from EV
## deltanSmr - HR of detection if smr- (nb cdr is smr+)

## CLASSES of parm -
## background common : BC
## background specific 2 : BS
## fit common :FC
## fit specific 2: FS


## NB will need to have changed SInfClean to accomodate new filestructure


## FC
FCnmz <- c('Pprotn','r0','fscale')
FCvlz <- c(0.391691635070732,0.0090966668965246,1.0)
FCbot <- c(0.2,0.0002,0.8); FCtop <- c(0.9,0.006,1.2)
FCs1 <- c(1,1,1);FCs2 <- c(1,1,1);
## FCs1 <- c(2,4);FCs2 <- c(2,16);
FC <- list(names=FCnmz,bot=FCbot,top=FCtop,s1=FCs1,s2=FCs2)


## FS
Fitnmz <- c('beta','h_peak','tbd_betan','tbd_betap','betaH',  'hhhivRR','pmaxval')
safv <- c(6.68186713015108,2.50302685365914,0.405820022937214,0.545160918857013,0.713166251784168,6.16085407980799,0.1493)
zmfv <- c(4.76860790096107,1.9,0.223438290551665,0.400144898156517,0.68463986496968,4.34413895672601,0.02300104)
FSbotSA <- safv * 0.9;FStopSA <- safv * 1.1; #20% variation
FSbotZM <- zmfv * 0.9;FStopZM <- zmfv * 1.1; #20% variation
FSs1SA <- FSs2SA <- FSs1ZM <- FSs2ZM <- rep(1,length(Fitnmz))
## now change the detection ones - 3,4
FSbotSA[3:4] <- 0.1;FStopSA[3:4] <- 0.9; 
FSbotZM[3:4] <- 0.1;FStopZM[3:4] <- 0.9; 
FSs1SA[3:4] <- c(1,1); FSs2SA[3:4] <- c(1,1);
FSs1ZM[3:4] <- c(1,1); FSs2ZM[3:4] <- c(1,1);
## FSs1SA[3:4] <- c(4,4); FSs2SA[3:4] <- c(6,4);
## FSs1ZM[3:4] <- c(2,4); FSs2ZM[3:4] <- c(4,6);
## and beta:
FSbotSA[1] <- 5;FStopSA[1] <- 15; #20% variation
FSbotZM[1] <- 5;FStopZM[1] <- 15; #20% variation
## join
FSSA <- list(names=Fitnmz,bot=FSbotSA,top=FStopSA,s1=FSs1SA,s2=FSs2SA)
FSZM <- list(names=Fitnmz,bot=FSbotZM,top=FStopZM,s1=FSs1ZM,s2=FSs2ZM)



## BC
## BCnmz <- c('tbsurvn_k','tbsurvn_L','tbmortn','r1','f','smrfac','fscale','smrnmort','fSmr','deltanSmr','thiv_alpha','thiv_beta','hivdurf','rho','hhmoverate')
## BCvlz <- c(1.5,3.25,0.7,0.067,0.45,1.0,1.0,0.3,0.23,0.7,2.3,13.3,0.1,0.36,0.1)
BCnmz <- c('tbsurvn_k','tbsurvn_L','tbmortn','r1','f','smrfac', 'smrnmort','fSmr','deltanSmr','thiv_alpha','thiv_beta','hivdurf','rho','hhmoverate')
BCvlz <- c(1.5,3.25,0.7,0.067,0.45,1.0,  0.3,0.23,0.7,2.3,13.3,0.1,0.36,0.1)
BCbot <- BCvlz * 0.975 # NB these need checking to make sure meaningful
BCtop <- BCvlz * 1.025
BCs1 <- BCs2 <- rep(1,length(BCvlz))
BC <- list(names=BCnmz,top=BCtop,bot=BCbot,s1=BCs1,s2=BCs2)
if(!file.exists('BC.Rdata')) save(file='BC.Rdata',BC)
print('BC done!')

## BS
BSnmz <- c('D0','CDRn','RxL','Rxp','Rxpd','h_peakiness','h_peaktime','h_theta','CDRp','Rxpp','Rxppd','pscale','pn')
BSvlzSA <- c(600,0.72,0.5,0.7,0.25,27.1921095091006,1995.66102039002,0.572849904558087,0.72,0.7,0.25,10.21932597,1.868635566)
BSvlzZM <- c(400,0.72,0.5,0.7,0.25,2.5,1990,0.5,0.72,0.7,0.25,3.563549,2.687766)
BSbotSA <- BSvlzSA * 0.975; BSbotZM <- BSvlzZM * 0.975
BStopSA <- BSvlzSA * 1.025; BStopZM <- BSvlzZM * 1.025
BSs1SA <- BSs2SA <- BSs1ZM <- BSs2ZM <- rep(1,length(BSnmz))
BSSA <- list(names=BSnmz,top=BStopSA,bot=BSbotSA,s1=BSs1SA,s2=BSs2SA)
BSZM <- list(names=BSnmz,top=BStopZM,bot=BSbotZM,s1=BSs1ZM,s2=BSs2ZM)
## different for h_peaktime:
BSSA$bot[7] <- BSvlzSA[7]-3; BSSA$top[7] <- BSvlzSA[7]+3;
BSZM$bot[7] <- BSvlzZM[7]-3; BSZM$top[7] <- BSvlzZM[7]+3;

## join with old fit parms - all country specific
BSSA$names <- c(BSSA$names,FC$names,FSSA$names); BSZM$names <- c(BSZM$names,FC$names,FSZM$names);
BSSA$top <- c(BSSA$top,FC$top,FSSA$top); BSZM$top <- c(BSZM$top,FC$top,FSZM$top);
BSSA$bot <- c(BSSA$bot,FC$bot,FSSA$bot);  BSZM$bot <- c(BSZM$bot,FC$bot,FSZM$bot); 
BSSA$s1 <- c(BSSA$s1,FC$s1,FSSA$s1); BSZM$s1 <- c(BSZM$s1,FC$s1,FSZM$s1);
BSSA$s2 <- c(BSSA$s2,FC$s2,FSSA$s2);  BSZM$s2 <- c(BSZM$s2,FC$s2,FSZM$s2); 

if(!file.exists('BSSA.Rdata'))save(file='BSSA.Rdata',BSSA)
if(!file.exists('BSZM.Rdata'))save(file='BSZM.Rdata',BSZM)
print('BS done!')



## -----------------MAIN------------------
if( k==0 ){
  ## the initial writings out
  firstround( BC=BC,BS=BSSA,wkdir=wdir, country='SA' )
  firstround( BC=BC,BS=BSZM,wkdir=wdir, country='ZM' )
} else {
  ## ## run it - on PC
  ## setwd('SA')
  ## for(i in 1:L) system(paste0('./TB work',i,'/'))
  ## setwd('../ZM')
  ## for(i in 1:L) system(paste0('./TB work',i,'/'))
  ## setwd('..')

  ## analyse
  getanalysepropose(wdir, k)
}

