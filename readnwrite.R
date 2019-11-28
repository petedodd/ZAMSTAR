## a function for writing out a new parz.dat file
writeparz <- function(Intfile,Outputloc){
  
  para<-Intfile

  pp<-"//run parameters: noruns, starttime, stoptime, dt, S0,L0,D0,T0;"
  cat(pp,'\n',file=Outputloc)
  nn<-8 # Number at end of first row
  cat(para[1:nn,2],sep='\t',file=Outputloc,append=TRUE)

  pp<-"//population parameters: fertage_alpha, fertage_beta, mf_deathratio, AGEflag;"
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  i<-4 # Number at end of first row
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
    
  pp<-"//non-HIV TB parameters:beta, tbsurvn_k, tbsurvn_L, tbmortn, r0,r1,smrfac, fscale, smrnmort, CDRfile;"
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  i<-10 # Number at end of first row
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
    
  pp<-"//more non-HIV TB: CDRn,RxL,Rxp,Rxpd,Iprotn,IprotnS,Pprotn,PprotnS, fSmr, deltanSmr, tbd_betan;"
  i<-11 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
    
  pp<- "//HIV parameters:thiv_alpha,thiv_beta,mhiv_alpha,mhiv_beta,fhiv_alpha,fhiv_beta,fm_hivratio, cd40;"
  i<-8 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i

  pp<-"//HIV incidence parameters:h_t0,h_peak,h_peakiness, h_peaktime, h_theta, decline2010, hincF, hincFN;"
  i<-8 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
  
  pp<-"//HIV-TB parameters: hivdurf, f,g,rho,CDRp,Rxpp, Rxppd,deltapSmr, tbd_betap, H2Lhr;"
  i<-10 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
  
  pp<-"//household parameters: hhFLAG, nophh0, hhmoverate, betaH, hhhivRR;"
  i<-5 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
  
  pp<-"//intervention parameters:ECFflag, ECFst, ECFet, CDRn2, CDRp2, tbd_betan2, tbd_betap2;"
  i<-7 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
  
  pp<-"//Rx improvement: Rxflag, Rxtime, Rxp2, Rxpp2;"
  i<-4 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
    
  pp<-"//intervention parameters:IPTflag, IPTst, IPTet, IPTcov, IPTT, IPTart, IPTshiv, IPTcovP, IPTcovP2, IPTcovPt1, IPTcovPt2, TSTspec, TSTsens;"
  i<-13 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
  
  pp<-"//IPT parameters: IPTpprot, IPTdurnN, IPTdurnP, IPThrN, IPThrP, IPThrA, IPTcprobN, IPTcpropP, IPTmultiple;"
  i<-9 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i

  pp<-"//intervention parameters: ARTflag, ARTst, ARTrt, ARTcov, tbart, cd4h, cd4g, cd4g2, ARTf, cd42a, cd42b, cd42c, expET, artTk, artTl, dCD4st, dCD4et, dCD4ep, art750, artDR, artCU ;"
  i<-19+2 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
  
  pp<-"//intervention parameters:HHIflag, HIPflag, HARflag, HHIst, HHIet, HHIcov, TSTu16sens, hhdOR, hhdF, deltaSmrI; "
  i<-10 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
  
  pp<-"//intervention parameters: tECFflag, tECFst, tECFet, tecfdOR, tecfdF, tECFhaz ;"
  i<-6 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
  
  pp<-"//external data filenames: birthFN, deathFN, migrateFN, amFN, popfile, LTbFN, abFN, rpopFN;"
  i<-8 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
  
  pp<-"//additional ART modifications:ART2TB, ART2TBst, ART2HHhx, ARThxT, shivt, shivp ;"
  i<-6 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
  
  pp<-"//TnT parameters: massart, artstart, artT, artP, partrefusenik;"
  i<-5 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
  
  pp<-"//HIV screening parameters: HIVflag, HIVst, HIVet, HIVT, HIVeff;"
  i<-5 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
  
  pp<-"//PACF parameters: PACFlag, PACFst, PACFet, PACFT, PACFeffsp, PACFeffsn, PACFeffart, PACFhivOR;"
  i<-8 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
  
  pp<-"//resistance parameters: Rflag, Rfrac0;"
  i<-2 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
  
  pp<-"//dynamic HIV: DHflag, sbFracH, sbFracM, sbc1, sbc2, sbc3, sba1, sba2, sba3, sbA, lam, riA, riE, durA, durE, ARTeps, hiv0 ;"
  i<-17 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
  
  pp<-"//HIVMAC parms: hivmac, artsig, pscale, pn, pmaxval, kq, gg, gi, glrst, nglstr, hseed, SQ, SQGLprop;"
  i<-13 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i
  
  pp<-"//switches:qlog, light, snaps, bigSS;"
  i<-4 # Number at end of first row
  cat('\n',pp,'\n',file=Outputloc,append=TRUE)
  cat(para[(nn+1):(nn+i),2],sep='\t',file=Outputloc,append=TRUE);nn <- nn+i

}





trim <- function (x) gsub("^\\s+|\\s+$", "", x)
## a function that converts parz.dats to input.txts
convertdat2txt <- function(datfn,txtfn){
  nolines <- 48
  con <- file(datfn)
  psa <- readLines(con)
  close(con)
  tmpm <- c('parm\tvalue')
  for (i in seq(from=1,to=nolines,by=2)){
    tmp <- strsplit(psa[i],':')[[1]][2]
    tmp <- strsplit(tmp,';')[[1]][1]
    tmp <- strsplit(tmp,',')[[1]]
    tmp <- trim(tmp)                     #the list of names
    tmpv <- strsplit(psa[i+1],'\\s+')[[1]] #the values
    n <-length(tmp) 
    if(n!=length(tmpv)){print(i);print(tmp);print(tmpv);}
    tmpmtmp <- paste(tmp,rep('\t',n),tmpv,sep='')
    tmpm <- c(tmpm,tmpmtmp)
  }
  con <- file(txtfn)
  writeLines(tmpm,con)
  close(con)
}

## convertdat2txt('../data/parz.dat','../data/ZOrig.txt') #generate a base set of parameters in txt form
## para <- read.table('../data/ZOrig.txt', header = TRUE, stringsAsFactors = FALSE )
## writeparz(para,'pSCtest.dat')

## this reads in a txt tab-file and outputs a  useable parz file
## changes must be a named vector/list
changeparz <- function(basetxtfile,outdatfile,changes){
    input <- read.table(basetxtfile, header = TRUE, stringsAsFactors = FALSE )
    if( length(changes) > 0){
        row.names(input) <- input[,1]
        for(nm in names(changes)){
            if(! nm %in% input[,1])stop(paste0('Trying to change parameter not there! i.e. ',nm))
            input[nm,2] <- changes[nm]        #change the inputs
        }
    }
    writeparz(input,outdatfile)  
}

load('shared/Stargs.Rdata')
load('shared/Ztargs.Rdata')

## function to collect the target errors
getErr <- function(country='ZM'){
    setwd('data')
    ## check results are there!
    E1 <- file.exists(file='results1.dat')
    E2 <- file.exists(file='results2.dat')
    E3 <- file.exists(file='results3.dat')
    E4 <- file.exists(file='results4.dat')
    E5 <- file.exists(file='results5.dat')
    if( ! E1&E2&E3&E4&E5) stop('results missing!')                     #all ok
    ## read in relevant bits
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
    setwd('..')
    ## computations
    if(country=='SA')                      #choose country
        targs <- Stargs
    else
        targs <- Ztargs
    SSE <- rep(0,9)                 #vector of errors
    yr <- 1990:2010                     #years
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
    SSE
}

## a function to change parameter, run the model and return some error statistics
ChangeAndRun <- function(changes){
    ## --- make the changes
    changeparz('data/ZOrig.txt','data/parz.dat',changes=changes)
    ## --- run the model
    system('./TB data/')                #Un*x
    ## system('TB.exe data/')                #windows
    ## --- get the target errors
    getErr('ZM')
}
