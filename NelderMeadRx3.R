## this is going to be a Nelder-Mead fitting routine, aimed at approaching the HIV parameters

## generic NM routine:


## make initial simplex
initSimplex <- function( x0, lb, ub, w=.1){ #'width' as a propn of width
    np <- length(lb)
    X <- matrix(0,nrow=(np+1),ncol=(np+1))       #each row a point, each col a par
    for(i in 1:np){
        X[,i] <- runif(np+1) * w * (ub[i]-lb[i]) + lb[i]
    }
    X[,(np+1)] <- NA                    #the f(x)
    return(X)
}

## ## trial function
## F <- function(x){
##     return(sum((x-c(0.1,0.2,0.2,0.3,0.5))^2))
## }

## volume of simplex
simpVol <- function(X,np){
    Y <- X[1:np,1:np]
    Y <- Y-matrix(X[np+1,1:np],ncol=np,nrow=np,byrow = T)
    return(det(Y))
}

## fractional volume  made one-d
simpSize <- function(X,lb,ub){
    ans <- simpVol(X,length(lb)) / prod(ub-lb)
    ans <- abs(ans)
    return(ans^(1/length(lb)))
}

## pointwrapping
internalize <- function(X,lb,ub){
    for(i in 1:length(lb)){                     #a matrix
        if(!is.null(dim(X)[1])){
            X[X[,i]>ub[i],i] <- (X[X[,i]>ub[i],i] - lb[i])%%(ub[i]-lb[i]) + lb[i] #wrapping exceeders
            X[X[,i]<lb[i],i] <- ub[i] - (lb[i] - X[X[,i]<lb[i],i])%%(ub[i]-lb[i]) #wrapping undershoots
        } else {                        #a vector
            if(X[i]>ub[i]){
                X[i] <- (X[i] - lb[i])%%(ub[i]-lb[i]) + lb[i] #wrapping exceeders
            }
            if(X[i]<lb[i]){
                X[i] <- ub[i] - (lb[i] - X[i])%%(ub[i]-lb[i]) #wrapping undershoots
            }
        }
    }
    return(X)
}

## counter--
inc <- function(){
    assign("nstep", nstep+1, envir = .GlobalEnv)
}




## ---------------stepping function--------------------
NMstep <- function(X,F,lb,ub,tol=1e-10,tol2=1e-2,maxstep=1000,alph=1,gam=2,rho=-.5,sig=.5){
    np <- length(lb)
    fxr <- 10
    while( nstep < maxstep  & simpSize(X,lb=lb,ub=ub) > tol & fxr > tol2){
        print(paste('nstep:',nstep,sep=''))
        X <- X[order(X[,np+1]),] ## order
        xo <- colMeans(X[1:np,1:np]) ## CofM
        print(c('xo=',paste(xo)))
        xr <- xo + alph*(xo-X[np+1,1:np])        ## reflect
        xr <- internalize(xr,lb=lb,ub=ub)                    #wrap
        fxr <- F(xr); inc()
        ## history <<- rbind(history,colMeans(X[,1:np]))
        if( fxr >= X[1,np+1] & fxr < X[np,np+1]){
            X[np+1,] <-c(xr,fxr)
        } else {
            if ( fxr < X[1,np+1]){               #best so far
                xe <- xo + gam*(xo-X[np+1,1:np])  ## expand
                xe <- internalize(xe,lb=lb,ub=ub)                    #wrap
                fxe <- F(xe); inc()
                if( fxe < fxr ){
                    X[np+1,] <- c(xe,fxe)
                } else {
                    X[np+1,] <- c(xr,fxr)
                }
            } else {                            #fxr>fxn
                xc <- xo + rho*(xo-X[np+1,1:np]) #contract
                xc <- internalize(xc,lb=lb,ub=ub)                    #wrap
                fxc <- F(xc); inc()
                if( fxc < X[np+1,np+1] ){
                    X[np+1,] <- c(xc,fxc)
                } else {
                    for(i in 2:(np+1)){         #reduce
                        X[i,1:np] <- X[1,1:np] + sig*(X[i,1:np]-X[1,1:np])
                        X[i,1:np] <- internalize(X[i,1:np],lb=lb,ub=ub)                    #wrap
                        X[i,np+1] <- F(X[i,1:np]); inc()
                    }
                }                               #end reduce
            }                                   #end contract or reduce
        }                                       #end contract, expand or reduce
        X <- NMstep(X,F,lb=lb,ub=ub, tol=tol,maxstep=maxstep,alph=alph,gam=gam,rho=rho,sig=sig)
    }                                   #end while
    return(X)
}


## ------- to the steps
NMrun <- function(F,lb,ub,w=.5,tol=1e-10,tol2=1e-2,maxstep=1000,alph=1,gam=2,rho=-.5,sig=.5){
    np <- length(lb)
    X <- initSimplex(x0,lb=lb,ub=ub,w=w) ## initiate
    for(i in 1:(np+1)){
        X[i,(np+1)] <- F(X[i,1:np]) ## get values
    }
    assign("nstep", 0, envir = .GlobalEnv)
    X <- NMstep(X,F,lb=lb,ub=ub,tol=tol,maxstep=maxstep ,alph=alph,gam=gam,rho=rho,sig=sig)
    ans <- colMeans(X[,1:np])
    return(list(ans=ans,val=mean(X[,np+1]),err=simpSize(X,lb=lb,ub=ub),nosteps=nstep))
}


source('readnwrite.R')
load('data_rx/Stargs.Rdata')            #NB this is to do with ZAMSTAR
P <- read.csv(file='data_rx/TBtargetsSA.csv')
H <- read.csv(file='../shared/AIDSinfoHIV.csv') #SA national
zth <- H[H$country=='ZA',c('year','hivp','hivplo','hivphi')]
ht <- zth$hivp[1:21]                    #drop 2011
print(ht)

## ub <- c( 20,    1.5e-2,    .9,     .9,    3.2,50,2000,0.9,0.4)   #upper bounds
## lb <- c(  4,   1e-6,   .05,     0.2,    1.8,10,1995,0.4,0.1)  #lower bounds


## HIV only
ubh <- c(  3.2,50,2000,0.9,0.4)   #upper bounds
lbh <- c(   1.8,10,1995,0.4,0.1)  #lower bounds
ubt <- c( 20,    1.5e-2,    1, 0.4)   #upper bounds
lbt <- c(  4,   1e-6,    0.2, 0.3)  #lower bounds


## nph <- 5
## np <- nph
## history <- 0.5*(ubh+lbh)


## ub <- c( 20,    1.5e-2,    1,     1, 0.4,   3.2,50,2000,0.9,0.4)   #upper bounds
## lb <- c(  4,   1e-6,   .05,     0.2, 0.3,   1.8,10,1995,0.4,0.1)  #lower bounds
## np <- length(lb)                        #number of parameters
## nstep <- 0



## function for getting the error given the parameters
getErr <- function(x,err=rep(1,5)){
    x <- c(x[1:3],x[3:length(x)])       #dfrac 2x
    ##x are the parameters to go into the HIV line
    names(x) <- c('beta','r0', 'tbd_betan', 'tbd_betap','Pprotn', 'h_peak','h_peakiness','h_peaktime','h_theta','pmaxval')
    ## write out -  parz=parzo , wkdir=wkdir, IFname='SInfClean.txt' )
    changeparz('data_rx/SInfClean.txt','data_rx/parz.dat',changes=x)
    ## run
    system('./TBrx data_rx/')
    ## read results & compute error
    tmp <- geterrors('data_rx/')
    res <- tmp[[1]]
    res[is.na(res)] <- 0
    names(res) <- c('prev','hivinITB','SSE','SSEp','PSSE','PSSEp','HSSE','HSSEp','whSSE','aSSE')
    fz <- err[1]*res['SSE'] + err[2]*res['whSSE'] + err[3]*res['aSSE'] +  err[4]*res['HSSE']  +  err[5]*res['PSSE']
    print(res)
    print(paste('->',fz))
    ov <- c(x,res['SSE'],res['whSSE'],res['aSSE'],res['HSSE'],res['PSSE'])
    cat(paste(ov,collapse=','),file='output.csv',append=T);cat('\n',file='output.csv',append=T)
    return(fz)
}


## ---------------------functions taken and adapted from SISsa.R
## doesn't appear to need changing
geterrors <- function(bdr){
    setwd(bdr)
    E1 <- file.exists(file='results1.dat')
    E2 <- file.exists(file='results2.dat')
    E3 <- file.exists(file='results3.dat')
    E4 <- file.exists(file='results4.dat')
    E5 <- file.exists(file='results5.dat')

    if(E1){
        R1 <- read.table(file='results1.dat',sep='\t',header=FALSE) #has THIV is 6
        hmf <- R1[,4]                                               #HIV 15-49
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
        SSE <- rep(0,9)
        yr <- 1990:2010

        ## ht <- Stargs$whohiv                 #HIV prevalence in 15-49
        hm <- ht
        for(i in 1:length(hm)){
            hm[i] <- 100*hmf[which(R1[,1]==yr[i])]
        }
        print(hm)
        print(ht)
        tn <- length(ht)

        v1 <- P[,'inc100k']
        v1h <- P[,'inchi'];v1l <- P[,'inclo']
        v2 <- v1; for(i in 1:length(v2)){ v2[i] <- R4[which(R4[,1]==P[i,'year']),5];}
        w1 <- P[,'prev100k']
        w1h <- P[,'prevhi'];w1l <- P[,'prevlo']
        w2 <- w1; for(i in 1:length(w2)){ w2[i] <- R4[which(R4[,1]==P[i,'year']),2];}
        z1 <- P[,'tbhivinc100k']
        z1h <- P[,'tbhinchi'];z1l <- P[,'tbhinclo']
        z2 <- z1;
        for(i in 1:length(z2)){
            z2[i] <- R1[which(R1[,1]==P[i,'year']),6]
        }
        z2 <- z2*v2                           #hivTB inc
        wts2 <- rep(1,length(v1));
        weights <- wts2#wts2^1.5/mean(wts2^1.5)#rep(1,length(v1));
        ## weights<-rep(1,length(v1))*seq(from=1,to=2,along=v1);#
        SSE <- mean(weights*(v1-v2)^2/(v1h-v1l)^2)
        SSEp <- mean(weights*(v1-v2)^2/v1)
        PSSE <- mean(weights*(w1-w2)^2/(w1h-w1l)^2)
        PSSEp <- mean(weights*(w1-w2)^2/w1)
        HSSE <- mean(weights*(z1-z2)^2/(z1h-z1l)^2)
        HSSEp <- mean(weights*(z1-z2)^2/z1)
        whSSE <- mean(weights*(1-hm/ht)^2 / Stargs$whohivw^2) #WHO bit - corrected for length too
        aSSE <- mean(weights)*(1-100*R5[dim(R5)[1],7]/Stargs$zsart)^2 / Stargs$zsartw^2 #Stargs$zsartw^2 #ZS ART coverage
        ans1 <- c(prev,hivinITB,SSE,SSEp,PSSE,PSSEp,HSSE,HSSEp,whSSE,aSSE)
        ans2 <- R4[,5]
        ans3 <- R4[,2]
        ans4 <- R2[,7]
    } else {
        ans1 <- c(0,0,1e10,1e10,1e10,1e10,1e10,1e10,1e10,1e10)
        ans2 <- rep(0,300) ; ans3 <- ans2; ans4 <- ans2;
    }
    setwd('..')
    ans <- list(ans1,ans2,ans3,ans4)

    ## ans <- list(SSE=SSE,inc,prev,hii,hiv,art,rx,ari,
    ##             c(hcopm[length(hcopm)]/hmf[length(hcopm)],
    ##               tbor,mean(hrxor[(length(hrxor)-10):length(hrxor)])))
    return(ans)
}


## -----------------actual work
## first try had c(1,0,0,1,1) for TB error (data 4)
## second try ditching prevalence: c(1,0,0,1,0)
## 3rd try ditching HIVinTB c(1,0,0,0,1)

tbervec <- c(1,0,0,1,1)

## first HIV step
nstep <- 0
mtb <- 0.5*(lbt+ubt)
getErrh <- function(x){
    return(getErr(c(mtb,x),c(0,1,1,0,0)))
}
nmoutH <- NMrun(getErrh,lbh,ubh,w=.5,tol2=0.2,maxstep=100)

## first TBstep
nstep <- 0
getErrtb <- function(x){
    return(getErr(c(x,nmoutH$ans),tbervec))
}
nmoutTB <- NMrun(getErrtb,lbt,ubt,w=.5,maxstep=300)

## second HIV step
lbh <- nmoutH$ans*0.8
ubh <- nmoutH$ans*1.2
lbh[3] <- 1995;ubh[3] <- 2000;
nstep <- 0
getErrh <- function(x){
    return(getErr(c(nmoutTB$ans,x),c(0,1,1,0,0)))
}
nmoutH <- NMrun(getErrh,lbh,ubh,w=.5,tol2=0.2,maxstep=100)

## second TBstep
lbh <- nmoutTB$ans*0.8
ubh <- nmoutTB$ans*1.2
nstep <- 0
getErrtb <- function(x){
    return(getErr(c(x,nmoutH$ans),tbervec))
}
nmoutTB <- NMrun(getErrtb,lbt,ubt,w=.5,maxstep=300)


## 3rd HIV step
lbh <- nmoutH$ans*0.8
ubh <- nmoutH$ans*1.2
lbh[3] <- 1995;ubh[3] <- 2000;
nstep <- 0
getErrh <- function(x){
    return(getErr(c(nmoutTB$ans,x),c(0,1,1,0,0)))
}
nmoutH <- NMrun(getErrh,lbh,ubh,w=.5,tol2=0.2,maxstep=100)

## 4rd TBstep
lbh <- nmoutTB$ans*0.8
ubh <- nmoutTB$ans*1.2
nstep <- 0
getErrtb <- function(x){
    return(getErr(c(x,nmoutH$ans),tbervec))
}
nmoutTB <- NMrun(getErrtb,lbt,ubt,w=.5,maxstep=300)

save(nmoutH,file='nmoutH.Rdata')
save(nmoutTB,file='nmoutTB.Rdata')

