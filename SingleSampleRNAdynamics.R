library(rootSolve)
library(parallel)
library(data.table)


singleSampleRateEstimation  <- function(TPM.data){

    ## some specific helper functions
    refactor <- function(mat){
        if (! is.null(dim(mat))){
            return (list(t(mat[1,]),t(mat[2,])))
        }else{
            return(list(unlist(lapply(mat,"[[",1)),unlist(lapply(mat,"[[",2))))
        }
    }

    refactor2 <- function(mat){return(list(mat[,1],mat[,2]))}
    refactor4 <- function(mat){return(list(mat[1,],mat[2,],mat[3,],mat[4,]))}

    swap.elmt <- function(a1,a2,idx){
        tmp <- a1[idx]
        a1[idx] <- a2[idx]
        a2[idx] <- tmp
        return(list(a1,a2))
    }
    swap <- function(a1,b1,c1,a2,b2,c2){
        idx <- which(a1>a2)
        la <- swap.elmt(a1,a2,idx)
        lb <- swap.elmt(b1,b2,idx)
        lc <- swap.elmt(c1,c2,idx)
        return(list(la[[1]],lb[[1]],lc[[1]],la[[2]],lb[[2]],lc[[2]]))
    }
    EPS = exp(-30)
    data <- data.table(TPM.data)
    data[,pre.frac:=unlabeled.intron/unlabeled.exon]
    data[,lab.frac:=labeled.intron/labeled.exon]
    data[,solvable:=ifelse(lab.frac>1/(2-pre.frac),"C","A")]
    data[pre.frac < EPS,solvable:="L"]
    data[lab.frac>1 | pre.frac >= lab.frac,solvable:="B"]
    data[solvable=="C",c("deg.rate","proc.rate"):=refactor(mapply(solve.rates2,pre.frac,lab.frac))]
    data[solvable=="A",c("deg.rate","proc.rate","deg.rateb","proc.rateb"):=refactor4(mapply(solve.rates3,pre.frac,lab.frac))]
    data[,twosols:= !is.na(deg.rateb)]
    data[solvable=="A" & !twosols, c("deg.rate","proc.rate"):=refactor2(eval.rates.vec(pre.frac,lab.frac,cl=NULL))]
    data[,prod.rate:= get.production.rate(log(unlabeled.exon-unlabeled.intron),deg.rate,proc.rate)]
    ## data[,prod.rate:= get.production.rate.ss(unlabeled.exon,unlabeled.intron,labeled.exon, labeled.intron,deg.rate,proc.rate)]
    data[twosols==T,prod.rateb:= get.production.rate(log(unlabeled.exon),deg.rateb,proc.rateb)]
    data[twosols==T,c("prod.rate","deg.rate","proc.rate","prod.rateb","deg.rateb","proc.rateb"):=swap(prod.rate,deg.rate, proc.rate,prod.rateb,deg.rateb, proc.rateb)]
    return(data)
}


func <- function(logk,pre.frac,lab.frac){
    k  <- exp(logk)
    return(k/(k-1)*(log(pre.frac+k-1)-2*log(k)-log(pre.frac))-log(lab.frac-pre.frac)+log(pre.frac)+log(lab.frac*(k+1)-1))
}


dom.bound <- function(llab.frac){
    return(log(1-exp(llab.frac))-llab.frac)
}


c.from.k <- function(logk,pre.frac,lab.frac,tt){
    logc <- log(log(lab.frac-pre.frac)-log(pre.frac)-log(lab.frac*(exp(logk)+1)-1))-log(tt)
    return(logc)
}

solve.rates2 <- function(pre.frac,lab.frac,t=1){
    solve.rates(c(pre.frac,lab.frac),t)
}

solve.rates3 <- function(pre.frac,lab.frac,t=1){
    sol <- solve.rates(c(pre.frac,lab.frac),t)
    if (length(sol)==2){
        return(c(sol,NA_real_,NA_real_))
    }
    if (length(sol)==4){
        return(sol)
    }
}
##' Extracts processing and degradation rates from unlabled (pre-existing) and labeled intro to exon ratios for a single transcript. 
##'
##' See reference paper for the details of the method. Ratios expected to be smaller that one and labeled ratio is expected to be larger than the unlabeled ratio. 
##' @title 
##' @param frac a vector of size 2 containing :
##' 1. the prexisting intron to exon ratio
##' 2. the labeled intron to exon ratio
##' @param t the labeling time
##' @return a vector containing:
##' 1. the log degradation rate
##' 2. the log processing rate
##' the rates in units given by the inverse of the units of t
##' If their is only one solution the vector is of size 2. If there are two
##' solutions, the vector is of
##' size 4 (degr.rate1, proc.rate1, degr.rate2, proc.rate2).
##' If there are no solution the vector c(NA,NA) is returned. 
##' 
##'
##' @author Micha Hersch
solve.rates <- function(frac,t=1){
    if(is.na(sum(frac))){
        return(c(NA_real_,NA_real_))
    }
    pre.frac <- frac[1]
    lab.frac <- frac[2]
    eps <- 0.0000001
    lb <- max(dom.bound(log(lab.frac)),log(1-pre.frac))
    ub <- dom.bound(log(pre.frac))
    if(!is.finite(ub)){
        print(c(frac,lb,ub))
    }
    lbm  <- (1+sign(lb)*eps)*lb
    ubm  <- (1-sign(ub)*eps)*ub
    if(abs(lb)<eps){
        lbm  <- eps
    }
    if(abs(ub)<eps){
        ubm <- -eps
    }
    all <- uniroot.all(func,c(lbm,ubm),pre.frac=pre.frac,lab.frac = lab.frac)
    if(length(all)==1){
        logc <- c.from.k(all,pre.frac,lab.frac,t)
        return(c(logc-all,logc))
    }else{
        if(length(all)==2){
            logc1 <- c.from.k(all[1],pre.frac,lab.frac,t)
            logc2 <- c.from.k(all[2],pre.frac,lab.frac,t)
            return(c(logc1-all[1],logc1,logc2-all[2],logc2))
        }else{
            return(c(NA,NA))
        }
    }
}









r1 <- function(b,c,t=1){
    tryCatch({
        eps  <-  0.000000001
    thresh <- 5
    if(abs(b-c)<eps){
        return(-log(c*t+2))
    }
        if(b>c){
            lognum <- log(b-c)+log(b)-c*t
            logden  <- 2*log(b)-c*t + log(1-exp(2*(log(c)-log(b))+(c-b)*t))
        }else{
            lognum <- log(c-b)+log(b)-c*t
            logden <- 2*log(c)-b*t + log(1-exp(2*(log(b)-log(c))+(b-c)*t))
        }
        return(lognum-logden)
    },
    error=function(e){
        print(e)
        print(c(b,c))
        stop()
    })
}



t.from.c <- function(logc,log.pre.frac,log.lab.frac,nt=1){## nt is not used, here for modularity
    pf  <- exp(log.pre.frac)
    lf  <- exp(log.lab.frac)
    c  <- exp(logc)
    m1 <- rep(1,length(c))
    nat.t <- (-log(pf)-log(lf*(m1+c)-m1) +log(lf-pf))/c
    if(any(is.nan(nat.t))){
        print("nan")
        print(cbind(pf,lf,logc))
    }
    return(nat.t)
}

dt.dc <- function(logc,log.pre.frac,log.lab.frac){
    pf  <- exp(log.pre.frac)
    lf  <- exp(log.lab.frac)
    cc <- exp(logc)
    t1 <- t.from.c(logc,log.pre.frac,log.lab.frac)
    t2  <-  pf*lf/(cc*pf*(lf*(1+cc)-1))
    ret  <-  -(t1/cc+t2)
    return(ret)
}



r2 <- function(b,c,t=1){
    if(is.na(b) | is.na(c)){
        return(NA)
    }
    eps <- 0.000000001
    thresh <- 0
    if(abs(b-c)<eps){
        return(r2lim(c,t))
    }else{
        if(t<0){ ## never
            num <- b*(b-c)*(1-exp(-c*t))
            den <- b*b*(1-exp(-c*t))-c*c*(1-exp(-b*t))
            return(log(num/den))
        }else{
            if(t<eps){
                frac = sqrt(c/b)
            }else{
                frac <- sqrt(1-exp(-b*t))/sqrt(1-exp(-c*t))
            }
            lognum <- log(b)+log(abs(b-c))
            logden <-log(b+c*frac)+log(abs(b-c*frac))
            return(lognum-logden)
        }
    }
}


dr1.dc <- function(c,t){
   
    eps  <-  0.000000001
    if(abs(c-1)<eps){  # takes care of limit when c=1
        return(-t/(c*t+2)^2)
    }else{
        ec  <- exp(c*t)
        et  <- exp(t)
        num <- et+c*ec*(c-2+(c-1)*c*t)
        den <- (c-1)*(c*c*ec-et)
        ret  <- -num/den
        return(ret)
    }
}

dr1.dt <- function(c,t){
    eps  <-  0.000000001
    if(abs(c-1)<eps){  # takes care of limit when c=1
        return(c/(c*t+2)^2)
    }else{
        ec  <-  c*c*exp(c*t)
        num <- (c-1)*ec
        den  <- exp(t)-ec
        return(num/den)
    }
}


dr2.dc <- function(c,t){
    ec <- exp(c*t)
    ecm  <- ec-1
    et  <- exp(t)
    etm  <- et -1
    num  <- et + ec*ec*(et+(c-2)*c*etm)-ec*(2*et+c*etm*(c-2+(c-1)*c*t))
    den  <- (c-1)*ecm*(et+ec*(c*c*etm-et))
    return(-num/den)
    
}
dr2.dt <- function(c,t){
    ec <- exp(c*t)
    ecm  <- ec-1
    et  <- exp(t)
    etm  <- et -1
    num  <- c*c*ec*(ecm-c*etm)
    den  <- ecm*(et+ec*(c*c*etm-et))
    return(-num/den)
}

sq.error.line.nat <- function(logc,log.pre.frac,log.lab.frac,t){
    nat.t  <-  t.from.c(logc,central(log.pre.frac),central(log.lab.frac))
    if(is.na(nat.t)){
        nat.t <- t.from.c(logc,central(log.pre.frac),central(log.lab.frac))
    }
    cc <- exp(logc)
    bb <- 1
    lr1 <- r1(bb,cc,nat.t)
    lr2 <- r2(bb,cc,nat.t)
    ret <- (lr1-central(log.pre.frac))^2+(lr2-central(log.lab.frac))^2
    return (ret)# + sq.error.line(logc-log(nat.t),log.pre.frac,log.lab.frac,t))
}



ddc.sq.error.line.nat <- function(logc,log.pre.frac,log.lab.frac,t){
    nat.t  <-  t.from.c(logc,central(log.pre.frac),central(log.lab.frac))
    cc <- exp(logc)
    bb <- 1
    lr1 <- r1(bb,cc,nat.t)
    lr2 <- r2(bb,cc,nat.t)
    dtdc  <- dt.dc(logc,log.pre.frac,log.lab.frac)
    ret1 <- 2*(lr1-central(log.pre.frac))*(dr1.dc(cc,nat.t)+dr1.dt(cc,nat.t)*dtdc)
    ret2 <- 2*(lr2-central(log.lab.frac))*(dr2.dc(cc,nat.t)+dr2.dt(cc,nat.t)*dtdc)
    d.dc  <- ret1+ret2
    if(!is.finite(d.dc)){
        print("stop")
        stop("bad gradient")
    }
    return(d.dc*cc)
}

central <- function(x){
    return (mean(x,na.rm=T))
}

##' Performs function minimization on a bounded domain using a gradient method.
##'
##' Looks for two minimas and returns them both.
##' @title 
##' @param lpre.frac log unlabeled intron to exon ratio
##' @param llab.frac log labeled intron to exon ratio
##' @param init.par initial parameter value for searcu
##' @param lbm lower bound
##' @param ubm upper bound
##' @return Locations of two mimimas
##' @author Micha Hersch
moptim <- function(lpre.frac, llab.frac,init.par,lbm,ubm){
    ## try first minimum
    op <- tryCatch(optim(init.par,fn=sq.error.line.nat,gr=ddc.sq.error.line.nat,log.pre.frac=lpre.frac,## use gradient if possible
                         log.lab.frac=llab.frac,t=tp,method="Brent",lower=lbm,upper=ubm), 
                   error = function(e){optim(init.par,fn=sq.error.line.nat,log.pre.frac=lpre.frac, ## otherwise try without gradient
                                             log.lab.frac=llab.frac,t=1,method="Brent",lower=lbm,upper=ubm)})
    opt.par <- op$par
    marg  <- (ubm-lbm)/50

    ## look for second minimum below the first one
    op2 <- tryCatch(optim(lbm,fn=sq.error.line.nat,gr=ddc.sq.error.line.nat, log.pre.frac=lpre.frac,
                          log.lab.frac=llab.frac,t=tp,method="L-BFGS-B",
                          lower=lbm,upper=opt.par-marg),
                    error = function(e){print("trying again low ")
                        op2 <- optim(lbm+0.5*marg,fn=sq.error.line.nat,gr=NULL, log.pre.frac=lpre.frac,
                                     log.lab.frac=llab.frac,t=1,method="L-BFGS-B",
                                     lower=lbm,upper=opt.par-marg)
                    })
    ## look for second minimum above the first one.
    op3 <- tryCatch(optim(ubm-0.5*marg,fn=sq.error.line.nat,gr=ddc.sq.error.line.nat,log.pre.frac=lpre.frac,
                          log.lab.frac=llab.frac,t=1,method="L-BFGS-B",
                          lower=opt.par+marg,upper=ubm),
                    error = function(e){print("trying again up")
                        optim(lbm,fn=sq.error.line.nat,gr=NULL, log.pre.frac=lpre.frac,
                              log.lab.frac=llab.frac,t=tp,method="L-BFGS-B",
                              lower=ubm,lower=opt.par-marg)})
    
    if(op2$par + 1.1*marg < opt.par){
        opt.par <- c(op2$par,opt.par)
    }else if( op3$par - 1.1*marg > opt.par){ 
        opt.par <- c(opt.par,op3$par)
    }
    return(opt.par)
}


##' .. vectorized version of function optimization
##'
##' .. content for \details{} ..
##' @title 
##' @param pre.frac unlabled intron to exon ratios for multiple transcripts
##' @param lab.frac labled intron to exon ratios for the same transcripts
##' @param tp labeling time
##' @param cl cluster for parallelization (NULL if no parallelization)
##' @return a matrix where each row corresponds to a transcripts. The first column corresponds to the degradation rate, and the second column corresponds to the processing rate. A third and fourth column containing degradation and processing rates for a second mimima are provided if necessary.  
##' @author Micha Hersch
eval.rates.vec <- function(pre.frac,lab.frac,tp=1,cl=NULL){
    keep <-  which(!(pre.frac>1 | lab.frac>1 | is.na(pre.frac) | is.na(lab.frac) | pre.frac>=lab.frac))
    lpre.frac <- log(pre.frac[keep])
     llab.frac <- log(lab.frac[keep])
    eps <- 0.0000001
    marg.frac <- 50
    lb <- dom.bound(llab.frac)
    lb  <- pmax(lb, -10)
    ub <- dom.bound(lpre.frac)
    ub <- pmin(ub,10)
    lbm  <- lb+(ub-lb)/marg.frac
    ubm  <- ub-(ub-lb)/marg.frac
    init.par <- 0.5*(ub+lb)
    if(is.null(cl)){
        allcl  <- ifelse(lb < ub*(1-sign(lb)*eps)/(1+sign(ub)*eps),apply(cbind(lpre.frac,llab.frac,init.par,lbm,ubm),1,function(x)moptim(x[1],x[2], x[3],x[4],x[5])))
    }else{
        allcl  <- ifelse(lb < ub*(1-sign(lb)*eps)/(1+sign(ub)*eps),parApply(cl,cbind(lpre.frac,llab.frac,init.par,lbm,ubm),1,function(x)moptim(x[1],x[2], x[3],x[4],x[5])))
    }
    opt.par <- sapply(allcl,"[[",1)
    
    nat.t  <-log(t.from.c(opt.par,lpre.frac,llab.frac))
    ret  <- matrix(NA,nrow=length(pre.frac),ncol=4)
    ret[keep,1:2]  <- cbind(nat.t,opt.par+nat.t)
    ambig <- which(sapply(allcl,length)==2)
    if(length(ambig)>0){
        opt.par2 <- rep(NA,length(allcl))
        opt.par2[ambig] <- sapply(allcl[ambig],"[[",2)   
        nat.t2 <- rep(NA,length(allcl))
        nat.t2[ambig]  <-log(t.from.c(opt.par2[ambig],lpre.frac[ambig],llab.frac[ambig]))
        ret[keep,3:4]<- cbind(nat.t2,opt.par2+nat.t2)
    }
    return(ret)
}     
## lpre.m = log(pre.exon-pre.intron)
get.production.rate <- function(lpre.m,logb,logc,t=1){
    bb <- exp(logb)
    cc <- exp(logc)
    eb <- exp(-bb*t)
    ec <- exp(-cc*t)
    dd <- bb-cc
    si <- ifelse(dd>0,1,-1)## change sign
    logfac <- log(si*(ec - cc*eb/bb))- log(si*dd)
#    logfac <- log(si*(1 - cc/bb *exp(-dd*t))*ec)- log(si*dd)
    return(lpre.m -logfac)
}


get.norm.factor <- function(llab.exon,loga,logb,logc,t=1){
    aa <- exp(loga)
    bb <- exp(logb)
    cc <- exp(logc)
    eb <- exp(-bb*t)
    ec <- exp(-cc*t)
    dd <- bb-cc
    si <- ifelse(dd>0,1,-1) #change sign 
    logfac <- log(si*(dd/bb + cc/bb*eb-ec))+loga-log(si*dd)
    return(llab.exon-logfac)
}

get.production.rate.labeled <- function(llab.m, logb,logc,t=1){
    bb <- exp(logb)
    cc <- exp(logc)
    eb <- exp(-bb*t)
    ec <- exp(-cc*t)
    dd <- bb-cc
    si <- ifelse(dd>0,1,-1)## change sign
    loga <- llab.m + log(si*dd)-log(si*(cc*(eb-1)+bb-ec))
    return (loga)  
}

get.production.rate.ss <- function(pre.exon,pre.intron,lab.exon,lab.intron,logb,logc){
    bb <- exp(logb)
    cc <- exp(logc)
    dd <- cc+bb
    fac  <- -(dd*pre.intron-bb*pre.exon)/(dd*lab.intron-bb*lab.exon)
##    print(fac)
    loga = logb+log(pre.exon-pre.intron+fac*(lab.exon-lab.intron))
    return(loga)
}
