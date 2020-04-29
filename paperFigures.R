source("SingleSampleRNAdynamics.R")
library("ggplot2")
library("data.table")


### genrating simulated data 

nsim <- 50000  # number of data points

tp <- t(exp(seq(-5,5))) # time points for trajectories

ntp  <- length(tp)

## generating random rates
ra <- exp(rnorm(nsim,0,2)) # production rate
rb <- exp(rnorm(nsim,0,2)) # degradation rate
rc <- exp(rnorm(nsim,0,2)) # processing rate

pre.p <- ra/rc*exp(-rc %*% tp) # pre-existing (or unlabeled) premature RNA
pre.m <- ra/(rb-rc)*exp(-rc %*% tp) - ra*rc/(rb*(rb-rc))*exp(-rb %*% tp) # pre-existing mature RNA

lab.p <- ra/rc*(1-exp(-rc %*% tp)) # labeled premature RNA
lab.m <- ra/rb*(1+rc/(rb-rc)*exp(-rb %*% tp))-ra/(rb-rc)*exp(-rc %*% tp)# labeled mature RNA

pre.frac <- pre.p/(pre.p+pre.m) ## unlabeled observables 
lab.frac <- lab.p/(lab.p+lab.m) ## labeled observables


ersd = 0.00 ## simulated noise level (not used because we also have real data)

pre.frac <- exp(log(pre.frac)+rnorm(nsim,sd=ersd))
lab.frac <- exp(log(lab.frac)+rnorm(nsim,sd=ersd))

pre.frac[pre.frac==0] <-NA
lab.frac[lab.frac==0] <- NA
sc <- which(tp==1) 
t <- tp[sc]
frac <-cbind(pre.frac[,sc],lab.frac[,sc])

## estimating rates from observables
my.ratesn  <-  apply(frac,1,solve.rates,t=t)

#### figure 2 accuracy of solutions  ################
show.estimates <- function(lrb,lrc,rates){
  #  par(mfcol=c(1,2))
    plot(lrb,rates[1,],pch=20,col=rgb(0,0,1,0.05),ylim=c(-6,6),xlim=c(-6,6),xlab="true degradation rate", ylab="estimated degradation rate")
    abline(a=0,b=1)
    plot(lrc,rates[2,],pch=20,col=rgb(1,0,0,0.05),ylim=c(-6,6),xlim=c(-6,6),xlab="true processing rate", ylab="estimated processing rate")
    abline(a=0,b=1)
}


par(mfcol=c(2,1))
show.estimates(log(rb),log(rc),my.ratesn)


#### figure 3 number of solutions  ######################

get.table <- function(lrb,lrc,rates,frac){
    diffb = sqrt((log(rb)-rates[1,])^2)
    diffc = sqrt((log(rc)-rates[2,])^2)
    res.data <- data.table(rb=rb,rc=rc,lrb=log(rb), lrc=log(rc),diffb=diffb,diffc=diffc,preex.ratio=frac[,1],label.ratio=frac[,2])
    res.data[,col:=ifelse(is.na(diffc),"A",ifelse(diffc>0.1,"B","C"))]
    res.data[,k:= lrc-lrb]
    return(res.data)
}
resn  <- get.table(log(rb),log(rc),my.ratesn,frac)

ggplot(data=resn,aes(x=preex.ratio,y=label.ratio,color=col)) + geom_point(aes(stroke=0))+geom_line(data=data.frame(x=seq(0,1,0.05),col="green"),aes(x=x,y=1/(2-x),colour=as.factor(col)),show.legend=FALSE,guide=FALSE) + theme_classic() + scale_color_manual(name ="rate estimate", labels=c("ambiguous","wrong","correct","b=1/(2-a)"), values=c(rgb(0,0,1,0.2),rgb(1,0,0,0.2),rgb(0,0,0,0.2),"green"))+labs(x="unlabeled ratio (a)",y = "labeled ratio (b)") + scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1))


### figure 4 trajectories

trajs <- data.table(lra=log(ra),lrb=log(rb),lrc=log(rc),k=log(rc)-log(rb),preex.ratio=pre.frac,label.ratio=lab.frac)
p1  <- ggplot(trajs,aes(x=preex.ratio.V2,y=label.ratio.V2,z=k))+theme_classic()
for (i in seq(-6,6)){
    for (j in seq_along(tp)){
        cx  <- paste0("preex.ratio.V",j)
        cy  <- paste0("label.ratio.V",j)
        topl <- trajs[abs(2*k-i)<0.01]
        topl$col <- i
        p1  <- p1+ geom_point(topl,mapping=aes_string(x=cx,y=cy,col="col"),shape=20,size=1)
    }
}
p1  <- p1 + scale_color_gradient2(low="blue",high="red",mid="black",midpoint=0,name= "k [log]")
p1 <- p1 +labs(x="unlabeled ratio (a)",y = "labeled ratio (b)")
p2  <-  p1 +geom_line(data=data.frame(x=seq(0,1,0.05),col=1.0,k=100),aes(x=x,y=1/(2-x)),color="green",show.legend=FALSE)
kk  <- seq(-4,4,2)
p2  <- p2 + annotate("text",x = 1/(1+exp(kk/2)),y =1.05,label=paste("log(k)=",kk/2))
p2 
