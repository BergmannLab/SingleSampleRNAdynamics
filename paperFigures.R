source("SingleSampleRNAdynamics.R")
library("ggplot2")
library("data.table")
library("EnsDb.Mmusculus.v79")
library("org.Mm.eg.db")

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
sc <- 5 # using a single time point to infer rates.
t <- tp[sc]
frac <-cbind(pre.frac[,sc],lab.frac[,sc])

## estimating rates from observables
my.ratesn  <-  apply(frac,1,solve.rates,t=t)

#### figure 2 accuracy of solutions  ################
show.estimates <- function(lrb,lrc,rates){
  #  par(mfcol=c(1,2))
    plot(lrb,rates[,1],pch=20,col=rgb(0,0,1,0.05),ylim=c(-6,6),xlim=c(-6,6),xlab="true degradation rate", ylab="estimated degradation rate")
    abline(a=0,b=1)
    plot(lrc,rates[,2],pch=20,col=rgb(1,0,0,0.05),ylim=c(-6,6),xlim=c(-6,6),xlab="true processing rate", ylab="estimated processing rate")
    abline(a=0,b=1)
}

to.show <- sapply(my.ratesn,function(x)length(x)==2)
rates <- cbind(sapply(my.ratesn,"[[",1),sapply(my.ratesn,"[[",2))
par(mfcol=c(2,1))
show.estimates(log(rb[to.show]),log(rc[to.show]),rates[to.show,])
dev.copy2pdf(file="simul_rates.pdf",onefile=T)

#### figure 3 number of solutions  ######################

get.table <- function(lrb,lrc,rates,ambig,frac){
    diffb = sqrt((log(rb)-rates[,1])^2)
    diffc = sqrt((log(rc)-rates[,2])^2)
    res.data <- data.table(rb=rb,rc=rc,lrb=log(rb), lrc=log(rc),diffb=diffb,diffc=diffc,preex.ratio=frac[,1],label.ratio=frac[,2])
    res.data[,col:=ifelse(!ambig,"A",ifelse(diffc>0.1,"B","C"))]
    res.data[,k:= lrc-lrb]
    return(res.data)
}
resn  <- get.table(log(rb),log(rc),rates,to.show,frac)

ggplot(data=resn,aes(x=preex.ratio,y=label.ratio,color=col)) + geom_point(aes(stroke=0))+geom_line(data=data.frame(x=seq(0,1,0.05),col="green"),aes(x=x,y=1/(2-x),colour=as.factor(col)),size=1.5,show.legend=FALSE,guide=FALSE) + theme_classic() + scale_color_manual(name ="rate estimate", labels=c("ambiguous","wrong","correct","b=1/(2-a)"), values=c(rgb(0,0,1,0.2),rgb(1,0,0,0.2),rgb(0.5,0.5,0.5,0.2),"green"))+labs(x="unlabeled ratio (a)",y = "labeled ratio (b)") + scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) + theme(text = element_text(size = 20)) 

dev.copy2pdf(file="simul_phase.pdf",onefile=T)


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
dev.copy2pdf(file="trajectories.pdf",onefile=T)

## loading  and preparing data ####

data.WT10 <- fread("./data/transcripts_tpm.csv")

## Figure 5. plotting the raw data
p1 <- ggplot(data.WT10,aes(x=log(P_1_WT_10_8d.intron/P_1_WT_10_8d.exon),y = log(L_1_WT_10_8d.intron/L_1_WT_10_8d.exon)))
##p1 <- p1+ geom_point(alpha=0.01,pch=20)+scale_x_continuous(limits = c(-7, 3)) + scale_y_continuous(limits = c(-7, 3))
p1 <- p1 + geom_rect(xmin=-10,xmax=-7,ymin=0,ymax=10,fill="#CCFFFF",alpha=0.1)
p1 <- p1 + geom_rect(xmin=-10,xmax=-7,ymin=-10,ymax=0,fill="#FFCCCC",alpha=0.1)
p1 <- p1+ geom_point(aes(alpha=log(1+P_1_WT_10_8d.exon + P_1_WT_10_8d.intron)),pch=20)+scale_x_continuous(limits = c(-7, 3)) + scale_y_continuous(limits = c(-7, 3)) + scale_alpha(range=c(0,0.1),guide=F)

p1 <- p1 + geom_line(data=data.frame(x=seq(0,1,0.05)),aes(x=log(x),y=log(1/(2-x))),colour="green") + geom_abline(slope=1,colour=rgb(1,0,0,0.5)) + geom_abline(slope=0,colour=rgb(0,0,1,0.5))+theme_classic()+theme(text = element_text(size = 20))+ ylab("labeled ratio (b) [log]")+ xlab("unlabeled ratio (a) [log]") + xlab(bquote("observed unlabeled ratio "*r[u]~ "[log]")) + ylab(bquote("observed labeled ratio "*r[l]~ "[log]"))
p2 <- ggplot(data.WT10[biotype=="protein_coding" & is.finite(lab.frac.1)],aes(x=lab.frac.1>1,y=log(P_1_WT_10_8d.exon),fill=lab.frac.1>1))+geom_boxplot(notch=T)+theme_classic()
p2 <- p2 + xlab("")+ylab("expression [log TPM]")+scale_x_discrete(labels=c(expression(r[l] <= 1),expression(r[l] > 1))) +theme(legend.position = "none",text = element_text(size = 16))+scale_fill_manual(values=c("#FFCCCC", "#CCFFFF"))
p1+ annotation_custom(ggplotGrob(p2),xmin=-1.5,xmax=3,ymin=-8,ymax=-3)
dev.copy2pdf(file="data_phase.pdf",onefile=T)


## getting transcript annotations
edb  <- EnsDb.Mmusculus.v79
biot <- select(edb, keys=data.WT10$rn, columns=c("TXID", "TXBIOTYPE"),keytype="TXID")
idx <- match(biot$TXID,data.WT10$rn)
data.WT10$biotype <- NA
data.WT10$biotype[idx] <- biot$TXBIOTYPE

rownames(data.WT10) <- data.WT10$rn
data.WT10[,rn:=NULL]

## getting gene symbols


gene.info <- select(edb, keys=rownames(data.WT10), columns=c("TXID", "TXBIOTYPE", "SYMBOL"),keytype="TXID")
data.WT10$txid = rownames(data.WT10)
data.WT10[,symbol:=gene.info$SYMBOL[match(txid,gene.info$TXID)]]


nrep  <- 3

pre.frac <- data.WT10[,seq(nrep+1,2*nrep),with=F]/data.WT10[,seq(1:nrep),with=F]
lab.frac  <- data.WT10[,seq(3*nrep+1,4*nrep),with=F]/data.WT10[,seq(2*nrep+1,3*nrep),with=F]
data.WT10[,paste("pre.frac",seq(nrep),sep=".") := pre.frac]
data.WT10[,paste("lab.frac",seq(nrep),sep=".") := lab.frac]

## classifying transcripts:
## C: solvable with a single solution
## A: potentially ambiguous solutions
## B: bad and discarded
## L: limit cases, k may be infered but not the rates (rates are too fast)

data.WT10[,solvable1:=ifelse(lab.frac.1>1/(2-pre.frac.1),"C","A")]
data.WT10[,solvable2:=ifelse(lab.frac.2>1/(2-pre.frac.2),"C","A")]
data.WT10[,solvable3:=ifelse(lab.frac.3>1/(2-pre.frac.3),"C","A")]


EPS = exp(-30)
data.WT10[pre.frac.1 < EPS,solvable1:="L"]
data.WT10[pre.frac.2 < EPS,solvable2:="L"]
data.WT10[pre.frac.3 < EPS,solvable3:="L"]

data.WT10[lab.frac.1>1 | pre.frac.1 >= lab.frac.1,solvable1:="B"]
data.WT10[lab.frac.2>1 | pre.frac.2 >= lab.frac.2,solvable2:="B"]
data.WT10[lab.frac.3>1 | pre.frac.3 >= lab.frac.3,solvable3:="B"]


pl <- length(setdiff(data.WT10[biotype=="protein_coding" &  P_1_WT_10_8d.exon>1,symbol],data.WT10[biotype=="protein_coding" &  P_1_WT_10_8d.exon>1 & solvable1 %in% c("A","C","L"),symbol]))/length(unique(data.WT10[biotype=="protein_coding" &  P_1_WT_10_8d.exon>1,symbol]))
print(paste("proportion of genes with bad ratio",pl))

table(data.WT10$solvable1)
table(data.WT10$solvable2)
table(data.WT10$solvable3)


tt <- table(data.WT10[biotype=="protein_coding"& L_1_WT_10_8d.exon>10 & L_1_WT_10_8d.intron>10,solvable1])
(tt[2]+tt[4])/sum(tt)

tt <- table(data.WT10[ L_1_WT_10_8d.exon>10 & L_1_WT_10_8d.intron>10,solvable1])
(tt[2]+tt[4])/sum(tt)

cl <- NULL ## change this to cl  <- init.cluster(10) if you want to work in parallel

## computing rates ####

refactor <- function(mat){return(list(mat[1,],mat[2,]))}
refactor2 <- function(mat){return(list(mat[,1],mat[,2]))}
refactor4 <- function(mat){return(list(mat[1,],mat[2,],mat[3,],mat[4,]))}
refactor3 <- function(mat){return(list(mat[1,],mat[2,],mat[3,]))}

## computing processing and degradation rates for solvable cases
data.WT10[solvable1=="C",c("deg.rate1","proc.rate1"):=refactor(mapply(solve.rates2,pre.frac.1,lab.frac.1))]
data.WT10[solvable2=="C",c("deg.rate2","proc.rate2"):=refactor(mapply(solve.rates2,pre.frac.2,lab.frac.2))]
data.WT10[solvable3=="C",c("deg.rate3","proc.rate3"):=refactor(mapply(solve.rates2,pre.frac.3,lab.frac.3))]

## computing processing and rates for potentially ambiguous cases
data.WT10[solvable1=="A",c("deg.rate1","proc.rate1","deg.rate1b","proc.rate1b"):=refactor4(mapply(solve.rates3,pre.frac.1,lab.frac.1))]
data.WT10[solvable2=="A",c("deg.rate2","proc.rate2","deg.rate2b","proc.rate2b"):=refactor4(mapply(solve.rates3,pre.frac.2,lab.frac.2))]
data.WT10[solvable3=="A",c("deg.rate3","proc.rate3","deg.rate3b","proc.rate3b"):=refactor4(mapply(solve.rates3,pre.frac.3,lab.frac.3))]

## finding out ambiguous cases
data.WT10[,twosols.1:= !is.na(deg.rate1b)]
data.WT10[,twosols.2:= !is.na(deg.rate2b)]
data.WT10[,twosols.3:= !is.na(deg.rate3b)]

## optimization procedure for unsolvable cases
data.WT10[solvable1=="A" & !twosols.1, c("deg.rate1","proc.rate1"):=refactor2(eval.rates.vec(pre.frac.1,lab.frac.1,cl=cl))]
data.WT10[solvable2=="A" & !twosols.2, c("deg.rate2","proc.rate2"):=refactor2(eval.rates.vec(pre.frac.2,lab.frac.2,cl=cl))]
data.WT10[solvable3=="A" & !twosols.3, c("deg.rate3","proc.rate3"):=refactor2(eval.rates.vec(pre.frac.3,lab.frac.3,cl=cl))]


## computing production rates
data.WT10[,prod.rate1:= get.production.rate(log(P_1_WT_10_8d.exon),deg.rate1,proc.rate1)]
data.WT10[,prod.rate2:= get.production.rate(log(P_2_WT_10_8d.exon),deg.rate2,proc.rate2)]
data.WT10[,prod.rate3:= get.production.rate(log(P_3_WT_10_8d.exon),deg.rate3,proc.rate3)]

## also for ambiguous cases
data.WT10[twosols.1==T,prod.rate1b:= get.production.rate(log(P_1_WT_10_8d.exon),deg.rate1b,proc.rate1b)]
data.WT10[twosols.2==T,prod.rate2b:= get.production.rate(log(P_2_WT_10_8d.exon),deg.rate2b,proc.rate2b)]
data.WT10[twosols.3==T,prod.rate3b:= get.production.rate(log(P_3_WT_10_8d.exon),deg.rate3b,proc.rate3b)]

## a simple heuristics to guess which of the two solutions is correct, take the one with the smaller the production rate
swap.ambiguous <- TRUE
if(swap.ambiguous){
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
    data.WT10[twosols.1==T,c("prod.rate1","deg.rate1","proc.rate1","prod.rate1b","deg.rate1b","proc.rate1b"):=swap(prod.rate1,deg.rate1, proc.rate1,prod.rate1b,deg.rate1b, proc.rate1b)]
    data.WT10[twosols.2==T,c("prod.rate2","deg.rate2","proc.rate2","prod.rate2b","deg.rate2b","proc.rate2b"):=swap(prod.rate2,deg.rate2, proc.rate2,prod.rate2b,deg.rate2b, proc.rate2b)]
    data.WT10[twosols.3==T,c("prod.rate3","deg.rate3","proc.rate3","prod.rate3b","deg.rate3b","proc.rate3b"):=swap(prod.rate3,deg.rate3, proc.rate3,prod.rate3b,deg.rate3b, proc.rate3b)]
}

k.name <- paste("k",seq(3),sep=".")
k.name  <- c(k.name,paste0(k.name,"b"))
data.WT10[, c("k.1","k.2","k.3","k.1b","k.2b","k.3b") := list(proc.rate1-deg.rate1,proc.rate2-deg.rate2,proc.rate3-deg.rate3,proc.rate1b-deg.rate1b,proc.rate2b-deg.rate2b,proc.rate3b-deg.rate3b)]

## gathering some stats to compare across replicates
data.WT10[, sd.k:= apply(cbind(proc.rate1-deg.rate1,proc.rate2-deg.rate2,proc.rate3-deg.rate3),1,sd)]
data.WT10[, sd.proc:= apply(cbind(proc.rate1,proc.rate2,proc.rate3),1,sd)]
data.WT10[, sd.deg:= apply(cbind(deg.rate1,deg.rate2,deg.rate3),1,sd)]
data.WT10[, sd.prod:= apply(cbind(prod.rate1,prod.rate2,prod.rate3),1,sd)]
data.WT10[, avg.deg:= apply(cbind(deg.rate1,deg.rate2,deg.rate3),1,mean)]
data.WT10[, avg.k:= apply(cbind(proc.rate1-deg.rate1,proc.rate2-deg.rate2,proc.rate3-deg.rate3),1,mean)]
data.WT10[, avg.proc:= apply(cbind(proc.rate1,proc.rate2,proc.rate3),1,mean)]
data.WT10[, avg.prod:= apply(cbind(prod.rate1,prod.rate2,prod.rate3),1,mean)]

### preparing Fig. 6


## projecting rates on the abundance and reactivity axes for 
data.WT10[biotype=="protein_coding" &  !is.na(deg.rate1),concentration:= (prod.rate1-deg.rate1)/2]
data.WT10[biotype=="protein_coding" & !is.na(deg.rate1),reactivity:= (prod.rate1+deg.rate1)/2]

print("correlations:")
data.WT10[biotype=="protein_coding" & !is.na(deg.rate1),cor(concentration,reactivity)]
data.WT10[biotype=="protein_coding" &  !is.na(deg.rate1),cor(deg.rate1,prod.rate1)]


go.cat <- c("transcription", "monosaccharide metabolism")
go.fnames <- list()
go.fnames[["transcription"]] <- "transcription_go.csv"
go.fnames[["monosaccharide metabolism"]] <- "monosaccharide_met_go.csv"

go.dir <- "./data/"
evid <- c("EXP","IMP","IDA")

cat.pos <- data.table(cat=go.cat)
data.WT10[,go_func:=as.character(NA)]

for (gcat in go.cat){
#    data.WT10[,go_func:=as.character(NA)]
    gene.list <- fread(paste0(go.dir,go.fnames[[gcat]]))
    keep <- which(gene.list[["GO EVIDENCE CODE"]] %in% evid)
    data.WT10[symbol %in% gene.list$SYMBOL[keep],go_func:=ifelse(is.na(go_func),gcat,"pleiotropic")]

    cat.pos[cat==gcat,reactivity:= data.WT10[biotype=="protein_coding" & go_func==gcat,median(reactivity,na.rm=T)]]
    ss <- wilcox.test(data.WT10[biotype=="protein_coding"  & go_func==gcat,reactivity],data.WT10[biotype=="protein_coding" & is.na(go_func),reactivity])
    cat.pos[cat==gcat,reactivity.pval:=ss$p.value]
    cat.pos[cat==gcat,concentration:= data.WT10[biotype=="protein_coding" & go_func==gcat,median(concentration,na.rm=T)]]
    ss <- wilcox.test(data.WT10[biotype=="protein_coding" & go_func==gcat,concentration],data.WT10[biotype=="protein_coding"  & is.na(go_func),concentration])
    cat.pos[cat==gcat,concentration.pval:=ss$p.value]
    cat.pos[cat==gcat,n:=data.WT10[biotype=="protein_coding" & go_func==gcat,.N]]
}
data.WT10[go_func=="pleiotropic",go_func:=NA]

## plotting Fig. 6
ax2 <- data.table(from.x=c(0,0),from.y=c(0,0),to.x=c(3,4),to.y=c(3,-4),y.off=c(1,-1),label=c("responsiveness","abundance"))
p1 <- ggplot(data.WT10[biotype=="protein_coding"],aes(x=prod.rate1,y=deg.rate1)) + scale_x_continuous(limits = c(-3, 10)) + scale_y_continuous(limits = c(-5.5, 4.5))
p1 <- p1 + geom_abline(slope=1,intercept=seq(-8,0,2),alpha=0.1,linetype="solid",color="brown")+geom_abline(slope=-1,intercept=seq(-4,6,2),alpha=0.1,linetype="dashed",color="brown")
p1 <- p1 + geom_point(alpha = 0.02,color="black")
p1 <- p1 + geom_point(data=data.WT10[biotype=="protein_coding"& !is.na(go_func)],aes(color=go_func),pch=19,alpha = 1,cex=1)
p1 <- p1 + labs(x="synthesis rate [log]", y="degradation rate[log]",color="GO categories")
p1 <- p1 + theme_classic() +theme(legend.position=c(0.8,1),text = element_text(size=20),legend.text = element_text(size = 12),legend.title = element_text(size = 16))+coord_fixed(ratio = 1)
p1 <- p1 + geom_segment(aes(x=from.x,y=from.y,xend=to.x,yend=to.y),color="black",data=ax2,arrow = arrow(length = unit(0.2, "cm")),lineend = "round") + geom_text(aes(x=to.x,y=to.y+0.5*sign(to.y),label=label),data=ax2,size=6)
p1
med = data.WT10[biotype=="protein_coding",list(concentration=mean(concentration,na.rm=T),reactivity=mean(reactivity,na.rm=T),prod.rate=mean(prod.rate1,na.rm=T),deg.rate=mean(deg.rate1,na.rm=T) ) ,by=go_func]
p1 <- p1 + geom_point(data=med,aes(x=concentration,y=-concentration,color=go_func),pch=15,cex=3)
p1 <- p1 + geom_point(data=med,aes(x=reactivity,y=reactivity,color=go_func),pch=15,cex=3)
p1
dev.copy2pdf(file="real_rates.pdf",onefile=T)

p2 <- ggplot(data.WT10[biotype=="protein_coding"],aes(x=concentration,y=reactivity)) + scale_x_continuous(limits = c(-0.5, 4.5)) + scale_y_continuous(limits = c(-5, 5))
p2 <- p2 + geom_point(data=data.WT10[biotype=="protein_coding"& !is.na(go_func)],aes(color=go_func),pch=19,alpha = 0.5,cex=1)
p2 <- p2 + labs(x="steady-state abundance [log]", y="responsiveness [log]")
p2 <- p2 + theme_classic() +theme(legend.position = "none",text = element_text(size=20),legend.text = element_text(size = 12),legend.title = element_text(size = 16))+coord_fixed(ratio = 1,clip="off",expand=F)
p2 <- p2 + geom_point(data=med,aes(x=concentration,y=reactivity,color=go_func),pch=17,cex=3)
##p2 <- p2 + geom_segment(data=med,aes(xend=-1,yend=reactivity,color=go_func))
p2 <- p2 + geom_hline(aes(yintercept=reactivity,color=go_func),data=med,linetype="dotted")+ geom_vline(aes(xintercept=concentration,color=go_func),data=med,linetype="dotted")
p2 <- p2 + geom_point(data=med,aes(y=reactivity,x=-Inf,color=go_func),pch=15,cex=3)+geom_point(data=med,aes(y=-Inf,x=concentration,color=go_func),pch=15,cex=3)
p2
dev.copy2pdf(file="real_rates_rot.pdf",onefile=T)



wilcox.test(data.WT10[biotype=="protein_coding" & go_func=="monosaccharide metabolism",reactivity],data.WT10[biotype=="protein_coding" & is.na(go_func),reactivity])
wilcox.test(data.WT10[biotype=="protein_coding" & go_func=="monosaccharide metabolism",concentration],data.WT10[biotype=="protein_coding" & is.na(go_func),concentration])
wilcox.test(data.WT10[biotype=="protein_coding" & go_func=="transcription",reactivity],data.WT10[biotype=="protein_coding" & is.na(go_func),reactivity])
wilcox.test(data.WT10[biotype=="protein_coding" & go_func=="transcription",concentration],data.WT10[biotype=="protein_coding" & is.na(go_func),concentration])


p1 <- ggplot(data.WT10[biotype=="protein_coding"],aes(x=prod.rate1,y=proc.rate1))+ scale_x_continuous(limits = c(-2, 8)) + scale_y_continuous(limits = c(-4, 1.5))
p1 <- p1 +  geom_point(pch=19,alpha = 0.05,cex=1) +scale_color_gradient2(low="black",mid="blue",high="red",midpoint=5,na.value = "grey")
p1 <- p1 + theme_classic() +theme(text = element_text(size=20))+coord_fixed(ratio = 1)
p1 <- p1 + labs(x="synthesis rate [log]", y="processing rate [log]",color="intron size [log]")
p1
dev.copy2pdf(file="real_proc.pdf",onefile=T)

## p1 <- ggplot(data.WT10[biotype=="protein_coding"],aes(x=reactivity,y=proc.rate1))+ scale_x_continuous(limits = c(-4, 5)) + scale_y_continuous(limits = c(-4, 1.5))
## p1 <- p1 +  geom_point(pch=19,alpha = 0.1,cex=1) 
## p1 <- p1 + theme_classic() +theme(text = element_text(size=20))+coord_fixed(ratio = 1)
## p1 <- p1 + labs(x="responsiveness [log]", y="processing rate [log]")
## p1

dev.copy2pdf(file="real_proc.pdf",onefile=T)
data.WT10[biotype=="protein_coding",cor(prod.rate1,proc.rate1,use="p")]

## comparing with Herzog et al.

data.dir <- "./data/"
am.fname <- "rates_herzog.txt"
adata <- fread(paste(data.dir,am.fname,sep="/"),fill=T)
colnames(adata)[which(colnames(adata)=="Half-life (h)")] <- "Half.life"
res <- select(org.Mm.eg.db, keys=adata$Name , columns=c("SYMBOL","ENSEMBLTRANS"), keytype="SYMBOL")
ma.idx <- match(rownames(data.WT10),res$ENSEMBLTRANS)
data.WT10[,gene.name:= res$SYMBOL[ma.idx]]
data.WT10b = merge(data.WT10[biotype=="protein_coding"],adata,by.x="gene.name",by.y="Name",all.y=T)
data.WT10b[,sum1:=P_1_WT_10_8d.exon+L_1_WT_10_8d.exon]
data.WT10b[,sum2:=P_2_WT_10_8d.exon+L_2_WT_10_8d.exon]
data.WT10b[,sum3:=P_3_WT_10_8d.exon+L_3_WT_10_8d.exon]


aa1 <- data.WT10b[solvable1 %in% c("A","C"),list(fin.deg=sum(deg.rate1*sum1,na.rm=T)/sum(sum1,na.rm=T),slam.deg=mean(log(log(2))-log(Half.life),na.rm=T),msum=sum(sum1,na.rm=T),solvable=solvable1[which.max(sum1)],frac.solv=sum(sum1[solvable1=="A"])/sum(sum1),lab.frac=sum(lab.frac.1*sum1)/sum(sum1),pre.frac=sum(pre.frac.1*sum1)/sum(sum1),k=sum(k.1*sum1)/sum(sum1)),by=gene.name]

aa2 <- data.WT10b[solvable2 %in% c("A","C"),list(fin.deg=sum(deg.rate2*sum2,na.rm=T)/sum(sum2,na.rm=T),slam.deg=mean(log(log(2))-log(Half.life),na.rm=T),msum=sum(sum2,na.rm=T),solvable=solvable2[which.max(sum2)],frac.solv=sum(sum2[solvable2=="A"])/sum(sum2),lab.frac=sum(lab.frac.2*sum2)/sum(sum2),pre.frac=sum(pre.frac.2*sum2)/sum(sum2),k=sum(k.2*sum2)/sum(sum2)),by=gene.name]

aa3 <- data.WT10b[solvable3 %in% c("A","C"),list(fin.deg=sum(deg.rate3*sum3,na.rm=T)/sum(sum3,na.rm=T),slam.deg=mean(log(log(2))-log(Half.life),na.rm=T),msum=sum(sum3,na.rm=T),solvable=solvable3[which.max(sum3)],frac.solv=sum(sum3[solvable3=="A"])/sum(sum3),lab.frac=sum(lab.frac.3*sum3)/sum(sum3),pre.frac=sum(pre.frac.3*sum3)/sum(sum3),k=sum(k.3*sum3)/sum(sum3)),by=gene.name]

data.WT10b[,asum:=(sum1+sum2+sum3)/3]
data.WT10[,rn := rownames(data.WT10)]

# Fig 7

thresh <- 200
ss <- summary(lm(fin.deg+log(6) ~ slam.deg,data=aa1[msum>thresh],weights=1-lab.frac))

p1 <- ggplot(aa1[msum>thresh],aes(y=fin.deg+log(6),x=slam.deg))+ geom_point(aes(alpha=1-lab.frac))
p1 + geom_abline(slope=ss$coefficient[2,1],intercept=ss$coefficient[1,1],colour="red")+ theme_classic()+annotate("text",x=-0.75,y=-4,label=paste("R =",format(100*sqrt(ss$r.squared),digits=0),"%"),size=8)+ labs(x="slam-seq estimate", y="single sample estimate ", title="degradation rate [log]")+ guides(alpha=FALSE)+theme(text = element_text(size=20))

dev.copy2pdf(file="valid_t100.pdf",onefile=T)

#### Figure 8 ##########################333

cor.func <- function(th,aa){
    ss <- summary(lm(fin.deg+log(6) ~ slam.deg,data=aa[msum>th & is.finite(fin.deg)],weights=1-lab.frac))
    return(c(sqrt(ss$r.squared),ss$df[2]))
}
cor.data <- data.table(thresh=c(0,100,200,500,1000))
cort <- matrix(NA,5,3)
dft  <- matrix(NA,5,3)
##laa  <- list(aa1,aa2,aa3,aa4,aa5)
laa  <- list(aa1,aa2,aa3)
for (j in seq(length(laa))){
    var =paste0(c("cor","df"),j)
    cor.data[,(var) := refactor(sapply(seq(5),function(i)cor.func(cor.data$thresh[i],laa[[j]])))]
}
cor.data[,quant1:=df1/df1[1]]
cor.data[,quant2:=df2/df2[1]]
cor.data[,quant3:=df3/df3[1]]

p1 <- ggplot(cor.data,aes(x=thresh/2))+geom_line(aes(y=cor1,color="replicate 1"))+geom_line(aes(y=cor2,color="replicate 2"))+geom_line(aes(y=cor3,color="replicate 3"))+scale_y_continuous(lim=c(0,0.7))+theme_classic()+theme(text = element_text(size=20),legend.position=c(0.8,0.5),axis.text=element_text(angle=0,size=12)) +geom_point(data=cor.data[thresh==200,],mapping=aes(y=cor1, color="replicate 1"),size=2)
p1+ scale_x_continuous(breaks=cor.data[,thresh/2],label=cor.data[,paste0(thresh/2,"\n(",format(100-100*(quant1+quant2+quant3)/3,digits=2),"%)")],name="expression threshold [TPM (and quantiles)]")+ labs(x="expression threshold [TPM (and quantiles)]",y="correlation",color="")


dev.copy2pdf(file="correlations.pdf",onefile=T)
