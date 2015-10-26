##Kosambi mapping function
mf.k<-function (d) 
    0.5 * tanh(d/50)

##Inverse of Kosambi mapping function
imf.k<-function (r)
{
    r[r >= 0.5] <- 0.5 - 1e-14
    50 * atanh(2 * r)
}

##Simulates one diploid gamete
sim.gam<-function(n.mrk, r, init=0)
{
    g<-rep(init, n.mrk)
    for(i in 2:n.mrk)
        g[i]<-sample(c(g[i-1],0^g[i-1]), size=1, prob=c(1-r, r))
    g
}

##Simulates a backcross population
sim.pop.bc<-function(n.ind, n.mrk, ch.len, missing=0, n.ch=1, verbose=TRUE)
{
    n.mrk.per.ch<-ceiling(n.mrk/n.ch)
    mrk.names<-NULL
    for(j in 1:n.ch)
        mrk.names<-c(mrk.names, paste("ch", j, "M",1: n.mrk.per.ch   , sep=""))
    dat<-matrix(0,n.ind, n.mrk.per.ch)
    dat.final<-NULL
    r<-mf.k(ch.len/n.mrk)
    if(verbose) pb <- txtProgressBar(min = 1, max = n.ind, style = 3)
    for(j in 1:n.ch)
    {
        for(i in 1:n.ind)
        {
            if(verbose) setTxtProgressBar(pb, i)
            init<-sample(c(0,1), size=1)
            dat[i,]<-sim.gam(n.mrk.per.ch,r,init)
        }
        dat.final<-cbind(dat.final, dat)
    }
    if(verbose) close(pb)
    colnames(dat.final)<-mrk.names
    dat.final<-dat.final+1
    dat.final[sample(1:length(dat.final), size=length(dat.final)*missing/100)]<-0
    dat.mmk<-dat.final
    dat.mmk[dat.mmk==0]<-NA
    dat.mmk<-list(dat.mmk, type="bc")
    structure(list(geno = dat.final, geno.mmk = dat.mmk, n.ind = n.ind, n.mar = n.mrk,
                   segr.type = NA, segr.type.num=NA, phase=NA,
                   input="none.txt", n.phe=0, pheno = NA),  class = "bc.onemap")    
}

##Check type of markers in f2 populations
check.type<-function(x)
{
    if(sum(is.na(match(na.omit(unique(x)), c(1,2,3))))==0) return(1)
    else if (sum(is.na(match(na.omit(unique(x)), c(4,3))))==0) return(2)
    else if (sum(is.na(match(na.omit(unique(x)), c(5,1))))==0) return(3)
    else return(4)
}

##Simulates a F2 population 
sim.pop.f2<-function(n.ind, n.mrk, ch.len, dom43=0, dom51=0, missing=0, n.ch=1, verbose=TRUE)
{
    n.mrk.per.ch<-ceiling(n.mrk/n.ch)
    mrk.names<-NULL
    for(j in 1:n.ch)
        mrk.names<-c(mrk.names, paste("ch", j, "M",1: n.mrk.per.ch   , sep=""))
    dat<-matrix(0,n.ind, n.mrk.per.ch)
    dat.final<-NULL
    r<-mf.k(ch.len/n.mrk)
    if(verbose) pb <- txtProgressBar(min = 1, max = n.ind, style = 3)
    for(j in 1:n.ch)
    {
        for(i in 1:n.ind)
        {       
            if(verbose) setTxtProgressBar(pb, i)
            dat[i,]<-sim.gam(n.mrk.per.ch,r,sample(c(0,1), size=1))+sim.gam(n.mrk.per.ch,r,sample(c(0,1), size=1))
        }
        dat.final<-cbind(dat.final, dat)
    }
    if(verbose) close(pb)
    dat.final<-dat.final+1
    colnames(dat.final)<-mrk.names
    dom<-sample(1:ncol(dat.final), size=ncol(dat.final)*(dom43+dom51)/100)
    d1<-sample(dom, size=length(dom)*(dom43/(dom43+dom51)))
    d2<-dom[is.na(match(dom,d1))]
    for(i in d1)
        dat.final[dat.final[,i] < 3, i] <- 4
    for(i in d2)
        dat.final[dat.final[,i] > 1, i] <- 5
    dat.final[sample(1:length(dat.final), size=length(dat.final)*missing/100)]<-0
    dat.mmk<-dat.final
    dat.mmk[dat.mmk==0]<-NA
    type<-apply(dat.mmk, 2, check.type)
    dat.mmk<-list(dat.mmk, type="f2")
    code<-c(-2,-2,-3)
    structure(list(geno = dat.final, geno.mmk = dat.mmk, n.ind = n.ind, n.mar = n.mrk,
                   segr.type = NA, segr.type.num=type, phase=code[type],
                   input="none.txt", n.phe = 0, pheno = NA),  class = "f2.onemap")    
}

##Simulates an Outcross population
sim.pop.out<-function(n.ind, n.mrk, ch.len, missing=0, prob=c(1,1,1,1,1,1,1), n.ch=1, verbose=TRUE)
{
    n.mrk.per.ch<-ceiling(n.mrk/n.ch)
    mrk.names<-NULL
    for(j in 1:n.ch)
        mrk.names<-c(mrk.names, paste("ch", j, "M",1: n.mrk.per.ch   , sep=""))
    dat<-matrix(0,n.ind, n.mrk.per.ch)
    dat.final<-NULL
    type<-numeric(length(n.mrk.per.ch))
    type.final<-NULL
    r<-mf.k(ch.len/n.mrk)
    if(verbose) pb <- txtProgressBar(min = 1, max = n.ind, style = 3)
    for(j in 1:n.ch)
    {
        for(i in 1:n.ind)
        {
            if(verbose) setTxtProgressBar(pb, i)
            x1<-sim.gam(n.mrk.per.ch,r,sample(c(0,1), size=1))
            x2<-sim.gam(n.mrk.per.ch,r,sample(c(0,1), size=1))
            dat[i,] <- as.numeric(paste(x1,x2, sep=""))
        }
        for(i in 1:n.mrk.per.ch)
        {
            dt<-dat[,i]
            type[i]<-sample(1:7, 1, prob=prob)
            if(type[i]==1)
            {
                dt[dt==1]<-2
                dt[dt==0]<-1
                dt[dt==10]<-3
                dt[dt==11]<-4
            }
            else if(type[i]==2)
            {
                dt[dt==0]<-1
                dt[dt==10]<-2
                dt[dt==11]<-3
            }
            else if(type[i]==3)
            {
                dt[dt==1]<-2
                dt[dt==0]<-1
                dt[dt==10]<-1
                dt[dt==11]<-3
            }
            else if(type[i]==4)
            {
                dt[dt==1]<-2
                dt[dt==0]<-1
                dt[dt==10]<-2
                dt[dt==11]<-3
            }
            else if(type[i]==5)
            {
                dt[dt==0]<-1
                dt[dt==10]<-1
                dt[dt==11]<-2
            }
            else if(type[i]==6)
            {
                dt[dt==0]<-1
                dt[dt==10]<-2
                dt[dt==11]<-2
            }
            else if(type[i]==7)
            {
                dt[dt==1]<-2
                dt[dt==0]<-1
                dt[dt==10]<-1
                dt[dt==11]<-2
            }
            dat[,i]<-dt
        }
        dat.final<-cbind(dat.final, dat)
        type.final<-cbind(type.final, type)
    }
    if(verbose) close(pb)
    dat.final[sample(1:length(dat.final), size=length(dat.final)*missing/100)]<-0
    colnames(dat.final)<-mrk.names
    structure(list(geno = dat.final, n.ind = n.ind, n.mar = n.mrk, segr.type = as.character(type.final),
                   segr.type.num=type.final, input="none.txt", n.phe=0, pheno = NA), class = "outcross")
}
