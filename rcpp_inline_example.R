src<-'
    Rcpp::NumericMatrix geno(genoR);
    int n_ind = geno.nrow();
    int n_mar = geno.ncol();
    Rcpp::NumericMatrix rec(n_mar, n_mar);
    double ct;
    for (int i = 0; i < (n_mar-1); i++) 
    {
        for (int j = (i+1); j < n_mar; j++) 
	{
            ct=0.0;
            for(int k=0; k < n_ind; k++)
	    {
                if(geno(k,i)!=geno(k,j))
                    ct++;
                    }
            rec(j,i)=rec(i,j)=ct/n_ind;
	}
    }
    return(rec);
    '

require(Rcpp)
require(inline)
require(fields)
require(onemap)

est_bc <- cxxfunction(signature(genoR = "numeric"), body = src, plugin="Rcpp")
dat<-as.matrix(read.table("mouse.txt"))
rec<-est_bc(dat)
image.plot(rec, col=rev(tim.colors()))

##More markers
source("simulate_diploid_populations.R")
dat.bc<-sim.pop.bc(n.ind = 250, n.mrk = 5000, ch.len = 200, missing = 0, n.ch = 1, verbose = FALSE)
dat.bc
dat.bc.t<-dat.bc$geno
system.time(rec<-est_bc(dat.bc.t))
choose(5000, 2)
image.plot(rec, col=rev(tim.colors()))
