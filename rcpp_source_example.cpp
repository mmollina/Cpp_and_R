#include <Rcpp.h>

// [[Rcpp::export]]
SEXP est_bc_source(Rcpp::NumericMatrix genoR) {
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
}
