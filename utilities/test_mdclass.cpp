// testing multidim arrays in Rcpp

#include <Rcpp.h>
#include <R.h>
#include "mdarrclass.h"

using namespace Rcpp;

//[[Rcpp::export]]
RcppExport SEXP mdarr(SEXP dms_,SEXP vin_) 
{
  
  IntegerVector dm = as<Rcpp::IntegerVector>(dms_);
  NumericVector vin = as<Rcpp::NumericVector>(vin_);
  int dd1 = dm[0];
  int dd2 = dm[1];
  int dd3 = dm[2];
  int dd4 = dm[3];
  NumericVector resv(dd1*dd2*dd3*dd4);

  arr4d res;
  res.arr4(dd1,dd2,dd3,dd4);

  int i,j,k,l,elem;
  for(l=0;l<dd4;l++) {
    for(k=0;k<dd3;k++) {
      for(j=0;j<dd2;j++) { 
        for(i=0;i<dd1;i++) {
    
          elem = dd3*dd2*dd1*l+dd2*dd1*k+dd1*j+i;
          res(i,j,k,l) = vin(elem);
          resv(elem) = res(i,j,k,l);

        }
      }
    }
  }

  return Rcpp::wrap(resv);

}

