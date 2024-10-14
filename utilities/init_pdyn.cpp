////////////////////////////////////////////////////////// 
// generate initial population-per-recruit for ABC code //
//////////////////////////////////////////////////////////
// R. Hillary & I. Mosqueira 2023 ////////////////////////
//////////////////////////////////////////////////////////

#include <Rcpp.h>
#include <R.h>
#include "mdarrclass.h"

using namespace Rcpp;

//[[Rcpp::export]]
RcppExport SEXP initpdyn(SEXP dm_,SEXP srec_,SEXP psi_,SEXP M_,SEXP mata_,SEXP wta_,SEXP sela_,SEXP hinit_) 
{

  IntegerVector dm = as<Rcpp::IntegerVector>(dm_);
  int ns = dm[0];
  int na = dm[1];
  int nf = dm[2];
  int a,f,g,s;

  int srec = as<int>(srec_)-1;
  double psi = as<double>(psi_);
  double M = as<double>(M_);

  // vector objects

  NumericVector matv = as<Rcpp::NumericVector>(mata_);
  NumericVector wtv = as<Rcpp::NumericVector>(wta_);
  NumericVector selv = as<Rcpp::NumericVector>(sela_); 
  NumericVector hinitv = as<Rcpp::NumericVector>(hinit_); 

  // reconstruct back into multi-dim arrays
  
  arr3d mata;
  arr3d wta;
  arr2d hinit;
  arr2d C;
  arr4d sela;
  arr3d N;
  
  mata.arr3(na,ns,2); // age + season + sex
  wta.arr3(na,ns,2);  // age + season + sex
  hinit.arr2(ns,nf); // season + fishery
  C.arr2(ns,nf); // season + fishery
  sela.arr4(na,ns,2,nf); // age + season + sex + fishery
  N.arr3(na,ns,2); // age + season + sex
  NumericVector cvec(ns*nf); // return catch object
  NumericVector nvec(na*ns*2); // return numbers object

  int elem;

  for(g=0;g<2;g++) {
    for(s=0;s<ns;s++) {
      for(a=0;a<na;a++) {

        elem = na*ns*g+na*s+a;
        mata(a,s,g) = matv(elem);
        wta(a,s,g) = wtv(elem);

      }
    }
  }

  for(f=0;f<nf;f++) {
    for(s=0;s<ns;s++) {
      
      elem = ns*f+s;
      hinit(s,f) = hinitv(elem);

    }
  } 
   
  for(f=0;f<nf;f++) {
    for(g=0;g<2;g++) {
      for(s=0;s<ns;s++) { 
        for(a=0;a<na;a++) {

          elem = na*ns*2*f+na*ns*g+na*s+a;
          sela(a,s,g,f) = selv(elem);

        }
      }
    }
  }

  // sex set up:
  // 0: female
  // 1: male

  // unexploited equilbrium

  int spwn = srec == 0 ? ns-1 : srec-1;
  double hsum;
  double spr0;

  for(g=0;g<2;g++) {
    for(s=0;s<ns;s++) {

      if(s < srec) N(0,s,g) = 0.;
      if(s == srec) {

        if(g == 0) N(0,s,g) = psi;
        if(g == 1) N(0,s,g) = 1.-psi;

      }
      if(s > srec) N(0,s,g) = N(0,s-1,g)*exp(-M);
      
    }

    for(a=1;a<na;a++) {

      for(s=0;s<ns;s++) {

        if(s == 0) {

          N(a,0,g) = N(a-1,ns-1,g)*exp(-M);

        } else {

          N(a,s,g) = N(a,s-1,g)*exp(-M);
      
        }
      }
    }
  }

  for(spr0=0.,a=0;a<na;a++) spr0 += N(a,spwn,0)*mata(a,spwn,0)*wta(a,spwn,0); 

  // exploited equilibrium

  double sprf;

  for(g=0;g<2;g++) {
    for(s=0;s<ns;s++) {

      if(s < srec) N(0,s,g) = 0.;
      if(s == srec) {

        if(g == 0) N(0,s,g) = psi;
        if(g == 1) N(0,s,g) = 1.-psi;

      } 

      if(s > srec) {
        
        for(hsum=0.,f=0;f<nf;f++) hsum += hinit(s-1,f)*sela(0,s-1,g,f);
        hsum = hsum > 0.9 ? 0.9 : hsum;
        N(0,s,g) = N(0,s-1,g)*exp(-M)*(1.-hsum);

      }
    }

    for(a=1;a<na;a++) {

      for(s=0;s<ns;s++) {

        if(s == 0) {

          for(hsum=0.,f=0;f<nf;f++) hsum += hinit(ns-1,f)*sela(a-1,ns-1,g,f); 
          hsum = hsum > 0.9 ? 0.9 : hsum;
          N(a,0,g) = N(a-1,ns-1,g)*exp(-M)*(1.-hsum);

        } else {

          for(hsum=0.,f=0;f<nf;f++) hsum += hinit(s-1,f)*sela(a,s-1,g,f); 
          hsum = hsum > 0.9 ? 0.9 : hsum;
          N(a,s,g) = N(a,s-1,g)*exp(-M)*(1.-hsum);
      
        }
      }
    }
  }

  for(sprf=0.,a=0;a<na;a++) sprf += N(a,spwn,0)*mata(a,spwn,0)*wta(a,spwn,0); 

  for(g=0;g<2;g++) {
    for(s=0;s<ns;s++) { 
      for(a=0;a<na;a++) {

        elem = na*ns*g+na*s+a;
        nvec(elem) = N(a,s,g);

      }
    }
  }

  // catch biomass-by-fleet

  for(s=0;s<ns;s++) {
    for(f=0;f<nf;f++) {

      C(s,f) = 0.; 
      for(g=0;g<2;g++) { 
       
        for(a=0;a<na;a++) C(s,f) += N(a,s,g)*wta(a,s,g)*sela(a,s,g,f)*hinit(s,f);

      }
    }
  }

  for(f=0;f<nf;f++) {
    for(s=0;s<ns;s++) {

       elem = ns*f+s;
       cvec(elem) = C(s,f);

    }
  }


  double rho = sprf/spr0;
  List res = Rcpp::List::create(Named("rho")=rho,Named("C")=cvec,Named("N")=nvec,Named("spr0")=spr0);

  // delete arrays

  mata.del();
  wta.del();
  hinit.del();
  C.del();
  sela.del();
  N.del();
  
  return Rcpp::wrap(res);

}
