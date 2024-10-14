////////////////////////////////////////////////////////// 
// generate eqm popn for MSY estimation IOTC ABC code ////
//////////////////////////////////////////////////////////
// R. Hillary & I. Mosqueira 2023 ////////////////////////
//////////////////////////////////////////////////////////

#include <Rcpp.h>
#include <R.h>
#include "mdarrclass.h"

using namespace Rcpp;

//[[Rcpp::export]]
RcppExport SEXP msypdyn(SEXP dm_,SEXP srec_,SEXP R0_,SEXP hh_,SEXP psi_,SEXP M_,SEXP mata_,SEXP wta_,SEXP sela_,SEXP hinit_) 
{

  IntegerVector dm = as<Rcpp::IntegerVector>(dm_);
  int ns = dm[0];
  int na = dm[1];
  int nf = dm[2];
  int a,f,g,s;

  int srec = as<int>(srec_)-1;
  double R0 = as<double>(R0_);
  double hh = as<double>(hh_); 
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

  double rho = sprf/spr0; 
  double B0 = R0*spr0;
  double alp = 4.*hh/(spr0*(1.-hh));
  double bet = (5.*hh-1.)/(B0*(1.-hh));
  double Bbar = (alp*sprf-1.)/bet > 0. ? (alp*sprf-1.)/bet : 0.;
  double Rbar = Bbar/sprf;

  // scale the exploited population by Rbar

  for(g=0;g<2;g++) 
    for(s=0;s<ns;s++) 
      for(a=0;a<na;a++) N(a,s,g) *= Rbar; 

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

  List res = Rcpp::List::create(Named("rho")=rho,Named("C")=cvec,Named("spr0")=spr0,Named("Bmsy")=Bbar,Named("Rmsy")=Rbar);

  mata.del();
  wta.del();
  hinit.del();
  C.del();
  sela.del();
  N.del(); 
  
  return Rcpp::wrap(res);

}
