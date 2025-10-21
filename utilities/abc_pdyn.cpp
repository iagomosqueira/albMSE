//////////////////////////////////////////////////////
// abcpdyn.cpp ///////////////////////////////////////
//////////////////////////////////////////////////////
// R. Hillary CSIRO 2022 /////////////////////////////
//////////////////////////////////////////////////////

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <iostream>
#include <R.h>
#include "assert.h"

using namespace Rcpp;

/*
 * Project forward in time for fixed catch, get CPUE and CAL
 */
//[[Rcpp::export]]
RcppExport SEXP abcpdyn(SEXP  B0_,SEXP hh_,SEXP M_,SEXP hinit_,SEXP cshare_,SEXP mat_,SEXP wt_,SEXP epsr_,SEXP sigmar_,SEXP pla_,SEXP self_,SEXP dm_,SEXP C_,SEXP fidx_)
{
  IntegerVector dm = as<Rcpp::IntegerVector>(dm_);
  int ny = dm[0];
  int na = dm[1];
  int nf = dm[2];
  int nl = dm[3];
  int fidx = as<int>(fidx_)-1;
  int a,f,i,l,y;

  // read ins

  NumericVector M = as<Rcpp::NumericVector>(M_);
  NumericVector wt = as<Rcpp::NumericVector>(wt_); 
  NumericVector mat = as<Rcpp::NumericVector>(mat_);
  NumericMatrix sel = as<Rcpp::NumericMatrix>(self_); 
  NumericVector cshare = as<Rcpp::NumericVector>(cshare_); 
  NumericVector epsr = as<Rcpp::NumericVector>(epsr_); 
  NumericMatrix C = as<Rcpp::NumericMatrix>(C_); 
  NumericMatrix pla = as<Rcpp::NumericMatrix>(pla_); 
  double B0 = as<double>(B0_);
  double hh = as<double>(hh_);
  double hinit = as<double>(hinit_);
  double sigmar = as<double>(sigmar_);
  double xisum,xsum;
  double alp,bet,SPR0,SPR,phi,Rbar,Sbar;

  // memory allocation

  NumericMatrix N(na,ny);
  NumericVector SSB(ny);
  NumericVector Ninit(na);
  NumericMatrix X(nf,ny);
  NumericMatrix LF(nl,ny); 
  NumericVector hy(ny); 
  NumericVector I(ny);

  // S-R pars
  
  Ninit(0) = 1.; // per recruit
  for(a=1;a<na;a++) Ninit(a) = Ninit(a-1)*exp(-M(a-1));
  Ninit(na-1) /= (1.-exp(-M(na-1)));
  for(SPR0=0.,a=0;a<na;a++) SPR0 += Ninit(a)*wt(a)*mat(a);
  alp = 4.*hh/(SPR0*(1.-hh));
  bet = (5.*hh-1.)/(B0*(1.-hh));

  // init popn @ exploited equilibirium

  Ninit(0) = 1.;
  for(a=1;a<na;a++) {
    
    for(xisum=0.,f=0;f<nf;f++) xisum += hinit*sel(f,a-1)*cshare(f);
    Ninit(a) = Ninit(a-1)*exp(-M(a-1))*(1.-xisum);

  }

  for(xisum=0.,f=0;f<nf;f++) xisum += hinit*sel(f,na-1)*cshare(f);
  Ninit(na-1) /= (1.-exp(-M(na-1))*(1.-xisum));
  for(SPR=0.,a=0;a<na;a++) SPR += Ninit(a)*wt(a)*mat(a); 
  phi = SPR/SPR0;
  Sbar = (alp*SPR-1.)/bet;
  Sbar = Sbar > 0. ? Sbar : 0.;
  Rbar = Sbar/SPR;
  N(0,0) = Rbar*exp(epsr(0));
  for(a=1;a<na;a++) N(a,0) = Rbar*Ninit(a);
  for(SSB(0)=0.,a=0;a<na;a++) SSB(0) += N(a,0)*wt(a)*mat(a);

  // set up the fixed catch-derived exploitation rates

  hy(0) = 0.;
  for(f=0;f<nf;f++) {
      
    xsum = 0.; 
    for(a=0;a<na;a++) xsum += N(a,0)*sel(f,a)*wt(a);  
    X(f,0) = cshare(f)*C(f,0)/xsum;
    X(f,0) = X(f,0) < 0.9 ? X(f,0) : 0.9; 
    hy(0) += X(f,0);

  }  

  // loop through the years

  for(y=1;y<ny;y++) {

    // recruitment

      N(0,y) = alp*SSB(y-1)/(1.+bet*SSB(y-1))*exp(epsr(y)); 

    // rest of the ages

    for(a=1;a<na;a++) {

      for(xisum=0.,f=0;f<nf;f++) xisum += X(f,y-1)*sel(f,a-1);
      xisum = xisum < 0.95 ? xisum : 0.95;

      N(a,y) = N(a-1,y-1)*exp(-M(a-1))*(1.-xisum);

    }

    // plus group

    for(xisum=0.,f=0;f<nf;f++) xisum += X(f,y-1)*sel(f,na-1); 
    xisum = xisum < 0.99 ? xisum : 0.99; 
    N(na-1,y) += N(na-1,y-1)*exp(-M(na-1))*(1.-xisum);

    // SSB

    SSB(y) = 0.;
    for(a=0;a<na;a++) SSB(y) += N(a,y)*wt(a)*mat(a); 

    // set up the fixed catch-derived exploitation rates

    hy(y) = 0.;
    for(f=0;f<nf;f++) {
      
      xsum = 0.; 
      for(a=0;a<na;a++) xsum += N(a,y)*sel(f,a)*wt(a);
      X(f,y) = cshare(f)*C(f,y)/xsum;
      X(f,y) = X(f,y) < 0.9 ? X(f,y) : 0.9;
      hy(y) += X(f,y);

    }

  }

  ///////////////////////////////////////
  // calculate predicted observed data //
  ///////////////////////////////////////

  double lfsum;
  for(y=0;y<ny;y++) {

    // CPUE for index fleet

    for(I(y)=0.,a=0;a<na;a++) I(y) += N(a,y)*sel(fidx,a);

    // LF for index fleet
    
    for(lfsum=0.,l=0;l<nl;l++) {

      for(LF(l,y) = 0.,a=0;a<na;a++) LF(l,y) += N(a,y)*sel(fidx,a)*pla(l,a);
      lfsum += LF(l,y);
     
    }
     
    for(l=0;l<nl;l++) LF(l,y) /= lfsum;

  }

  List res = Rcpp::List::create(Named("N")=N,Named("SSB")=SSB,Named("hy")=hy,Named("I")=I,Named("LF")=LF);
 
  return Rcpp::wrap(res); 
}
