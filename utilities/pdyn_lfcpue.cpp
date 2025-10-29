////////////////////////////////////////////////////////// 
// dynamic non-eqm population dynamics for  ABC code /////
// :: ////////////////////////////////////////////////////
// generates model-predicted LF & CPUE data //////////////
//////////////////////////////////////////////////////////
// R. Hillary & I. Mosqueira 2023 ////////////////////////
//////////////////////////////////////////////////////////

#include <Rcpp.h>
#include <R.h>
#include "mdarrclass.h"

using namespace Rcpp;

//[[Rcpp::export]]

RcppExport SEXP pdynlfcpue(SEXP dm_,SEXP srec_,SEXP R0_,SEXP hh_,SEXP psi_,SEXP epsr_,SEXP spr0_,SEXP M_,SEXP mata_,SEXP wta_,SEXP sela_,SEXP Ninit_,SEXP Cb_,SEXP pla_,SEXP fref_) 
{

  IntegerVector dm = as<Rcpp::IntegerVector>(dm_);
  int ny = dm[0];
  int ns = dm[1];
  int na = dm[2];
  int nl = dm[3];
  int nf = dm[4];
  int a,f,g,l,s,y;

  int srec = as<int>(srec_)-1;
  int fref = as<int>(fref_)-1;
  double R0 = as<double>(R0_);
  double hh = as<double>(hh_); 
  double spr0 = as<double>(spr0_);
  double psi = as<double>(psi_);
  double M = as<double>(M_);
  double B0 = R0*spr0;
  double alp = 4.*hh/(spr0*(1.-hh));
  double bet = (5.*hh-1.)/(B0*(1.-hh));

   // vector objects

  NumericVector matv = as<Rcpp::NumericVector>(mata_);
  NumericVector wtv = as<Rcpp::NumericVector>(wta_);
  NumericVector selv = as<Rcpp::NumericVector>(sela_); 
  NumericVector ninitv = as<Rcpp::NumericVector>(Ninit_);
  NumericVector cbv = as<Rcpp::NumericVector>(Cb_); 
  NumericVector epsr = as<Rcpp::NumericVector>(epsr_);
  NumericVector plav = as<Rcpp::NumericVector>(pla_); 

  // reconstruct back into multi-dim arrays
  
  arr3d mata;
  arr3d wta;
  arr3d Cb;
  arr3d H;
  arr4d sela;
  arr4d N;
  arr2d S;
  arr2d I;
  arr4d LF;
  arr4d pla;

  mata.arr3(na,ns,2); // age + season + sex
  wta.arr3(na,ns,2);  // age + season + sex
  Cb.arr3(ny,ns,nf); // year + season + fishery
  H.arr3(ny,ns,nf); // year + season + fishery
  sela.arr4(na,ns,2,nf); // age + season + sex + fishery
  N.arr4(ny,na,ns,2); // year + age + season + sex
  S.arr2(ny,ns); // year + season (female SSB)
  I.arr2(ny,ns); // year + season (for reference fishery)
  LF.arr4(ny,nl,ns,nf); // year + length + season + fishery
  pla.arr4(nl,na,ns,2); // length + age + season + sex

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
      for(y=0;y<ny;y++) {

        elem = ny*ns*f+ny*s+y;
        Cb(y,s,f) = cbv(elem);

      }
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

  for(g=0;g<2;g++) {
    for(s=0;s<ns;s++) {
      for(a=0;a<na;a++) { 
        for(l=0;l<nl;l++) {
      
          elem = nl*na*ns*g+nl*na*s+nl*a+l;
          pla(l,a,s,g) = plav(elem);

        }
      }
    }
  } 

  int spwn = srec == 0 ? ns-1 : srec-1;
  int ysp;
  double hsum,lsum,xsum; 
  
  ////////////////////////
  // initial conditions //
  ////////////////////////

  for(g=0;g<2;g++) {
    for(s=0;s<ns;s++) {
      for(a=0;a<na;a++) {

        elem = na*ns*g+na*s+a;
        N(0,a,s,g) = ninitv(elem);

      }
    }
  }  

  for(s=0;s<ns;s++) {

    S(0,s) = 0.;
    for(a=0;a<na;a++) S(0,s) +=  N(0,a,s,0)*mata(a,s,0)*wta(a,s,0);
     
  } 

  for(f=0;f<nf;f++) {
    for(s=0;s<ns;s++) { 

      xsum=0.;
      for(g=0;g<2;g++) {
       
        for(a=0;a<na;a++) xsum += N(0,a,s,g)*sela(a,s,g,f)*wta(a,s,g);

      }

      hsum = Cb(0,s,f)/xsum;
      H(0,s,f) = hsum < 0.9 ? hsum : 0.9;

      // CPUE

      if(f == fref) {

        xsum = 0.;
        for(g=0;g<2;g++) 
          for(a=0;a<na;a++) xsum += N(0,a,s,g)*sela(a,s,g,f)*wta(a,s,g);

        I(0,s) = xsum;

      }

      // LF data

      for(lsum=0.,l=0;l<nl;l++) {
        
        LF(0,l,s,f) = 0.;
        for(a=0;a<na;a++) {
          for(g=0;g<2;g++) {
            
            LF(0,l,s,f) += N(0,a,s,g)*sela(a,s,g,f)*pla(l,a,s,g);

          }
        }

        lsum += LF(0,l,s,f);
      }
      
      for(l=0;l<nl;l++) LF(0,l,s,f) /= lsum;

    }
  } 

  ////////////////////////////
  // loop through the years //
  ////////////////////////////

  double Rtot;
  for(y=1;y<ny;y++) {

    ysp = srec == 0 ? y-1 : y;

    // season by season

    for(s=0;s<ns;s++) {

      // recruits

      if(s < srec) {

        for(g=0;g<2;g++) N(y,0,s,g) = 0.; 
      
      }

      if(s == srec) {

        Rtot = (alp*S(ysp,spwn)/(1.+bet*S(ysp,spwn)))*exp(epsr(y-1)); 
        N(y,0,s,0) = Rtot*psi; 
        N(y,0,s,1) = Rtot*(1.-psi); 

      }

      if(s > srec) { 
      
        for(g=0;g<2;g++) {

          for(hsum=0.,f=0;f<nf;f++) hsum += H(y,s-1,f)*sela(0,s-1,g,f);
          hsum = hsum < 0.9 ? hsum : 0.9;
          N(y,0,s,g) = N(y,0,s-1,g)*exp(-M)*(1.-hsum);

        }
      }

      // loop through ages

      for(a=1;a<na;a++) {
        for(g=0;g<2;g++) { 
        
          if(s == 0) {
          
            for(hsum=0.,f=0;f<nf;f++) hsum += H(y-1,ns-1,f)*sela(a-1,ns-1,g,f);
            hsum = hsum < 0.9 ? hsum : 0.9; 
            N(y,a,s,g) = N(y-1,a-1,ns-1,g)*exp(-M)*(1.-hsum); 

          } else {

            for(hsum=0.,f=0;f<nf;f++) hsum += H(y,s-1,f)*sela(a,s-1,g,f);
            hsum = hsum < 0.9 ? hsum : 0.9;
            N(y,a,s,g) = N(y,a,s-1,g)*exp(-M)*(1.-hsum); 

          }
        }
      }

      // female SSB

      S(y,s) = 0.;
      for(a=0;a<na;a++) S(y,s) +=  N(y,a,s,0)*mata(a,s,0)*wta(a,s,0);
 
      // harvest rates

      for(f=0;f<nf;f++) {

        xsum=0.;
        for(g=0;g<2;g++) {
       
          for(a=0;a<na;a++) xsum += N(y,a,s,g)*sela(a,s,g,f)*wta(a,s,g);

        }

        hsum = Cb(y,s,f)/xsum;
        H(y,s,f) = hsum < 0.9 ? hsum : 0.9;

        // CPUE

        if(f == fref) {

          xsum = 0.;
          for(g=0;g<2;g++) 
            for(a=0;a<na;a++) xsum += N(y,a,s,g)*sela(a,s,g,f)*wta(a,s,g);

          I(y,s) = xsum;

        } 

        // LF data

        for(lsum=0.,l=0;l<nl;l++) {
        
          LF(y,l,s,f) = 0.;
          for(a=0;a<na;a++) {
            for(g=0;g<2;g++) {
            
              LF(y,l,s,f) += N(y,a,s,g)*sela(a,s,g,f)*pla(l,a,s,g);

            }
          }

          lsum += LF(y,l,s,f);
        }
      
        for(l=0;l<nl;l++) LF(y,l,s,f) /= lsum; 
      }

    }
  }

  // return vectors

  NumericVector svec(ny*ns); // return female SSB object
  NumericVector hvec(ny*ns*nf); // return harvest rate object
  NumericVector nvec(ny*na*ns*2); // return abundances object
  NumericVector ivec(ny*ns); // return CPUE object
  NumericVector lvec(ny*nl*ns*nf); // returb LF object

  for(s=0;s<ns;s++) {
    for(y=0;y<ny;y++) {
      
      elem = ny*s+y;
      svec(elem) = S(y,s);

    }
  }

  for(f=0;f<nf;f++) { 
    for(s=0;s<ns;s++) {
      for(y=0;y<ny;y++) {
      
        elem = ny*ns*f+ny*s+y;
        hvec(elem) = H(y,s,f);

      }
    }
  } 

  for(g=0;g<2;g++) {
    for(s=0;s<ns;s++) {
      for(a=0;a<na;a++) {
        for(y=0;y<ny;y++) { 

          elem = ns*na*ny*g+na*ny*s+ny*a+y;
          nvec(elem) = N(y,a,s,g);
  
        }
      }
    }
  }

  for(s=0;s<ns;s++) {
    for(y=0;y<ny;y++) { 

      elem = ny*s+y;
      ivec(elem) = I(y,s);

    }
  }

  for(f=0;f<nf;f++) { 
    for(s=0;s<ns;s++) {
      for(l=0;l<nl;l++) {
        for(y=0;y<ny;y++) { 

          elem = ny*nl*ns*f+ny*nl*s+ny*l+y;
          lvec(elem) = LF(y,l,s,f);

        }
      }
    }
  }

  List res = Rcpp::List::create(Named("S")=svec,Named("N")=nvec,Named("H")=hvec,Named("I")=ivec,Named("LF")=lvec);

  mata.del();
  wta.del();
  Cb.del();
  H.del();
  sela.del();
  N.del();
  S.del();
  I.del();
  LF.del();
  pla.del(); 

  return Rcpp::wrap(res);
}
