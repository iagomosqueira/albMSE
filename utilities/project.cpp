////////////////////////////////////////////////////////// 
// dynamic non-eqm population dynamics for  ABC code /////
// :: ////////////////////////////////////////////////////
// has MCMC iters directly included //////////////////////
//////////////////////////////////////////////////////////
// R. Hillary & I. Mosqueira 2025 ////////////////////////
//////////////////////////////////////////////////////////

#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include "mdarrclass.h"

using namespace Rcpp;
using namespace std;

//[[Rcpp::export]]
RcppExport SEXP projv2(SEXP dm_,SEXP srec_,SEXP R0_,SEXP hh_,SEXP psi_,SEXP sigmar_,SEXP spr0_,SEXP M_,SEXP mata_,SEXP wta_,SEXP sela_,SEXP Ninit_,SEXP Cb_,SEXP q_,SEXP fref_) 
{

  IntegerVector dm = as<Rcpp::IntegerVector>(dm_);
  int ny = dm[0];
  int ns = dm[1];
  int na = dm[2];
  int nf = dm[3];
  int nits = dm[4];
  int a,f,g,i,l,s,y;

  int srec = as<int>(srec_)-1;
  int fref = as<int>(fref_)-1;

   // vector objects

  NumericVector matv = as<Rcpp::NumericVector>(mata_);
  NumericVector wtv = as<Rcpp::NumericVector>(wta_);
  NumericVector selv = as<Rcpp::NumericVector>(sela_); 
  NumericVector ninitv = as<Rcpp::NumericVector>(Ninit_);
  NumericVector cbv = as<Rcpp::NumericVector>(Cb_); 
  NumericVector sigmar = as<Rcpp::NumericVector>(sigmar_); 
  NumericVector R0 = as<Rcpp::NumericVector>(R0_);
  NumericVector hh = as<Rcpp::NumericVector>(hh_);
  NumericVector spr0 = as<Rcpp::NumericVector>(spr0_); 
  NumericVector M = as<Rcpp::NumericVector>(M_); 
  NumericVector qv = as<Rcpp::NumericVector>(q_);
  NumericVector B0(nits);
  NumericVector alp(nits);
  NumericVector bet(nits); 
  double psi = as<double>(psi_);  

  for(i=0;i<nits;i++) {

    B0(i) = R0(i)*spr0(i);
    alp(i) = 4.*hh(i)/(spr0(i)*(1.-hh(i))); 
    bet(i) = (5.*hh(i)-1.)/(B0(i)*(1.-hh(i))); 

  }  

  // reconstruct back into multi-dim arrays
  
  arr2d q;
  arr3d mata;
  arr3d wta;
  arr4d Cb;
  arr4d H;
  arr5d sela;
  arr5d N;
  arr3d S;
  arr3d I;

  q.arr2(ns,nits); // season + iteration
  mata.arr3(na,ns,2); // age + season + sex
  wta.arr3(na,ns,2);  // age + season + sex
  Cb.arr4(ny,ns,nf,nits); // year + season + fishery + iter
  H.arr4(ny,ns,nf,nits); // year + season + fishery + iter
  sela.arr5(na,ns,2,nf,nits); // age + season + sex + fishery + iter
  N.arr5(ny,na,ns,2,nits); // year + age + season + sex + iter
  S.arr3(ny,ns,nits); // year + season + iter (female SSB)
  I.arr3(ny,ns,nits); // year + season + iter (for reference fishery)
  
  int elem;

  for(i=0;i<nits;i++) {
    for(s=0;s<ns;s++) {

      elem = ns*i+s;
      q(s,i) = qv(elem);

    }
  }

  for(i=0;i<nits;i++) {
    for(g=0;g<2;g++) {
      for(s=0;s<ns;s++) {
        for(a=0;a<na;a++) {

          elem = na*ns*g+na*s+a;
          mata(a,s,g) = matv(elem);
          wta(a,s,g) = wtv(elem);

        }
      }
    }
  }

  for(i=0;i<nits;i++) {
    for(f=0;f<nf;f++) {
      for(s=0;s<ns;s++) {
        for(y=0;y<ny;y++) {

          elem = ny*ns*nf*i+ny*ns*f+ny*s+y;
          Cb(y,s,f,i) = cbv(elem);

        }
      }
    }
  }

  for(i=0;i<nits;i++) { 
    for(f=0;f<nf;f++) {
      for(g=0;g<2;g++) {
        for(s=0;s<ns;s++) { 
          for(a=0;a<na;a++) {

            elem = na*ns*2*nf*i+na*ns*2*f+na*ns*g+na*s+a;
            sela(a,s,g,f,i) = selv(elem);

          }
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

  for(i=0;i<nits;i++) { 
    
    for(g=0;g<2;g++) {
      for(s=0;s<ns;s++) {
        for(a=0;a<na;a++) {

          elem = na*ns*2*i+na*ns*g+na*s+a;
          N(0,a,s,g,i) = ninitv(elem);

        }
      }
    }
  
    for(s=0;s<ns;s++) {

      S(0,s,i) = 0.;
      for(a=0;a<na;a++) S(0,s,i) +=  N(0,a,s,0,i)*mata(a,s,0)*wta(a,s,0);
     
    } 

    for(f=0;f<nf;f++) {
      for(s=0;s<ns;s++) { 

        xsum=0.;
        for(g=0;g<2;g++) 
          for(a=0;a<na;a++) xsum += N(0,a,s,g,i)*sela(a,s,g,f,i)*wta(a,s,g);

        hsum = Cb(0,s,f,i)/xsum;
        H(0,s,f,i) = hsum < 0.9 ? hsum : 0.9;

        // CPUE

        if(f == fref) {

          xsum = 0.;
          for(g=0;g<2;g++) 
            for(a=0;a<na;a++) xsum += N(0,a,s,g,i)*sela(a,s,g,f,i)*wta(a,s,g);

          I(0,s,i) = q(s,i)*xsum;

        }
      }
    }
  }

  /////////////////////////////////////////
  // loop through the years & iterations //
  /////////////////////////////////////////

  double Rtot,epsr;
  for(i=0;i<nits;i++) {
    for(y=1;y<ny;y++) {

      ysp = srec == 0 ? y-1 : y;

      // season by season

      for(s=0;s<ns;s++) {

        // recruits

        if(s < srec) {

          for(g=0;g<2;g++) N(y,0,s,g,i) = 0.; 
      
        }

        if(s == srec) {

          epsr = R::rnorm(0.,sigmar(i));
          Rtot = (alp(i)*S(ysp,spwn,i)/(1.+bet(i)*S(ysp,spwn,i)))*exp(epsr); 
          N(y,0,s,0,i) = Rtot*psi; 
          N(y,0,s,1,i) = Rtot*(1.-psi); 

        }

        if(s > srec) { 
      
          for(g=0;g<2;g++) {

            for(hsum=0.,f=0;f<nf;f++) hsum += H(y,s-1,f,i)*sela(0,s-1,g,f,i);
            hsum = hsum < 0.9 ? hsum : 0.9;
            N(y,0,s,g,i) = N(y,0,s-1,g,i)*exp(-M(i))*(1.-hsum);

          }
        }

        // loop through ages

        for(a=1;a<na;a++) {
          for(g=0;g<2;g++) { 
        
            if(s == 0) {
          
              for(hsum=0.,f=0;f<nf;f++) hsum += H(y-1,ns-1,f,i)*sela(a-1,ns-1,g,f,i);
              hsum = hsum < 0.9 ? hsum : 0.9; 
              N(y,a,s,g,i) = N(y-1,a-1,ns-1,g,i)*exp(-M(i))*(1.-hsum); 

            } else {

              for(hsum=0.,f=0;f<nf;f++) hsum += H(y,s-1,f,i)*sela(a,s-1,g,f,i);
              hsum = hsum < 0.9 ? hsum : 0.9;
              N(y,a,s,g,i) = N(y,a,s-1,g,i)*exp(-M(i))*(1.-hsum); 

            }
          }
        }

        // female SSB

        S(y,s,i) = 0.;
        for(a=0;a<na;a++) S(y,s,i) +=  N(y,a,s,0,i)*mata(a,s,0)*wta(a,s,0);
 
        // harvest rates

        for(f=0;f<nf;f++) {

          xsum=0.;
          for(g=0;g<2;g++) 
            for(a=0;a<na;a++) xsum += N(y,a,s,g,i)*sela(a,s,g,f,i)*wta(a,s,g);

          hsum = Cb(y,s,f,i)/xsum;
          H(y,s,f,i) = hsum < 0.9 ? hsum : 0.9;

          // CPUE

          if(f == fref) {

            xsum = 0.;
            for(g=0;g<2;g++) 
              for(a=0;a<na;a++) xsum += N(y,a,s,g,i)*sela(a,s,g,f,i)*wta(a,s,g);

            I(y,s,i) = q(s,i)*xsum;

          } 
        }
      }
    }
  }

  // return vectors

  NumericVector svec(ny*ns*nits); // return female SSB object
  NumericVector hvec(ny*ns*nf*nits); // return harvest rate object
  NumericVector nvec(ny*na*ns*2*nits); // return numbers at age object
  NumericVector ivec(ny*ns*nits); // return CPUE object

  for(i=0;i<nits;i++) {
    for(s=0;s<ns;s++) {
      for(y=0;y<ny;y++) {
      
        elem = ny*ns*i+ny*s+y;
        svec(elem) = S(y,s,i);

      }
    }
  }

  for(i=0;i<nits;i++) { 
    for(f=0;f<nf;f++) { 
      for(s=0;s<ns;s++) {
        for(y=0;y<ny;y++) {
      
          elem = ny*ns*nf*i+ny*ns*f+ny*s+y;
          hvec(elem) = H(y,s,f,i);

        }
      }
    }
  }

  for(i=0;i<nits;i++) {
    for(g=0;g<2;g++) {
      for(s=0;s<ns;s++) {
        for(a=0;a<na;a++) {
          for(y=0;y<ny;y++) { 

            elem = ny*na*ns*2*i+ns*na*ny*g+na*ny*s+ny*a+y;
            nvec(elem) = N(y,a,s,g,i);
  
          }
        }
      }
    }
  }

  for(i=0;i<nits;i++) { 
    for(s=0;s<ns;s++) {
      for(y=0;y<ny;y++) { 

        elem = ny*ns*i+ny*s+y;
        ivec(elem) = I(y,s,i);

      }
    }
  }

  List res = Rcpp::List::create(Named("S")=svec,Named("N")=nvec,Named("H")=hvec,Named("I")=ivec);

  q.del();
  mata.del();
  wta.del();
  Cb.del();
  H.del();
  sela.del();
  N.del();
  S.del();
  I.del();

  return Rcpp::wrap(res);
}
