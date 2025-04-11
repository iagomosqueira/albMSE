// Farmed out implementation of array classes to separate file


#include "mdarrclass.h"


double** arr2d::arr2(int d1,int d2)
{
  dd1 = d1;
  dd2 = d2;
  t = new double*[d1];
  for(int ii=0;ii<d1;ii++) t[ii] = new double[d2];

  return t;

}

double& arr2d::operator()(int i1,int i2) 
{
  return(t[i1][i2]);
}

void arr2d::del()
{
   
  for(int ii=0;ii<dd1;ii++) delete[] t[ii];

  delete[] t;

}


// Do nothing constructor
arr2d::arr2d(){
//  Rprintf("In arr2d basic constructor.\n");
}

// <as> - pass in FLQuant as SEXP and make an arr2d from it
arr2d::arr2d(SEXP flq_sexp){
  Rprintf("In arr2d intrinsic as.\n");
  // Pull FLQ apart to get the bits we need: data and dim
  Rcpp::S4 flq_s4(flq_sexp);
  // Pulling out the .Data slot.
  // Takes time but cannot get round this extraction 
  // We need it to be a NumericVector (so we can get dim and dimnames)
  Rcpp::NumericVector data_nv = flq_s4.slot(".Data"); 
  // Other FLQuant bits
  std::vector<unsigned int> dim = Rcpp::as<std::vector<unsigned int>>(data_nv.attr("dim"));
  // Looking at 2D arrays in the pdyn_lfcpue.cpp script, dims of 2D arrays are year and season
  // But this is not true of the Xd arrays (some are age, season, sex, some are year, season, fishery)
  // So just use first two dims of FLQ - but will need to watch permutation of dims.
  dd1 = dim[0];
  dd2 = dim[1];
  // Make the arr2 array of correct size
  arr2(dd1, dd2);
  // Fill it up
  // All 2D arrays in Rich's code are results of calcs, not args passed in
  // So this whole <as> is unnecessary - good example though - can be used for 3D etc
  
  for (int d2=0; d2<dd2; d2++){
    for (int d1=0; d1<dd1; d1++){
      int elem = dd1 * d2 + d1;
      t[d1][d2] = data_nv(elem);
    }
  }
}

// <wrap>  
// Put an arr2D back into an FLQuant - missing dimnames though
// Vector of doubles
arr2d::operator SEXP() const{
  Rprintf("In arr2D wrap\n");
  // Make a new SEXP FLQuant
  Rcpp::S4 flq_s4("FLQuant");
  // Make and fill the NumericVector that will be the 'data' slot 
  // Fill up array (really a vector)
  int size = dd1 * dd2;
  Rcpp::NumericVector data_nv(size);
  for (int d2=0; d2<dd2; d2++){
    for (int d1=0; d1<dd1; d1++){
      int elem = dd1*d2+d1;
      data_nv(elem) = t[d1][d2];
    }
  }
  
  // Delete t to avoid memory leaks
  for(int ii=0;ii<dd1;ii++) delete[] t[ii];
  delete[] t;
  
  // Sort out dims and dimnames
  std::vector<int> dim(6, 1);
  dim[0] = dd1;
  dim[1] = dd2;
  data_nv.attr("dim") = dim;
  // Have to give it something for dimnames - empty for now
  data_nv.attr("dimnames") = Rcpp::List(6);
  flq_s4.slot(".Data") = data_nv;
  //flq_s4.slot("units") = units;
  return Rcpp::wrap(flq_s4);
}




/* ---------------------------------------------------------------- */

double*** arr3d::arr3(int d1,int d2,int d3)
{
   
  dd1 = d1;
  dd2 = d2; 
  dd3 = d3;
  t = new double**[d1];
  for(int ii=0;ii<d1;ii++) {

    t[ii] = new double*[d2];
    for(int jj=0;jj<d2;jj++) t[ii][jj] = new double[d3];

  }

  return t;

}

double& arr3d::operator()(int i1,int i2,int i3) 
{
  return(t[i1][i2][i3]);
}

void arr3d::del()
{

  for(int ii=0;ii<dd1;ii++) 
    for(int jj=0;jj<dd2;jj++) delete[] t[ii][jj];

  for(int ii=0;ii<dd1;ii++) delete[] t[ii];
 
  delete[] t;

}

// Do nothing constructor
arr3d::arr3d(){
 // Rprintf("In arr3d basic constructor.\n");
}

// <as> - pass in FLQuant as SEXP and make an arr3d from it
arr3d::arr3d(SEXP flq_sexp){
  Rprintf("In arr3d intrinsic as.\n");
  // Pull FLQ apart to get the bits we need: data and dim
  Rcpp::S4 flq_s4(flq_sexp);
  // Pulling out the .Data slot.
  Rcpp::NumericVector data_nv = flq_s4.slot(".Data"); 
  // Other FLQuant bits
  std::vector<unsigned int> dim = Rcpp::as<std::vector<unsigned int>>(data_nv.attr("dim"));
  
  // Just use first three dims of FLQ - will need to permutate in R before sending to C++
  dd1 = dim[0];
  dd2 = dim[1];
  dd3 = dim[2];
  // Make the arr3 array of correct size
  arr3(dd1, dd2, dd3); 
  
  // Fill it up
  for (int d3=0; d3<dd3; d3++){
    for (int d2=0; d2<dd2; d2++){
      for (int d1=0; d1<dd1; d1++){
        int elem = dd1*dd2*d3 + dd1*d2 + d1;
        t[d1][d2][d3] = data_nv(elem);
      }
    }
  }
}

// <wrap>  
// Put an arr3D back into an FLQuant - missing dimnames though
// Vector of doubles
arr3d::operator SEXP() const{
  Rprintf("In arr3D wrap\n");
  // Make a new SEXP FLQuant
  Rcpp::S4 flq_s4("FLQuant");
  // Make and fill the NumericVector that will be the 'data' slot 
  
  // Fill up array (really a vector)
  int size = dd1 * dd2 * dd3; 
  Rcpp::NumericVector data_nv(size);
  for (int d3=0; d3<dd3; d3++){
    for (int d2=0; d2<dd2; d2++){
      for (int d1=0; d1<dd1; d1++){
        int elem = dd1*dd2*d3 + dd1*d2 + d1;
        data_nv(elem) = t[d1][d2][d3];
      }
    }
  }
  
  // Delete t to avoid memory leaks
  for(int ii=0;ii<dd1;ii++) 
    for(int jj=0;jj<dd2;jj++) delete[] t[ii][jj];

  for(int ii=0;ii<dd1;ii++) delete[] t[ii];
 
  delete[] t;
  
  // Apply dims and dimnames
  std::vector<int> dim(6, 1);
  dim[0] = dd1;
  dim[1] = dd2;
  dim[2] = dd3;
  // Fill the slots
  data_nv.attr("dim") = dim;
  // Have to give it something for dimnames - empty for now
  data_nv.attr("dimnames") = Rcpp::List(6);
  flq_s4.slot(".Data") = data_nv;
  //flq_s4.slot("units") = units;
  return Rcpp::wrap(flq_s4);
}




/* ---------------------------------------------------------------- */

double**** arr4d::arr4(int d1,int d2,int d3,int d4)
{

  dd1 = d1;
  dd2 = d2; 
  dd3 = d3;
  dd4 = d4;
  t = new double***[d1];
  for(int ii=0;ii<d1;ii++) {

    t[ii] = new double**[d2];
    for(int jj=0;jj<d2;jj++) {

      t[ii][jj] = new double*[d3];
      for(int kk=0;kk<d3;kk++) t[ii][jj][kk] = new double[d4];

    }
  }

  return t;
}

double& arr4d::operator()(int i1,int i2,int i3,int i4) 
{
  return(t[i1][i2][i3][i4]);
}

void arr4d::del()
{
   
  for(int ii=0;ii<dd1;ii++) 
    for(int jj=0;jj<dd2;jj++)  
      for(int kk=0;kk<dd3;kk++) delete[] t[ii][jj][kk];  

  for(int ii=0;ii<dd1;ii++) 
    for(int jj=0;jj<dd2;jj++) delete[] t[ii][jj];

  for(int ii=0;ii<dd1;ii++) delete[] t[ii]; 

  delete[] t;

}


double***** arr5d::arr5(int d1,int d2,int d3,int d4,int d5)
{

  dd1 = d1;
  dd2 = d2; 
  dd3 = d3;
  dd4 = d4;
  dd5 = d5;
  t = new double****[d1];
  for(int ii=0;ii<d1;ii++) {

    t[ii] = new double***[d2];
    for(int jj=0;jj<d2;jj++) {

      t[ii][jj] = new double**[d3];
      for(int kk=0;kk<d3;kk++) {

        t[ii][jj][kk] = new double*[d4];
        for(int ll=0;ll<d4;ll++) t[ii][jj][kk][ll] = new double[d5];
      
      }
    }
  }

  return t;
}

double& arr5d::operator()(int i1,int i2,int i3,int i4,int i5) 
{
  return(t[i1][i2][i3][i4][i5]);
}

void arr5d::del()
{
   
  for(int ii=0;ii<dd1;ii++) 
    for(int jj=0;jj<dd2;jj++)  
      for(int kk=0;kk<dd3;kk++)
        for(int ll=0;ll<dd4;ll++) delete[] t[ii][jj][kk][ll];  

  for(int ii=0;ii<dd1;ii++) 
    for(int jj=0;jj<dd2;jj++)  
      for(int kk=0;kk<dd3;kk++) delete[] t[ii][jj][kk]; 

  for(int ii=0;ii<dd1;ii++) 
    for(int jj=0;jj<dd2;jj++) delete[] t[ii][jj];

  for(int ii=0;ii<dd1;ii++) delete[] t[ii]; 

  delete[] t;

}

