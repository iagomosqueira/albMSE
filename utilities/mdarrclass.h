//////////////////////////////////////////////////////
// header classes for multidimensional arrays ////////
//////////////////////////////////////////////////////
// R. Hillary CSIRO 2023 /////////////////////////////
//////////////////////////////////////////////////////

#include <RcppCommon.h>


class arr2d
{

public:
 
  int dd1,dd2;
  double **t; 
  double **arr2(int,int);
  double& operator()(int,int);
  void del();
  
  // Constructors
	arr2d();
  arr2d(SEXP flq_sexp); // Used as intrusive 'as'
  operator SEXP() const; // Used as intrusive 'wrap'
  

};

class arr3d
{

public:
 
  int dd1,dd2,dd3; 
  double ***t; 
  double ***arr3(int,int,int);
  double& operator()(int,int,int);
  void del(); 
  
  // Constructors
	arr3d();
  arr3d(SEXP flq_sexp); // Used as intrusive 'as'
  operator SEXP() const; // Used as intrusive 'wrap'

};

class arr4d
{

public:
 
  int dd1,dd2,dd3,dd4; 
  double ****t; 
  double ****arr4(int,int,int,int);
  double& operator()(int,int,int,int);
  void del(); 
  
// Constructors
// arr4d();
// arr4d(SEXP flq_sexp); // Used as intrusive 'as'
// operator SEXP() const; // Used as intrusive 'wrap'

};

class arr5d
{

public:
 
  int dd1,dd2,dd3,dd4,dd5; 
  double *****t; 
  double *****arr5(int,int,int,int,int);
  double& operator()(int,int,int,int,int);
  void del(); 
  
// Constructors
// arr5d();
// arr5d(SEXP flq_sexp); // Used as intrusive 'as'
// operator SEXP() const; // Used as intrusive 'wrap'

};

// Must be seen after the specialisation or it is not seen by Rcpp types
// See Rcpp-Extending vignette page 2
// Not sure if necessary here as not specialising as or wrap but include it anyway
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

