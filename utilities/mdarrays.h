// header for multidimensional arrays

// function decs

double **arr2(int,int);
double ***arr3(int,int,int);
double ****arr4(int,int,int,int);

// function defs

double **arr2(int d1,int d2) 
{
 
  double **t = new double*[d1];
  for(int ii=0;ii<d1;ii++) t[ii] = new double[d2];

  return t;

}

double ***arr3(int d1,int d2,int d3)
{
  double ***t = new double**[d1];
  for(int ii=0;ii<d1;ii++) {

    t[ii] = new double*[d2];
    for(int jj=0;jj<d2;jj++) t[ii][jj] = new double[d3];

  }

  return t;
}

double ****arr4(int d1,int d2,int d3,int d4)
{

  double ****t = new double***[d1];
  for(int ii=0;ii<d1;ii++) {

    t[ii] = new double**[d2];
    for(int jj=0;jj<d2;jj++) {

      t[ii][jj] = new double*[d3];
      for(int kk=0;kk<d3;kk++) t[ii][jj][kk] = new double[d4];

    }
  }

  return t;
}

