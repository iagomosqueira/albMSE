//////////////////////////////////////////////////////
// multidim arrays as classes ////////////////////////
//////////////////////////////////////////////////////
// R. Hillary CSIRO 2023 /////////////////////////////
//////////////////////////////////////////////////////

#include <iostream>
#include "mdarrclass.h"

using namespace std;

int main ()
{
  arr2d x;
  arr3d y;
  arr4d z;

  x.arr2(1,1);
  x.t[0][0] = 1.;
  cout << x.t[0][0] << endl;
  x(0,0) += 1.;
  cout << x(0,0) << endl; 

  y.arr3(1,1,1);
  y.t[0][0][0] = 1.;
  cout << y.t[0][0][0] << endl;
  y(0,0,0) += 1.;
  cout << y(0,0,0) << endl; 

  z.arr4(1,1,1,1);
  z.t[0][0][0][0] = 1.;
  cout << z.t[0][0][0][0] << endl;
  z(0,0,0,0) += 1.;
  cout << z(0,0,0,0) << endl; 

  return 0;
}
