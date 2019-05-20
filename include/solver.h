#ifndef SOLVER_H
#define SOLVER_H

#include "boundary.h"
#include "utils.h"
#include "field.h"

#include <vector>

class solver
{
 public:
  int nPointsX;
  int nPointsY;
  int nPointsZ;
  double spacing;

  boundary * bound;
  
  field <double> * bval;
  field <double> * epVal;
  field <double> * sigVal;
  field <bool> * is_b;
  field <double> * Q;
  field <double> * dQdt;
  field <double> * potential;

  field <double> * resid;
  field <double> * a0;
  field <double> * a1;
  field <double> * a2;
  field <double> * a3;
  field <double> * a4;
  field <double> * a5;
  field <double> * a6;

  field <double> * b0;
  field <double> * b1;
  field <double> * b2;
  field <double> * b3;
  field <double> * b4;
  field <double> * b5;
  field <double> * b6;
  
  solver(boundary*, int, double);
  void solve_static();
  void accum_charge();
  void solve_charge();
};

#endif
