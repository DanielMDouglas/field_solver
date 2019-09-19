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
  double spacing = 0.01;

  boundary * bound;

  std::vector <double> x_axis;
  std::vector <double> y_axis;
  std::vector <double> z_axis;

  std::vector <double> x_axis_staggered;
  std::vector <double> y_axis_staggered;
  std::vector <double> z_axis_staggered;

  field <double> * bval;
  field <double> * epVal;
  field <double> * sigVal;
  field <bool> * is_in_volume;
  field <bool> * is_dirichlet;
  field <bool> * is_von_neumann;
  field <bool> * is_boundary;
  field <double> * von_neumann_dV;
  field <double> * Q;
  field <double> * dQdt;
  field <double> * E;
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
  void upscale(int);
  void initialize_axes();
  void initialize_geometry_fields();
  void initialize_permittivity_and_conductivity();
  void solve_static();
  void accum_charge();
  void solve_charge();
  void set_VN();
  void fill_empty_VN();
  void relax();
  int report(int);
};

#endif
