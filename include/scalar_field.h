#ifndef SCALAR_FIELD_H
#define SCALAR_FIELD_H

#include <functional>

#include "boundary.h"

class scalarField
{
 public:
  int xSize, ySize, zSize;

  std::vector <double> values;
  std::vector <double> x_space;
  std::vector <double> y_space;
  std::vector <double> z_space;

  scalarField(std::vector <double>,
	      std::vector <double>,
	      std::vector <double>,
	      double);
  scalarField(std::vector <double>,
	      std::vector <double>,
	      std::vector <double>,
	      std::function<double (double, double, double)>);
  scalarField(boundary, int, int, int, std::string);
  scalarField(std::string);
  void set(int, int, int, double);
  double get(int, int, int);
  void print_to_file(std::string);
  double interpolate(std::vector <double>);
};

double squared_diff(scalarField, scalarField);
double squared_diff(scalarField, scalarField, scalarField);

#endif
