#ifndef SCALAR_FIELD_H
#define SCALAR_FIELD_H

#include <functional>
#include <array>

#include "boundary.h"

template <typename T>
class scalarField
{
 public:
  int xSize, ySize, zSize;
  // const static int size = 21952000;
  
  /* std::vector <double> values; */
  // std::array <T, size> values;
  // std::array <T, 0> * values;
  T* values;
  std::vector <double> x_space;
  std::vector <double> y_space;
  std::vector <double> z_space;

  scalarField(std::vector <double>,
  	      std::vector <double>,
  	      std::vector <double>,
  	      T);
  scalarField(std::vector <double>,
	      std::vector <double>,
	      std::vector <double>,
	      std::function<T (double, double, double)>);
  scalarField(boundary, int, int, int, std::string);
  scalarField(std::string);
  void set(int, int, int, T);
  T get(int, int, int);
  void print_to_file(std::string);
  T interpolate(std::vector <double>);
  std::vector <T> interpolate_grad(std::vector <double>);
};

#include "../src/scalar_field.tpp"

double squared_sum(scalarField <double>);
double squared_sum(scalarField <double>, scalarField <bool>);
double squared_diff(scalarField <double>, scalarField <double>);
double squared_diff(scalarField <double>, scalarField <double>, scalarField <bool>);

#endif
