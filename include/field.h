#ifndef FIELD_H
#define FIELD_H

#include <functional>
#include <array>

#include "boundary.h"

template <typename T>
class field
{ 
 public:
  int xSize, ySize, zSize;
  T* values;
  std::vector <double> x_space;
  std::vector <double> y_space;
  std::vector <double> z_space;

  field(std::vector <double>,
	std::vector <double>,
	std::vector <double>,
	T);
  field(std::vector <double>,
	std::vector <double>,
	std::vector <double>,
	std::function<T (double, double, double)>);
  field(boundary, int, int, int, std::string);
  field(boundary,
	std::vector <double>,
	std::vector <double>,
	std::vector <double>,
	std::string);
  field(std::string);
  void set(int, int, int, T);
  T get(int, int, int);
  void print_to_file(std::string);
  T interpolate(std::vector <double>);
  std::vector <T> interpolate_grad(std::vector <double>);
};

#include "../src/field.tpp"

double squared_sum(field <double>);
double squared_sum(field <double>, field <bool>);
double squared_diff(field <double>, field <double>);
double squared_diff(field <double>, field <double>, field <bool>);

#endif
