#ifndef VECTOR_FIELD_H
#define VECTOR_FIELD_H

#include <functional>

#include "scalar_field.h"

class vectorField
{
 public:
  int xSize, ySize, zSize;

  std::vector <std::vector <double>> values;
  std::vector <double> x_space;
  std::vector <double> y_space;
  std::vector <double> z_space;

  vectorField(scalarField, std::string); 
  void set(int, int, int, std::vector <double>);
  std::vector <double> get(int, int, int);
  void print_to_file(std::string);
};

#endif
