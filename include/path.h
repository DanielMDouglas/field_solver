#ifndef PATH_H
#define PATH_H

#include <string>
#include <vector>

#include "utils.h"

class path
{
 public:
  double arrivalTime;
  double dt;
  std::string fate;
  std::vector <std::vector <double>> pos;
  std::vector <std::vector <double>> vel;
  std::vector <std::vector <double>> E;

  path(double);
  path(std::string);
  path * copy();
  void reflectX();
  void reflectY();
  void reflectZ();
  void shift(std::vector <double>);
  void print_to_file(std::string);
};

#endif
