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

  path(double);
  path(std::string);
  void shift(std::vector <double>);
  void print_to_file(std::string);
};

#endif
