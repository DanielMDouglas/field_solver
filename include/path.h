#ifndef PATH_H
#define PATH_H

#include <string>
#include <vector>

class path
{
 public:
  double arrivalTime;
  double dt;
  std::string fate;
  std::vector <std::vector <double>> steps;

  path(double);
  path(std::string);
  void print_to_file(std::string);
};

#endif
