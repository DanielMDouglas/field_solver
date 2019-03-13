#ifndef VOLUME_H
#define VOLUME_H

#include <functional>

class volume
{
 public:
  std::function<double (double, double, double)> V;
  double Xmin, Xmax;
  double Ymin, Ymax;
  double Zmin, Zmax;
 
  volume(double, double, double, double, double, double, std::function<double (double, double, double)>);
  double get_voltage(double, double, double);
  bool is_in_boundary(double, double, double);
};

#endif
