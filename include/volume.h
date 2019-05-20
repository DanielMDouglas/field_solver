#ifndef VOLUME_H
#define VOLUME_H

#include <functional>
#include <vector>
#include <string>

class volume
{
 public:
  std::function<double (double, double, double)> V;
  double Xmin, Xmax;
  double Ymin, Ymax;
  double Zmin, Zmax;
  std::vector <double> center;
  std::string type;
  std::function<double (double, double, double)> er;
  std::function<double (double, double, double)> sigma;
  bool isSensitive = false;
  
  volume(double, double, double,
	 double, double, double,
	 double, double);
  volume(double, double, double,
	 double, double, double,
	 std::function<double (double, double, double)>,
	 std::function<double (double, double, double)>);
  volume(double, double, double,
	 double, double, double,
	 std::function<double (double, double, double)>);
  double get_voltage(double, double, double);
  bool is_in_boundary(double, double, double);
};

#endif
