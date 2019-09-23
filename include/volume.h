#ifndef VOLUME_H
#define VOLUME_H

#include <functional>
#include <vector>
#include <string>

class volume
{
 public:
  std::function <double (double, double, double)> V;
  double Efield = 0;
  double Xmin, Xmax;
  double Ymin, Ymax;
  double Zmin, Zmax;
  std::vector <double> center;
  std::string type;
  std::function<double (double, double, double)> er;
  std::function<double (double, double, double)> sigma;
  bool isSensitive = false;

  // dielectric constructor
  volume(double, double, double,
	 double, double, double,
	 double, double);
  volume(double, double, double,
	 double, double, double,
	 std::function<double (double, double, double)>,
	 std::function<double (double, double, double)>);
  // conductor constructor
  volume(double, double, double,
	 double, double, double,
	 std::function<double (double, double, double)>);
  volume(double, double, double,
	 double, double, double,
	 double);
  double get_voltage(double, double, double);
  bool is_in_boundary(double, double, double);
};

#endif
