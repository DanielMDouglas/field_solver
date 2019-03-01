#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <vector>
#include <string>

#include "volume.h"
#include "utils.h"

class boundary
{
 public:
  const static int maxNVolumes = 200;
  int nVolumes = 0;
  volume * volumes[maxNVolumes];
  double Xmin, Xmax = 0;
  double Ymin, Ymax = 0;
  double Zmin, Zmax = 0;
  
  boundary();
  void make_pads(double, double,
		 double, double,
		 double, double);
  void make_wires(double, double,
		  double, double,
		  double, double);
  void make_field_cage(double, double,
		       double, double,
		       double, double,
		       double,
		       double,
		       double);
  bool is_in_boundary(double, double, double);
  double boundary_value(double, double, double);
};

#endif
