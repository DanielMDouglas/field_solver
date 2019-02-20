#include <vector>
#include <string>

#include <volume.h>

class boundary
{
 public:
  volume * volumes[6];
  double Xmin, Xmax = 0;
  double Ymin, Ymax = 0;
  double Zmin, Zmax = 0;
  
  boundary();
  bool is_in_boundary(double, double, double);
  double boundary_value(double, double, double);
};
