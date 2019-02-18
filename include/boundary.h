#include <vector>
#include <string>

#include <G4VisExtent.hh>

#include <volume.h>

class boundary
{
 public:
  std::vector <std::string> volNames = {"volFieldShell",
					"volCathode",
					"volLeftPixelPlane",
					"volRightPixelPlane"};
  std::vector <double> potentials = {0., 0., 0., 0.};
  volume * volumes [4];
  
  boundary(std::string);
  bool is_in_boundary(double, double, double);
  double boundary_value(double, double, double);
  G4VisExtent extent();
};
