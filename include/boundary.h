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
  std::vector <std::vector <double>> centers = {std::vector <double> {324.5, 1215.05, 659.95},
						std::vector <double> {3., 1215.1, 660.},
						std::vector <double> {-328.5, 0., 0.},
						std::vector <double> {328.5, 0., 0.}};
  std::vector <double> potentials = {0., -100., 100., 100.};
  volume * volumes [4];
  
  boundary(std::string);
  bool is_in_boundary(double, double, double);
  double boundary_value(double, double, double);
  G4VisExtent extent();
};
