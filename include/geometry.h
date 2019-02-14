#include <string>
#include <vector>

#include <G4LogicalVolume.hh>

class volume
{
 public:
  double V;
  G4LogicalVolume * log_vol;

  volume(G4LogicalVolume*, double);
  bool is_in_boundary(double, double, double);
};

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
};
