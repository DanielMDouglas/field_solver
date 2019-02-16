#include <string>
#include <vector>

#include <G4LogicalVolume.hh>
#include <G4VisExtent.hh>
#include <G4AffineTransform.hh>

class volume
{
 public:
  double V;
  G4LogicalVolume * log_vol;

  volume(G4LogicalVolume*, double);
  bool is_in_boundary(double, double, double);
  G4VisExtent extent();
};

class boundary
{
 public:
  std::vector <std::string> volNames = {"volFieldShell",
					"volCathode",
					"volLeftPixelPlane",
					"volRightPixelPlane"};
  /* std::vector <std::string> posNames = {"volModuleTopWall_pos_pos", */
  /* 					"volModuleTopWall_pos_pos", */
  /* 					"volModuleTopWall_pos_pos", */
  /* 					"volModuleTopWall_pos_pos"}; */
  std::vector <double> potentials = {0., 0., 0., 0.};
  volume * volumes [4];
  
  boundary(std::string);
  bool is_in_boundary(double, double, double);
  double boundary_value(double, double, double);
  G4VisExtent extent();
};
