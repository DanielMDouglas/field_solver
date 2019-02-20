#include <G4LogicalVolume.hh>
#include <G4VisExtent.hh>

class volume
{
 public:
  double V;
  G4LogicalVolume * log_vol;
  double Xmin, Xmax;
  double Ymin, Ymax;
  double Zmin, Zmax;
 
  volume(double, double, double, double, double, double, double);
  bool is_in_boundary(double, double, double);
};
