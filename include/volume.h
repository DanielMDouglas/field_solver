#include <G4LogicalVolume.hh>
#include <G4VisExtent.hh>

class volume
{
 public:
  double V;
  G4LogicalVolume * log_vol;

  volume(G4LogicalVolume*, double);
  bool is_in_boundary(double, double, double);
  G4VisExtent extent();
};
