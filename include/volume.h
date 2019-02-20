#include <functional>

#include <G4LogicalVolume.hh>
#include <G4VisExtent.hh>

class volume
{
 public:
  std::function<double (double, double, double)> V;
  G4LogicalVolume * log_vol;
  double Xmin, Xmax;
  double Ymin, Ymax;
  double Zmin, Zmax;
 
  volume(double, double, double, double, double, double, std::function<double (double, double, double)>);
  double get_voltage(double, double, double);
  bool is_in_boundary(double, double, double);
};
