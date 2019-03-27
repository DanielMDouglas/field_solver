#include "volume.h"
#include "utils.h"

volume::volume(double Xm = 0, double XM = 0,
	       double Ym = 0, double YM = 0,
	       double Zm = 0, double ZM = 0,
	       std::function<double (double, double, double)> voltage_function = constant(0))
{
  Xmin = Xm;
  Xmax = XM;
  Ymin = Ym;
  Ymax = YM;
  Zmin = Zm;
  Zmax = ZM;
  V = voltage_function;

  center = std::vector <double> {0.5*(Xm + XM),
				 0.5*(Ym + YM),
				 0.5*(Zm + ZM)};
}

double volume::get_voltage(double x, double y, double z)
{
  return V(x, y, z);
}

bool volume::is_in_boundary(double x, double y, double z)
{
  if ( ( Xmin <= x ) and
       ( x <= Xmax ) and
       ( Ymin <= y ) and
       ( y <= Ymax ) and
       ( Zmin <= z ) and
       ( z <= Zmax ) ) {
    return true;
  }
  else return false;
}
