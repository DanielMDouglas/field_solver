#include "volume.h"

volume::volume(double Xm, double XM,
	       double Ym, double YM,
	       double Zm, double ZM,
	       double voltage)
{
  Xmin = Xm;
  Xmax = XM;
  Ymin = Ym;
  Ymax = YM;
  Zmin = Zm;
  Zmax = ZM;
  V = voltage;
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
