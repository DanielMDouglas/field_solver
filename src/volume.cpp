#include "volume.h"
#include "utils.h"

volume::volume(double Xm = 0, double XM = 0,
	       double Ym = 0, double YM = 0,
	       double Zm = 0, double ZM = 0,
	       double relative_permittivity = 1,
	       double conductivity = 0)
{
  // dielectric constructor
  Xmin = Xm;
  Xmax = XM;
  Ymin = Ym;
  Ymax = YM;
  Zmin = Zm;
  Zmax = ZM;
  V = constant(0);
  type = "dielectric";
  er = constant(relative_permittivity);
  sigma = constant(conductivity);

  center = std::vector <double> {0.5*(Xm + XM),
				 0.5*(Ym + YM),
				 0.5*(Zm + ZM)};
}


volume::volume(double Xm = 0, double XM = 0,
	       double Ym = 0, double YM = 0,
	       double Zm = 0, double ZM = 0,
	       std::function<double (double, double, double)> relative_permittivity = constant(1),
	       std::function<double (double, double, double)> conductivity = constant(0))
{
  // dielectric constructor
  Xmin = Xm;
  Xmax = XM;
  Ymin = Ym;
  Ymax = YM;
  Zmin = Zm;
  Zmax = ZM;
  V = constant(0);
  type = "dielectric";
  er = relative_permittivity;
  sigma = conductivity;

  center = std::vector <double> {0.5*(Xm + XM),
				 0.5*(Ym + YM),
				 0.5*(Zm + ZM)};
}

volume::volume(double Xm = 0, double XM = 0,
	       double Ym = 0, double YM = 0,
	       double Zm = 0, double ZM = 0,
	       std::function<double (double, double, double)> voltage_function = constant(0))
{
  // conductor constructor
  Xmin = Xm;
  Xmax = XM;
  Ymin = Ym;
  Ymax = YM;
  Zmin = Zm;
  Zmax = ZM;
  V = voltage_function;
  type = "conductor";
  // er = constant(99999); // should be infinite, but that would be too many 9's
  // sigma = constant(99999);
  er = constant(1);
  sigma = constant(0);
  
  center = std::vector <double> {0.5*(Xm + XM),
				 0.5*(Ym + YM),
				 0.5*(Zm + ZM)};
}

volume::volume(double Xm = 0, double XM = 0,
	       double Ym = 0, double YM = 0,
	       double Zm = 0, double ZM = 0,
	       double E = 0)
{
  // von Neumann boundary condition constructor
  Xmin = Xm;
  Xmax = XM;
  Ymin = Ym;
  Ymax = YM;
  Zmin = Zm;
  Zmax = ZM;
  V = constant(0);
  Efield = E;
  type = "VN";
  // er = constant(99999); // should be infinite, but that would be too many 9's
  er = constant(0);
  sigma = constant(0);
  
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
