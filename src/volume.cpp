#include <G4LogicalVolume.hh>
#include <G4VSolid.hh>
#include <G4ThreeVector.hh>
#include <G4VisExtent.hh>

#include "volume.h"

volume::volume(G4LogicalVolume* vol, double voltage, std::vector <double> center)
{
  std::cout << "Initializing volume " << vol -> GetName()
	    << " with potential " << voltage
	    << std::endl;
  log_vol = vol;
  V = voltage;
  cen = center;
}

bool volume::is_in_boundary(double x, double y, double z)
{
  G4VSolid * solid = log_vol -> GetSolid();
  G4ThreeVector p(x + cen[0],
		  y + cen[1],
		  z + cen[2]);
  return ( solid -> Inside(p) == kInside );
}

G4VisExtent volume::extent()
{
  G4VSolid * solid = log_vol -> GetSolid();
  G4VisExtent ext = solid -> GetExtent();

  ext.SetXmin(ext.GetXmin() + cen[0]);
  ext.SetXmax(ext.GetXmax() + cen[0]);
  ext.SetYmin(ext.GetYmin() + cen[1]);
  ext.SetYmax(ext.GetYmax() + cen[1]);
  ext.SetZmin(ext.GetZmin() + cen[2]);
  ext.SetZmax(ext.GetZmax() + cen[2]);

  return ext;
}
