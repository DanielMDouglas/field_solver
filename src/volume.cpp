#include <G4LogicalVolume.hh>
#include <G4VSolid.hh>
#include <G4ThreeVector.hh>
#include <G4VisExtent.hh>

#include "volume.h"

volume::volume(G4LogicalVolume* vol, double voltage)
{
  std::cout << "Initializing volume " << vol -> GetName()
	    << " with potential " << voltage
	    << std::endl;
  log_vol = vol;
  V = voltage;
}

bool volume::is_in_boundary(double x, double y, double z)
{
  G4VSolid * solid = log_vol -> GetSolid();
  G4ThreeVector p(x, y, z);
  return ( solid -> Inside(p) == kInside );
}

G4VisExtent volume::extent()
{
  G4VSolid * solid = log_vol -> GetSolid();
  return solid -> GetExtent();
}
