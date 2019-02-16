#include <iostream>
#include <string>
#include <vector>

#include <G4GDMLParser.hh>
#include <G4LogicalVolume.hh>
#include <G4VisExtent.hh>
#include <G4DisplacedSolid.hh>
#include <G4AffineTransform.hh>

#include <G4VoxelLimits.hh>
// #include <G4double.hh>


#include "geometry.h"

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
  G4DisplacedSolid * disp_solid = solid -> GetDisplacedSolidPtr();
  G4VoxelLimits * pVoxelLimit = new G4VoxelLimits();
  G4AffineTransform * pTransform = new G4AffineTransform();
  G4double pMin;
  G4double pMax;
  
  return solid -> GetExtent();
}

// void volume::transform()
// {
//   G4VSolid * solid = log_vol -> GetSolid();
//   G4DisplacedSolid * disp_solid = solid -> GetDisplacedSolidPtr();
//   G4AffineTransform transf =  disp_solid -> GetTransform();
// }

boundary::boundary(std::string filename)
{
  G4GDMLParser parser;
  parser.Read(filename);

  for ( uint i = 0; i < volNames.size(); i++ ) {
    volumes[i] = new volume(parser.GetVolume(volNames[i]), potentials[i]);
    // std::cout << parser.GetPosition(posNames[i]) << std::endl;
  }
}

bool boundary::is_in_boundary(double x, double y, double z)
{
  bool is_in_any = false;
  for ( uint i = 0; i < volNames.size(); i++ ) {
    if ( volumes[i] -> is_in_boundary(x, y, z) ) {
      is_in_any = true;
    }
  }
  return is_in_any;
}

double boundary::boundary_value(double x, double y, double z)
{
  for ( uint i = 0; i < volNames.size(); i++ ) {
    if ( volumes[i] -> is_in_boundary(x, y, z) ) {
      return volumes[i] -> V;
    }
  }
  std::cout << "should not get here!" << std::endl;
  return 0.;
}

G4VisExtent boundary::extent()
{
  G4VisExtent extent;
  std::cout << extent.GetXmin() << std::endl;

  for ( uint i = 0; i < volNames.size(); i++ ) {
    G4VisExtent this_extent = volumes[i] -> extent();
    
    if ( this_extent.GetXmin() < extent.GetXmin() ) {
      extent.SetXmin(this_extent.GetXmin());
    }
    if ( this_extent.GetXmax() > extent.GetXmax() ) {
      extent.SetXmax(this_extent.GetXmax());
    }
    if ( this_extent.GetYmin() < extent.GetYmin() ) {
      extent.SetYmin(this_extent.GetYmin());
    }
    if ( this_extent.GetYmax() > extent.GetYmax() ) {
      extent.SetYmax(this_extent.GetYmax());
    }
    if ( this_extent.GetZmin() < extent.GetZmin() ) {
      extent.SetZmin(this_extent.GetZmin());
    }
    if ( this_extent.GetZmax() > extent.GetZmax() ) {
      extent.SetZmax(this_extent.GetZmax());
    }
  }
  return extent;
}
