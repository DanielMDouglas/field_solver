#include <iostream>
#include <string>
#include <vector>

#include <G4GDMLParser.hh>
#include <G4LogicalVolume.hh>

#include "geometry.h"

volume::volume(G4LogicalVolume* vol, double voltage) {
  std::cout << "Initializing volume " << vol -> GetName()
	    << " with potential " << voltage
	    << std::endl;
  log_vol = vol;
  V = voltage;
}

bool volume::is_in_boundary(double x, double y, double z) {
  G4VSolid * solid = log_vol -> GetSolid();
  G4ThreeVector p(x, y, z);
  return ( solid -> Inside(p) == kInside );
}

boundary::boundary(std::string filename) {
  G4GDMLParser parser;
  parser.Read(filename);

  for ( uint i = 0; i < volNames.size(); i++ ) {
    volumes[i] = new volume(parser.GetVolume(volNames[i]), potentials[i]);
  }
}

bool boundary::is_in_boundary(double x, double y, double z) {
  bool is_in_any = false;
  for ( uint i = 0; i < volNames.size(); i++ ) {
    if ( volumes[i] -> is_in_boundary(x, y, z) ) {
      is_in_any = true;
    }
  }
  return is_in_any;
}

double boundary::boundary_value(double x, double y, double z) {
  for ( uint i = 0; i < volNames.size(); i++ ) {
    if ( volumes[i] -> is_in_boundary(x, y, z) ) {
      return volumes[i] -> V;
    }
  }
  std::cout << "should not get here!" << std::endl;
  return 0.;
}
