#include <G4GDMLParser.hh>

#include <boundary.h>

boundary::boundary(std::string filename)
{
  G4GDMLParser parser;
  parser.Read(filename);

  for ( uint i = 0; i < volNames.size(); i++ ) {
    volumes[i] = new volume(parser.GetVolume(volNames[i]), potentials[i]);
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
