#include <G4GDMLParser.hh>

#include "boundary.h"

boundary::boundary()
{
  // make the field cage
  double height = 100;
  double width = 100;
  double length = 100;
  double wall_thickness = 2;
  
  volumes[0] = new volume(0., length,
			  0., width,
			  0., wall_thickness,
			  constant(0.));
  volumes[1] = new volume(0., length,
			  0., width,
			  height - wall_thickness, height,
			  constant(100.));
  volumes[2] = new volume(0., wall_thickness,
			  0., width,
			  0., height,
			  linear(0., 0., 0., 1.));
  volumes[3] = new volume(length - wall_thickness, length,
			  0., width,
			  0., height,
			  linear(0., 0., 0., 1.));
  volumes[4] = new volume(0., length,
			  0., wall_thickness,
			  0., height,
			  linear(0., 0., 0., 1.));
  volumes[5] = new volume(0., length,
			  width - wall_thickness, width,
			  0., height,
			  linear(0., 0., 0., 1.));

  nVolumes = 6;
  
  // make the pads
  double padSize = 9;
  double minSpacing = 5;
  int nPadsPerRow = (length - minSpacing)/(padSize + minSpacing);
  double spacing = (length - nPadsPerRow*padSize)/(nPadsPerRow + 1);
  for ( int i = 0; i < nPadsPerRow; i++ ) {
    for ( int j = 0; j < nPadsPerRow; j++ ) {
      volumes[6 + nPadsPerRow*j + i] = new volume(spacing + (padSize + spacing)*i,
							    (spacing + padSize) + (padSize + spacing)*i,
							    spacing + (padSize + spacing)*j,
							    (spacing + padSize) + (padSize + spacing)*j,
							    97.,
							    98.,
							    constant(200.));
      nVolumes++;
    }
  }
    
  for ( uint i = 0; i < nVolumes; i++ ) {
    if ( volumes[i] -> Xmin < Xmin ) {
      Xmin = volumes[i] -> Xmin;
    }
    if ( volumes[i] -> Xmax > Xmax ) {
      Xmax = volumes[i] -> Xmax;
    }
    if ( volumes[i] -> Ymin < Ymin ) {
      Ymin = volumes[i] -> Ymin;
    }
    if ( volumes[i] -> Ymax > Ymax ) {
      Ymax = volumes[i] -> Ymax;
    }
    if ( volumes[i] -> Zmin < Zmin ) {
      Zmin = volumes[i] -> Zmin;
    }
    if ( volumes[i] -> Zmax > Zmax ) {
      Zmax = volumes[i] -> Zmax;
    }
  }

  std::cout << Xmin << '\t' << Xmax << '\n'
	    << Ymin << '\t' << Ymax << '\n'
	    << Zmin << '\t' << Zmax << '\n'
	    << std::endl;
}

bool boundary::is_in_boundary(double x, double y, double z)
{
  bool is_in_any = false;
  for ( uint i = 0; i < nVolumes; i++ ) {
    if ( volumes[i] -> is_in_boundary(x, y, z) ) {
      is_in_any = true;
    }
  }
  return is_in_any;
}

double boundary::boundary_value(double x, double y, double z)
{
  for ( uint i = 0; i < nVolumes; i++ ) {
    if ( volumes[i] -> is_in_boundary(x, y, z) ) {
      return volumes[i] -> V(x, y, z);
    }
  }
  return 0;
}
