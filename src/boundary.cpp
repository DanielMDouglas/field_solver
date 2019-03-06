#include <G4GDMLParser.hh>

#include "boundary.h"

boundary::boundary()
{
  // make the field cage
  double xLow = -1.05;
  double xHigh = 1.05;
  double yLow = -1.05;
  double yHigh = 1.05;
  double zLow = -0.3;
  double zHigh = 1.8;
  
  double wall_thickness = 0.05;
  
  // nVolumes = 0;

  // far cathode plane  
  volumes[nVolumes] = new volume(xLow, xHigh,
				 yLow, yHigh,
				 zHigh - wall_thickness, zHigh,
				 constant(-0.3));
  nVolumes++;
  
  // make_field_cage(xLow, xHigh,
  // 		  yLow, yHigh,
  // 		  1.0, zHigh,
  // 		  wall_thickness,
  // 		  0.23, -0.273);
  
  make_wires(xLow, xHigh,
  	     yLow, yHigh,
  	     zLow, zHigh);  

  // make_pads(xLow, xHigh,
  // 	    yLow, yHigh,
  // 	    zLow, zHigh);
  
  for ( int i = 0; i < nVolumes; i++ ) {
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
}

void boundary::make_pads(double xLow, double xHigh,
			 double yLow, double yHigh,
			 double zLow, double zHigh)
{
  double padSize = 0.2;
  double padThickness = 0.05;
  double minSpacing = 0.1;
  int nPadsPerRow = (xHigh - xLow - minSpacing)/(padSize + minSpacing);
  double spacing = (xHigh - xLow - nPadsPerRow*padSize)/(nPadsPerRow + 1);
  double padPotential = 0.23;

  // make the pads
  for ( int i = 0; i < nPadsPerRow; i++ ) {
    for ( int j = 0; j < nPadsPerRow; j++ ) {
      volumes[nVolumes] = new volume(xLow + spacing + (padSize + spacing)*i,
				     (xLow + spacing + padSize) + (padSize + spacing)*i,
				     yLow + spacing + (padSize + spacing)*j,
				     (yLow + spacing + padSize) + (padSize + spacing)*j,
				     -padThickness/2, padThickness/2,
				     constant(padPotential));
      nVolumes++;
    }
  }

  double shaperPotential = 0;
  
  // make the field shaper
  volumes[nVolumes] = new volume(xLow, xHigh,
  				 yLow, yHigh,
  				 -0.025, 0.025,
  				 constant(shaperPotential));
  nVolumes++;
}

void boundary::make_wires(double xLow, double xHigh,
			  double yLow, double yHigh,
			  double zLow, double zHigh)
{
  // wires!
  double wire_rad = 0.025;
  double wire_pitch = 0.3;
  double wire_potential;
  for ( double x = -0.9; x < 0.9; x += wire_pitch ) {
    for ( int row = 0; row < 3; row ++ ) {
      double z = row*wire_pitch;
      if ( row == 0 ) {
	wire_potential = 0.23;
      }
      else if ( row == 1 ) {
	wire_potential = 0;
      }
      else if ( row == 2 ) {
	wire_potential = -0.11;
      }
      
      volumes[nVolumes] = new volume(x - wire_rad, x + wire_rad,
				     yLow, yHigh,
				     z - wire_rad, z + wire_rad,
				     constant(wire_potential));
      nVolumes++;
    }
  }
}

void boundary::make_field_cage(double xLow, double xHigh,
			       double yLow, double yHigh,
			       double zLow, double zHigh,
			       double wall_thickness,
			       double intercept,
			       double zSlope)
{
  // Field cage walls
  
  volumes[nVolumes] = new volume(xLow, xLow + wall_thickness,
  				 yLow, yHigh,
  				 zLow, zHigh,
  				 linear(intercept, 0., 0., zSlope));
  nVolumes++;
  volumes[nVolumes] = new volume(xHigh - wall_thickness, xHigh,
  				 yLow, yHigh,
  				 zLow, zHigh,
  				 linear(intercept, 0., 0., zSlope));
  nVolumes++;
  volumes[nVolumes] = new volume(xLow, xHigh,
  				 yLow, yLow + wall_thickness,
  				 zLow, zHigh,
  				 linear(intercept, 0., 0., zSlope));
  nVolumes++;
  volumes[nVolumes] = new volume(xLow, xHigh,
  				 yHigh - wall_thickness, yHigh,
  				 zLow, zHigh,
  				 linear(intercept, 0., 0., zSlope));
  nVolumes++;
}

bool boundary::is_in_boundary(double x, double y, double z)
{
  bool is_in_any = false;
  for ( int i = 0; i < nVolumes; i++ ) {
    if ( volumes[i] -> is_in_boundary(x, y, z) ) {
      is_in_any = true;
    }
  }
  return is_in_any;
}

double boundary::boundary_value(double x, double y, double z)
{
  for ( int i = 0; i < nVolumes; i++ ) {
    if ( volumes[i] -> is_in_boundary(x, y, z) ) {
      return volumes[i] -> V(x, y, z);
    }
  }
  return 0;
}
