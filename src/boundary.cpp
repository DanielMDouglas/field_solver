#include <iostream>

#include "boundary.h"

boundary::boundary(std::string which)
{
  if ( which == "bulkPix" ) {
    make_bulkPix();
  }
  else if ( which == "bulkPixWeighting" ) {
    make_bulkPixWeighting();
  }
  else if ( which == "bulkWires" ) {
    make_bulkWires();
  }
}

void boundary::make_bulkPix()
{
  periodicX = true;
  periodicY = true;
    
  double xLow = -1.;
  double xHigh = 1.;
  double yLow = -1.;
  double yHigh = 1.;
  double zLow = -0.2;
  double zHigh = 1.8;
  
  double wall_thickness = 0.01;
  
  // far cathode plane  
  volumes[nVolumes] = new volume(xLow, xHigh,
				 yLow, yHigh,
				 zHigh - wall_thickness, zHigh,
				 constant(-0.6));
  nVolumes++;
  
  double padSize = 0.2;
  double padThickness = 0.05;
  double spacing = 0.2;
  int nPadsPerRow = (xHigh - xLow)/(padSize + spacing);
  double padPotential = 0.;

  std::cout << "Setting up pads with size "
	    << padSize << " and " << padSize + spacing
	    << " pitch!" << std::endl;
  
  // make the pads
  for ( int i = 0; i < nPadsPerRow; i++ ) {
    for ( int j = 0; j < nPadsPerRow; j++ ) {
      volumes[nVolumes] = new volume(xLow + spacing/2 + (padSize + spacing)*i,
				     (xLow + spacing/2 + padSize) + (padSize + spacing)*i,
				     yLow + spacing/2 + (padSize + spacing)*j,
				     (yLow + spacing/2 + padSize) + (padSize + spacing)*j,
				     -padThickness/2, padThickness/2,
				     constant(padPotential));
      volumes[nVolumes] -> isSensitive = true;
      nVolumes++;
    }
  }

  Xmin = xLow;
  Xmax = xHigh;
  Ymin = yLow;
  Ymax = yHigh;
  Zmin = zLow;
  Zmax = zHigh;
}

void boundary::make_bulkPixWeighting()
{
  // periodicX = true;
  // periodicY = true;
    
  double xLow = -1.;
  double xHigh = 1.;
  double yLow = -1.;
  double yHigh = 1.;
  double zLow = -0.2;
  double zHigh = 1.8;
  
  double wall_thickness = 0.01;
  
  // far cathode plane  
  volumes[nVolumes] = new volume(xLow, xHigh,
				 yLow, yHigh,
				 zHigh - wall_thickness, zHigh,
				 constant(0));
  nVolumes++;
  
  double padSize = 0.2;
  double padThickness = 0.05;
  double spacing = 0.2;
  int nPadsPerRow = (xHigh - xLow)/(padSize + spacing);
  double padPotential;

  std::cout << "Setting up pads with size "
	    << padSize << " and " << padSize + spacing
	    << " pitch!" << std::endl;
  
  // make the pads
  for ( int i = 0; i < nPadsPerRow; i++ ) {
    for ( int j = 0; j < nPadsPerRow; j++ ) {
      if ( ( i == 1 ) and ( j == 1 ) ) {
	padPotential = 1;
      }
      else {
	padPotential = 0;
      }
      volumes[nVolumes] = new volume(xLow + spacing/2 + (padSize + spacing)*i,
				     (xLow + spacing/2 + padSize) + (padSize + spacing)*i,
				     yLow + spacing/2 + (padSize + spacing)*j,
				     (yLow + spacing/2 + padSize) + (padSize + spacing)*j,
				     -padThickness/2, padThickness/2,
				     constant(padPotential));
      volumes[nVolumes] -> isSensitive = true;
      nVolumes++;
    }
  }

  Xmin = xLow;
  Xmax = xHigh;
  Ymin = yLow;
  Ymax = yHigh;
  Zmin = zLow;
  Zmax = zHigh;
}

void boundary::make_bulkWires()
{
  periodicX = true;
  periodicY = true;
    
  double xLow = -1.05;
  double xHigh = 1.05;
  double yLow = -1.05;
  double yHigh = 1.05;
  double zLow = -0.25;
  double zHigh = 1.85;
  
  double wall_thickness = 0.01;
  
  // far cathode plane  
  volumes[nVolumes] = new volume(xLow, xHigh,
				 yLow, yHigh,
				 zHigh - wall_thickness, zHigh,
				 constant(-0.4376));
  nVolumes++;

  // backstop plane
  volumes[nVolumes] = new volume(xLow, xHigh,
				 yLow, yHigh,
				 zLow, zLow + wall_thickness,
				 constant(0.22425));
  nVolumes++;

  
  double wire_rad = 0.01;
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

  Xmin = xLow;
  Xmax = xHigh;
  Ymin = yLow;
  Ymax = yHigh;
  Zmin = zLow;
  Zmax = zHigh;
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
