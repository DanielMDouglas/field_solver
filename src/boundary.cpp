#include <iostream>
#include <vector>

#include "boundary.h"
#include "pad.h"

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
  else if ( which == "linear" ) {
    make_linear();
  }
  else if ( which == "cap" ) {
    make_capacitor();
  }
}

void boundary::make_linear()
{
  periodicX = true;
  periodicY = true;
  
  Xmin = 0;
  Xmax = 1;
  Ymin = 0;
  Ymax = 1;
  Zmin = 0;
  Zmax = 1;

  double wall_thickness = 0.01;

  // low voltage (0) at z = 0
  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymax,
			Zmin, Zmin + wall_thickness,
			constant(0.)));

  // high voltage (1) at z = 1
  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymax,
			Zmax - wall_thickness, Zmax,
			constant(1.)));
}

void boundary::make_capacitor()
{
  Xmin = -1;
  Xmax = 1;
  Ymin = -1;
  Ymax = 1;
  Zmin = -1;
  Zmax = 1;

  double wall_thickness = 0.02;
  double plate_separation = 0.4;
  double plate_width = 1.2;
  double plate_thickness = 0.04;
  
  // walls
  add_volume(new volume(Xmin, Xmin + wall_thickness,
			Ymin, Ymax,
			Zmin, Zmax,
			constant(0.)));
  add_volume(new volume(Xmax - wall_thickness, Xmax,
			Ymin, Ymax,
			Zmin, Zmax,
			constant(0.)));
  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymin + wall_thickness,
			Zmin, Zmax,
			constant(0.)));
  add_volume(new volume(Xmin, Xmax,
			Ymax - wall_thickness, Ymax,
			Zmin, Zmax,
			constant(0.)));
  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymax,
			Zmin, Zmin + wall_thickness,
			constant(0.)));
  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymax,
			Zmax - wall_thickness, Zmax,
			constant(0.)));

  // top plate
  add_volume(new volume(-plate_width/2, plate_width/2,
			-plate_width/2, plate_width/2,
			plate_separation - plate_thickness/2, plate_separation + plate_thickness/2,
			constant(1.)));
  // bottom plate
  add_volume(new volume(-plate_width/2, plate_width/2,
			-plate_width/2, plate_width/2,
			-plate_separation - plate_thickness/2, -plate_separation + plate_thickness/2,
			constant(-1.)));
}

void boundary::make_bulkPix()
{
  periodicX = true;
  periodicY = true;
    
  Xmin = -1.4;
  Xmax = 1.4;
  Ymin = -1.4;
  Ymax = 1.4;
  Zmin = -0.2;
  Zmax = 2.6;
  
  double wall_thickness = 0.01;
  
  // far cathode plane  
  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymax,
			Zmax - wall_thickness, Zmax,
			constant(-0.6)));
  
  double padSize = 0.2;
  double padThickness = 0.05;
  double spacing = 0.2;
  // int nPadsPerRow = (Xmax - Xmin)/(padSize + spacing);
  int nPadsPerRow = 7;
  double padPotential = 0.;

  std::cout << "# Setting up pads with size "
	    << padSize << " and " << padSize + spacing
	    << " pitch!" << std::endl;
  
  // make the pads
  for ( int i = 0; i < nPadsPerRow; i++ ) {
    for ( int j = 0; j < nPadsPerRow; j++ ) {
      add_volume(new volume(Xmin + spacing/2 + (padSize + spacing)*i,
			    (Xmin + spacing/2 + padSize) + (padSize + spacing)*i,
			    Ymin + spacing/2 + (padSize + spacing)*j,
			    (Ymin + spacing/2 + padSize) + (padSize + spacing)*j,
			    -padThickness/2, padThickness/2,
			    constant(padPotential)));
    }
  }
}

void boundary::make_bulkPixWeighting()
{
  // periodicX = true;
  // periodicY = true;
    
  Xmin = -1.4;
  Xmax = 1.4;
  Ymin = -1.4;
  Ymax = 1.4;
  Zmin = -0.2;
  Zmax = 2.6;
  
  double wall_thickness = 0.01;
  
  // far cathode plane  
  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymax,
			Zmax - wall_thickness, Zmax,
			constant(0)));
  
  double padSize = 0.2;
  double padThickness = 0.05;
  double spacing = 0.2;
  int nPadsPerRow = 7;
  double padPotential;

  std::cout << "# Setting up pads with size "
	    << padSize << " and " << padSize + spacing
	    << " pitch!" << std::endl;
  
  // make the pads
  for ( int i = 0; i < nPadsPerRow; i++ ) {
    for ( int j = 0; j < nPadsPerRow; j++ ) {
      if ( ( i == 3 ) and ( j == 3 ) ) {
	padPotential = 1;
      }
      else {
	padPotential = 0;
      }
      add_volume(new volume(Xmin + spacing/2 + (padSize + spacing)*i,
			    (Xmin + spacing/2 + padSize) + (padSize + spacing)*i,
			    Ymin + spacing/2 + (padSize + spacing)*j,
			    (Ymin + spacing/2 + padSize) + (padSize + spacing)*j,
			    -padThickness/2, padThickness/2,
			    constant(padPotential)));
    }
  }
}

void boundary::make_bulkWires()
{
  periodicX = true;
  periodicY = true;
    
  Xmin = -1.05;
  Xmax = 1.05;
  Ymin = -1.05;
  Ymax = 1.05;
  Zmin = -0.25;
  Zmax = 1.85;
  
  double wall_thickness = 0.01;
  
  // far cathode plane  
  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymax,
			Zmax - wall_thickness, Zmax,
			constant(-0.4376)));

  // backstop plane
  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymax,
			Zmin, Zmin + wall_thickness,
			constant(0.175)));
  
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
      
      add_volume(new volume(x - wire_rad, x + wire_rad,
			    Ymin, Ymax,
			    z - wire_rad, z + wire_rad,
			    constant(wire_potential)));
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
  
  add_volume(new volume(xLow, xLow + wall_thickness,
			yLow, yHigh,
			zLow, zHigh,
			linear(intercept, 0., 0., zSlope)));
  add_volume(new volume(xHigh - wall_thickness, xHigh,
			yLow, yHigh,
			zLow, zHigh,
			linear(intercept, 0., 0., zSlope)));
  add_volume(new volume(xLow, xHigh,
			yLow, yLow + wall_thickness,
			zLow, zHigh,
			linear(intercept, 0., 0., zSlope)));
  add_volume(new volume(xLow, xHigh,
			yHigh - wall_thickness, yHigh,
			zLow, zHigh,
			linear(intercept, 0., 0., zSlope)));
}

void boundary::add_volume(volume * newVol)
{
  volumes[nVolumes] = newVol;
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

double boundary::permittivity(double x, double y, double z)
{
  for ( int i = 0; i < nVolumes; i++ ) {
    if ( volumes[i] -> is_in_boundary(x, y, z) ) {
      return volumes[i] -> er;
    }
  }
  return 1;
}
