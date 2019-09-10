#include <iostream>
#include <vector>
#include <TRandom3.h>

#include "boundary.h"
#include "pad.h"

boundary::boundary(std::string which)
{
  if ( which == "bulkPix" ) {
    make_bulkPix();
  }
  if ( which == "box_uneven_res" ) {
    make_box_uneven_res();
  }
  if ( which == "bulkPixSingle" ) {
    make_bulkPix_single();
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
  else if ( which == "linear_cond" ) {
    make_linear_cond_defect();
  }
  else if ( which == "linear_diel" ) {
    make_linear_diel_defect();
  }
  else if ( which == "sheet" ) {
    make_sheet();
  }
  else if ( which == "sheet_cond" ) {
    make_sheet_cond_defect();
  }
  else if ( which == "sheet_rand_cond" ) {
    make_sheet_random_cond_defect();
  }
  else if ( which == "sheet_diel" ) {
    make_sheet_diel_defect();
  }
  else if ( which == "cap" ) {
    make_capacitor();
  }
  else if ( which == "cap_die" ) {
    make_capacitor_with_dielectric();
  }
  else if ( which == "test" ) {
    make_test();
  }
}

void boundary::make_test()
{
  Xmin = 0;
  Xmax = 1;
  Ymin = 0;
  Ymax = 1;
  Zmin = 0;
  Zmax = 1;

  double wall_thickness = 0.2;
  
  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymax,
			Zmin, Zmin + wall_thickness,
			-1));

  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymax,
			Zmax - wall_thickness, Zmax,
			constant(0)));

  add_volume(new volume(Xmin, Xmin + wall_thickness,
			Ymin, Ymax,
			Zmin, Zmax,
			0));

  add_volume(new volume(Xmax - wall_thickness, Xmax,
			Ymin, Ymax,
			Zmin, Zmax,
			0));

  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymin + wall_thickness,
			Zmin, Zmax,
			0));

  add_volume(new volume(Xmin, Xmax,
			Ymax - wall_thickness, Ymax,
			Zmin, Zmax,
			0));
}

void boundary::make_box_uneven_res()
{
  // periodicX = true;
  periodicY = true;

  Xmin = 0;
  Xmax = 3;
  Ymin = 0;
  Ymax = 0.5;
  Zmin = 0;
  Zmax = 5;

  double wall_thickness = 0.15;
  double nom_res = 1.e12;
  double percent_var = 0.1;
  
  // bottom wall
  // "anode"
  add_volume(new volume(Xmin, Xmax,
  			Ymin, Ymax,
  			Zmin, Zmin + wall_thickness,
  			// -0.1));
  			constant(0.)));

  // top wall
  // "cathode"
  // but really this a von Neumann b.c.
  // so let's place it slightly outside of the volume
  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymax,
			Zmax - wall_thickness, Zmax,
			// -0.5));
			constant(2.5)));

  // left wall
  // -10% resistivity variation
  add_volume(new volume(Xmin, Xmin + wall_thickness,
  			Ymin, Ymax,
  			Zmin, Zmax,
  			// constant(1.)));
  			constant(3.4),
  			// constant(1./(nom_res))));
			linear(0.5/(nom_res), 0, 0, 0.5/(nom_res)/2.5)));
			
  // right wall
  // +10% resistivity variation
  // add_volume(new volume(Xmax - wall_thickness, Xmax,
  // 			Ymin, Ymax,
  // 			Zmin, Zmax - wall_thickness,
  // 			// constant(1.)));
  // 			constant(3.4),
  // 			constant(1./(nom_res))));

  // back wall
  // nominal resistivity
  // add_volume(new volume(Xmin, Xmax,
  // 			Ymin, Ymin + wall_thickness,
  // 			Zmin, Zmax - wall_thickness,
  // 			// constant(1.)));
  // 			constant(3.4),
  // 			constant(1./(nom_res))));

  // front wall
  // nominal variation
  // add_volume(new volume(Xmin, Xmax,
  // 			Ymax - wall_thickness, Ymax,
  // 			Zmin, Zmax - wall_thickness,
  // 			// constant(1.)));
  // 			constant(3.4),
  // 			constant(1./(nom_res))));

  // // add a test block with really high permittivity...
  // add_volume(new volume(0.4, 0.6,
  // 			0.4, 0.6,
  // 			0.4, 0.6,
  // 			constant(100),
  // 			constant(0)));

  // main volume
  // Liquid Argon
  add_volume(new volume(Xmin, Xmax,
  			Ymin, Ymax,
  			Zmin, Zmax,
  			constant(1.504),
  			constant(0)));
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

  // main volume, no irregularities
  add_volume(new volume(Xmin, Xmax,
  			Ymin, Ymax,
  			Zmin, Zmax,
  			constant(1),
  			constant(1)));
}

void boundary::make_linear_cond_defect()
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

  // main volume, conductivity variation
  add_volume(new volume(Xmin, Xmax,
  			Ymin, Ymax,
  			Zmin, Zmax,
  			constant(1),
  			gaussian((Xmax - Xmin)/2,
  				 (Ymax - Ymin)/2,
  				 (Zmax - Zmin)/2,
  				 0.1, 1, -0.75)));
}

void boundary::make_linear_diel_defect()
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

  // main volume, permittivity variation
  add_volume(new volume(Xmin, Xmax,
  			Ymin, Ymax,
  			Zmin, Zmax,
  			gaussian((Xmax - Xmin)/2,
  				 (Ymax - Ymin)/2,
  				 (Zmax - Zmin)/2,
  				 0.1, 1, 5),
  			constant(1)));
}

void boundary::make_sheet()
{
  periodicX = true;
  periodicY = true;
  
  Xmin = 0;
  Xmax = 1;
  Ymin = 0;
  Ymax = 0.03;
  Zmin = 0;
  Zmax = 1;

  double wall_thickness = 0.01;

  // // low voltage (0) at z = 0
  // add_volume(new volume(Xmin, Xmax,
  // 			Ymin, Ymax,
  // 			Zmin, Zmin + wall_thickness,
  // 			constant(0.)));

  // TEST: VN boundary below, which defaults to V = 0 for now
  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymax,
			Zmin, Zmin + wall_thickness,
			-1));

  // high voltage (1) at z = 1
  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymax,
			Zmax - wall_thickness, Zmax,
			constant(1.)));

  // main volume, no irregularities
  add_volume(new volume(Xmin, Xmax,
  			Ymin, Ymax,
  			Zmin, Zmax,
  			constant(1),
  			constant(0)));
}

void boundary::make_sheet_cond_defect()
{
  periodicX = true;
  periodicY = true;
  
  Xmin = 0;
  Xmax = 1;
  Ymin = 0;
  Ymax = 0.03;
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

  // main volume, conductivity variation
  add_volume(new volume(Xmin, Xmax,
  			Ymin, Ymax,
  			Zmin, Zmax,
  			constant(1),
  			gaussian((Xmax - Xmin)/2,
  				 (Ymax - Ymin)/2,
  				 (Zmax - Zmin)/2,
  				 0.1, 1, -1)));
}

void boundary::make_sheet_random_cond_defect()
{
  // periodicX = true;
  periodicY = true;
  
  Xmin = 0;
  Xmax = 3;
  Ymin = 0;
  Ymax = 0.5;  
  Zmin = 0;
  Zmax = 5;

  double wall_thickness = 0.15;
  double nom_res = 1.e12;
  double percent_var = 0.1;
  
  // low voltage (0) at z = 0
  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymax,
			Zmin, Zmin + wall_thickness,
			constant(0.)));

  // high voltage (1) at z = 1
  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymax,
			Zmax - wall_thickness, Zmax,
			constant(2.5)));

  // main volume, conductivity variation
  double nomCond = 1.e-12; // not accurate
  double nomDiel = 3.4;
  int Ndefects = 500;
  double stdDefectSize = 5.e-13; // this determines the magnitude
  double stdDefectWidth = 0.1; // this determines the scale

  std::vector <std::function <double (double, double, double)>> defectList = {};

  TRandom3 * rng = new TRandom3();
  for ( int i = 0; i < Ndefects; i++ ) { 
    double height = rng -> Gaus(); // drawn from a normal distribution
    double locX = rng -> Uniform(Xmin, Xmax); // drawn from a uniform distribution
    double locY = (Ymax - Ymin)/2; // fixed
    double locZ = rng -> Uniform(Zmin, Zmax); // drawn from a uniform distribution

    defectList.push_back(gaussian(locX,
				  locY,
				  locZ,
				  stdDefectWidth, 0, height));
  }
  // sum up all of the defects
  // mean should be 0
  std::function <double (double, double, double)> cond_func = func_sum(defectList);
  // find the standard deviation from a statistical sample
  // maybe put this into a function?
  int Nsamples = 100000;
  // double mean = 0;
  double std = 0;
  for ( int i = 0; i < Nsamples; i++ ) {
    double x = rng -> Uniform(Xmin, Xmax);
    double y = rng -> Uniform(Ymin, Ymax);
    double z = rng -> Uniform(Zmin, Zmax);

    // mean += cond_func(x, y, z);
    std += pow(cond_func(x, y, z), 2);
  }
  // mean /= Nsamples;
  std /= (Nsamples - 1);
  std = pow(std, 0.5);
  // std::cout << "mean of the defects is " << mean << " (should be 0)" << std::endl;
  // std::cout << "std of the defects is " << std << " (should be something?)" << std::endl;

  cond_func = func_mul({cond_func, constant(stdDefectSize/(std*nomCond))});
  cond_func = func_sum({cond_func, constant(1.)});
  cond_func = func_mul({cond_func, constant(nomCond)}); 
  cond_func = func_rect(cond_func);
  
  // mean = 0;
  // std = 0;
  // for ( int i = 0; i < Nsamples; i++ ) {
  //   double x = rng -> Uniform(Xmin, Xmax);
  //   double y = rng -> Uniform(Ymin, Ymax);
  //   double z = rng -> Uniform(Zmin, Zmax);

  //   mean += cond_func(x, y, z);
  //   std += pow(cond_func(x, y, z) - nomCond, 2);
  // }
  // mean /= Nsamples;
  // std /= (Nsamples - 1);
  // std = pow(std, 0.5);
  // std::cout << "mean of the defects is " << mean << " (should be 0)" << std::endl;
  // std::cout << "std of the defects is " << std << " (should be something?)" << std::endl;

  std::function <double (double, double, double)> diel_func = constant(nomDiel);
  
  add_volume(new volume(Xmin, Xmin + wall_thickness,
  			Ymin, Ymax,
  			Zmin, Zmax,
  			diel_func,
			cond_func));

  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymax,
			Zmin, Zmax,
			constant(1.504),
			constant(0)));
}

void boundary::make_sheet_diel_defect()
{
  periodicX = true;
  periodicY = true;
  
  Xmin = 0;
  Xmax = 1;
  Ymin = 0;
  Ymax = 0.03;
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

  // main volume, permittivity variation
  add_volume(new volume(Xmin, Xmax,
  			Ymin, Ymax,
  			Zmin, Zmax,
  			gaussian((Xmax - Xmin)/2,
  				 (Ymax - Ymin)/2,
  				 (Zmax - Zmin)/2,
  				 0.1, 1, 5),
  			constant(1)));
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

void boundary::make_capacitor_with_dielectric()
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
  // dielectric in the gap
  add_volume(new volume(-plate_width/4, plate_width/4,
			-plate_width/4, plate_width/4,
			-plate_separation/2, plate_separation/2,
			4., 0));
  // material everywhere else
  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymax,
			Zmin, Zmax,
			1., 0));
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

  // // PCB is a dielectric with ep_r ~ 4.35
  // add_volume(new volume(Xmin, Xmax,
  // 			Ymin, Ymax,
  // 			-padThickness/2, padThickness/2,
  // 			4.35));

  // backstop
  add_volume(new volume(Xmin, Xmax,
  			Ymin, Ymax,
  			Zmin, Zmin + wall_thickness,
  			constant(-0.1)));
}

void boundary::make_bulkPix_single()
{ 
  // single pad version of bulkPix
  
  periodicX = true;
  periodicY = true;
    
  Xmin = -0.2;
  Xmax = 0.2;
  Ymin = -0.2;
  Ymax = 0.2;
  Zmin = -0.2;
  Zmax = 1.4;
  
  double wall_thickness = 0.01;

  // want the e field to be 500V/cm,
  // so "cathode" potential should be
  // roughly 500V*Zmax
  double cathode_potential = -0.5*Zmax;
  // far cathode plane  
  add_volume(new volume(Xmin, Xmax,
			Ymin, Ymax,
			Zmax - wall_thickness, Zmax,
			constant(cathode_potential)));
  
  double padSize = 0.2;
  double padThickness = 0.05;
  double padPotential = 0.;

  // make the pad
  add_volume(new volume(-padSize/2, padSize/2,
			-padSize/2, padSize/2,
			-padThickness/2, padThickness/2,
			constant(padPotential)));

  // PCB is a dielectric with ep_r ~ 4.35
  // and volume resistivity ~ 1.e8 MOhm cm
  add_volume(new volume(Xmin, Xmax,
  			Ymin, Ymax,
  			-padThickness/2, padThickness/2,
  			4.35, 1.e-14));

  // backstop
  add_volume(new volume(Xmin, Xmax,
  			Ymin, Ymax,
  			Zmin, Zmin + wall_thickness,
  			constant(-0.1)));
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

void boundary::add_volume(volume * newVol)
{
  volumes[nVolumes] = newVol;
  nVolumes++;
}

bool boundary::is_in_conductor(double x, double y, double z)
{
  bool is_in_any = false;
  for ( int i = 0; i < nVolumes; i++ ) {
    if ( ( volumes[i] -> type == "conductor" )
	 and ( volumes[i] -> is_in_boundary(x, y, z) ) ) {
      is_in_any = true;
    }
  }
  return is_in_any;
}

bool boundary::is_in_VN(double x, double y, double z)
{
  bool is_in_any = false;
  for ( int i = 0; i < nVolumes; i++ ) {
    if ( ( volumes[i] -> type == "VN" )
	 and ( volumes[i] -> is_in_boundary(x, y, z) ) ) {
      is_in_any = true;
    }
  }
  return is_in_any;
}

bool boundary::is_in_volume(double x, double y, double z)
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
    if ( ( volumes[i] -> type == "conductor" )
	 and ( volumes[i] -> is_in_boundary(x, y, z) ) ) {
      return volumes[i] -> V(x, y, z);
    }
  }
  return 0;
}

double boundary::Efield(double x, double y, double z)
{
  for ( int i = 0; i < nVolumes; i++ ) {
    if ( ( volumes[i] -> type == "VN" )
	 and ( volumes[i] -> is_in_boundary(x, y, z) ) ) {
      return volumes[i] -> Efield;
    }
  }
  return 0;
}

double boundary::permittivity(double x, double y, double z)
{
  for ( int i = 0; i < nVolumes; i++ ) {
    if ( volumes[i] -> is_in_boundary(x, y, z) ) {
      return volumes[i] -> er(x, y, z);
    }
  }
  return 0;
}

double boundary::conductivity(double x, double y, double z)
{
  for ( int i = 0; i < nVolumes; i++ ) {
    if ( volumes[i] -> is_in_boundary(x, y, z) ) {
      return volumes[i] -> sigma(x, y, z);
    }
  }
  return 0;
}
