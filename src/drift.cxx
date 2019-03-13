#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <math.h>
#include <sstream>

#include "drift.h"

const double dt = 1.e-3; // us

const int max_iter = 1e5;
const double ds = 0.01; // differentiation distance

std::vector <double> E(std::vector <double> pos, scalarField V)
{
  return -1*V.interpolate_grad(pos);
}

std::vector <double> driftV(std::vector <double> eField, double temperature)
{
  // Drift Velocity as a function of Electric Field and LAr Temperature
  // from : W. Walkowiak, NIM A 449 (2000) 288-294
  //
  // Efield should have units of kV/cm
  // Temperature should have units of Kelvin

  double magE = mag(eField);
  std::vector <double> dir = norm(eField);
  if ( magE > 4.0 ) {
    magE = 4.0;
  }

  if ( temperature < 87.0 || temperature > 94.0 ) {
    std::cout << "DriftVelocity Warning!: Temperature value of "
	      << temperature
	      << " K is outside of range covered by drift velocity"
	      << " parameterization.  Returned value may not be"
	      << " correct"
	      << std::endl;
  }

  double tShift = -87.203 + temperature;
  double xFit = 0.0938163 - 0.0052563*tShift - 0.0001470*tShift*tShift;
  double uFit = 5.18406 + 0.01448*tShift - 0.003497*tShift*tShift - 0.000516*tShift*tShift*tShift;
  double vd;
  
  // Icarus Parameter set
  // Use as default
  double P1 = -0.04640; // K^-1
  double P2 = 0.01712; // K^-1
  double P3 = 1.88125; // (kV/cm)^-1
  double P4 = 0.99408; // kV/cm
  double P5 = 0.01172; // (kV/cm)^-P6
  double P6 = 4.20214; 
  double T0 = 105.749; // K

  // Walkowiak Parameter set
  double P1W = -0.01481; // K^-1
  double P2W = 0.0075; // K^-1
  double P3W = 0.141; // (kV/cm)^-1
  double P4W = 12.4; // kV/cm
  double P5W = 1.627; // (kV/cm)^-P6
  double P6W = 0.317; 
  double T0W = 90.371; // K

  // From Craig Thorne . . . currently not documented
  // smooth transition from linear at small fields to
  // icarus fit at most fields to Walkowiak at very high fields
  if ( magE < xFit ) vd = magE*uFit;
  else if ( magE < 0.619 ) {
    vd = ((P1*(temperature-T0)+1)
	  *(P3*magE*log(1+P4/magE) + P5*pow(magE, P6))
	  +P2*(temperature-T0));
  }
  else if ( magE < 0.699 ) {
    vd = 12.5*(magE - 0.619)*((P1W*(temperature-T0W)+1)
			      *(P3W*magE*log(1+P4W/magE) + P5W*pow(magE, P6W))
			      +P2W*(temperature-T0W)) +
      12.5*(0.699 - magE)*((P1*(temperature-T0)+1)
			   *(P3*magE*log(1+P4/magE) + P5*pow(magE, P6))
			   +P2*(temperature-T0));
  }
  else {
    vd = ((P1W*(temperature-T0W)+1)
	  *(P3W*magE*log(1+P4W/magE) + P5W*pow(magE, P6W))
	  +P2W*(temperature-T0W));
  }

  vd /= 10;

  return -vd*dir; // cm/us
}

path drift_path(std::vector <double> init_pos, scalarField V, boundary b)
{
  double t = 0;
  path trajectory = path(dt);
  std::vector <double> pos = init_pos; // cm
  
  std::cout << "Initial pos: " << "(" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
  
  while ( true ) {
    trajectory.steps.push_back(pos);
    if ( (pos[0] < b.Xmin) or (pos[0] > b.Xmax) or
	 (pos[1] < b.Ymin) or (pos[1] > b.Ymax) or
	 (pos[2] < b.Zmin) or (pos[2] > b.Zmax) ) {
      std::cout << "particle left boundary!" << std::endl;
      trajectory.fate = "OOB";
      break;
    }
    else if ( t >= dt*max_iter ) {
      std::cout << "particle timed out!" << std::endl;
      trajectory.fate = "OOT";
      break;
    }
    else if ( b.is_in_boundary(pos[0], pos[1], pos[2]) ) {
      std::cout << "particle terminated in a volume!" << std::endl;
      trajectory.fate = "volume";
      break;
    }
    else {
      // assume T = boiling point of Ar for now
      pos = pos + driftV(E(pos, V), 87.302)*dt;
      t += dt;
    }
  }
  
  trajectory.arrivalTime = t;

  std::cout << "Final pos: " << "(" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
  std::cout << "Arrival time: " << t << std::endl;
  std::cout << "Fate: " << trajectory.fate << std::endl;
  std::cout << std::endl;

  return trajectory;
}

void drift_and_save(double xi, scalarField potential, boundary detector)
{
  path trajectory = drift_path(std::vector <double> {xi, 0.15, 0.8}, potential, detector);
  std::ostringstream stringStream;
  stringStream << "paths/drift_from_" << xi << ".dat";
  trajectory.print_to_file(stringStream.str());
}

int main()
{
  scalarField potential = scalarField("wire_cathode_300V.dat");
  boundary detector;
  detector.Zmin = -0.2;
  
  const int nPaths = 50;

  // std::vector <std::thread> threads {nPaths};
  std::thread * threads [nPaths];
  
  std::vector <double> xi_space = linspace(-0.9, 0.9, nPaths);
  // for ( double xi: linspace(-0.9, 0.9, nPaths) ) {
  for ( int i = 0; i < nPaths; i++ ) {
    double xi = xi_space[i];
    threads[i] = new std::thread (drift_and_save, xi, potential, detector);
  }

  for ( int i = 0; i < nPaths; i++ ) {
    threads[i] -> join();
  }
    
  return 0;
}
