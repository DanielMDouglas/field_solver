#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <math.h>
#include <sstream>

#include "drift.h"

const double dt = 1.e-4; // us
const int max_iter = 1e5;

path drift_path(std::vector <double> init_pos, scalarField * V, boundary b)
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

void drift_and_save(double xi, scalarField * potential, boundary detector)
{
  path trajectory = drift_path(std::vector <double> {xi, 0.15, 0.8}, potential, detector);
  std::ostringstream stringStream;
  stringStream << "paths/drift_from_" << xi << ".dat";
  trajectory.print_to_file(stringStream.str());
}

int main()
{
  // scalarField potential = scalarField("wire_cathode_300V.dat");
  scalarField * potential = new scalarField("final.dat");
  boundary detector ("bulkPix");
  
  const int nPaths = 50;

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
