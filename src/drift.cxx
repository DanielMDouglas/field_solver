#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <math.h>
#include <sstream>

#include "drift.h"

std::string fieldFileName = "none";
std::string geom = "bulkPix";

void handleOpts(int argc, char const * argv[])
{
  int opt = 0;
  while ( opt < argc ) {
    std::stringstream optValue;
    std::stringstream argValue;
    optValue << argv[++opt];
    argValue << argv[++opt];
    
    if ( optValue.str() == "-f" ) {
      argValue >> fieldFileName;
    }
    if ( optValue.str() == "-g" ) {
      argValue >> geom;
    }
  }

  if ( fieldFileName == "none" ) {
    std::cout << "Need a field file!" << std::endl;
    exit(1);
  }
  
  std::cout << "Using arguments: \n"
	    << "potential field:  " << fieldFileName << '\n'
	    << "geom:             " << geom << std::endl;
}

// void drift_and_save(double xi, field * potential, boundary detector, path * trajectory)
// {
//   drift_path(std::vector <double> {xi, 0., 1.75}, potential, detector);
// }

int main(int argc, char const * argv[])
{
  handleOpts(argc, argv);
  
  field <double> * potential = new field <double> (fieldFileName);
  boundary detector (geom);
  
  const int nPaths = 100;

  std::thread * threads [nPaths*nPaths];
  path * trajectories [nPaths*nPaths];

  // double xi;
  // double yi;

  std::vector <double> xi_space = linspace(0., 0.195, nPaths);
  std::vector <double> yi_space = linspace(0., 0.195, nPaths);

  std::ofstream outFile ( "final_pos.dat");

  for ( int i = 0; i < xi_space.size(); i++ ) {
    double xi = xi_space[i];
    for ( int j = 0; j < yi_space.size(); j++ ) {
      double yi = yi_space[j];
      threads[i*nPaths + j] = new std::thread ( drift_path,
						std::vector <double> {xi, yi, 1.2},
						potential,
						detector,
						std::ref(trajectories[i*nPaths + j]) );
    }
    for ( int j = 0; j < nPaths; j++ ) {
      threads[i*nPaths + j] -> join();
      std::vector <double> finalPos = trajectories[i*nPaths + j] -> pos.back();
      if ( trajectories[i*nPaths + j] -> fate == "volume" ) {
	outFile << finalPos[0] << '\t' << finalPos[1] << '\t'  << finalPos[2] << '\n';
      }
      delete trajectories[i*nPaths + j];
    }
    // std::ostringstream stringStream;
    // stringStream << "paths/drift_from_" << xi << ".dat";
    // trajectories[i] -> print_to_file(stringStream.str());
    
  }
  outFile.close();
  // for ( int i = 0; i < nPaths*nPaths; i++ ) {
  //   std::vector <double> finalPos = trajectories[i] -> pos.back();
  //   if ( trajectories[i] -> fate == "volume" ) {
  //     std::cout << finalPos[0] << '\t' << finalPos[1] << '\t'  << finalPos[2] << std::endl;
  //   }
  // }

  return 0;
}
