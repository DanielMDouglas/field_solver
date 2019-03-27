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

// void drift_and_save(double xi, scalarField * potential, boundary detector, path * trajectory)
// {
//   drift_path(std::vector <double> {xi, 0., 1.75}, potential, detector);
// }

int main(int argc, char const * argv[])
{
  handleOpts(argc, argv);
  
  scalarField * potential = new scalarField(fieldFileName);
  boundary detector (geom);
  
  const int nPaths = 102;

  // std::thread * threads [nPaths];
  // path * trajectories [nPaths];

  double xi;
  
  std::vector <double> xi_space = linspace(-1.2, 1.2, nPaths);
  for ( int i = 0; i < nPaths; i++ ) {
    xi = xi_space[i];
    path * driftPath;
    drift_path(std::vector <double> {xi, 0, 1.75}, potential, detector, driftPath);

    std::ostringstream stringStream;
    stringStream << "paths/drift_from_" << xi << ".dat";
    driftPath -> print_to_file(stringStream.str());

    delete driftPath;
    
    // threads[i] = new std::thread ( drift_path,
    // 				   std::vector <double> {xi, 0, 1.75},
    // 				   potential,
    // 				   detector,
    // 				   driftPath );
  }

  // for ( int i = 0; i < nPaths; i++ ) {
  //   xi = xi_space[i];
  //   threads[i] -> join();
  //   std::ostringstream stringStream;
  //   stringStream << "paths/drift_from_" << xi << ".dat";
  //   trajectories[i] -> print_to_file(stringStream.str());
  // }
    
  return 0;
}
