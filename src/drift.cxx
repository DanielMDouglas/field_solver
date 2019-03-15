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
    std::cout << "Need a weighting field file!" << std::endl;
    exit(1);
  }
  
  std::cout << "Using arguments: \n"
	    << "potential field:  " << fieldFileName << '\n'
	    << "geom:             " << geom << std::endl;
}

void drift_and_save(double xi, scalarField * potential, boundary detector)
{
  path trajectory = drift_path(std::vector <double> {xi, 0., 1.75}, potential, detector);
  std::ostringstream stringStream;
  stringStream << "paths/drift_from_" << xi << ".dat";
  trajectory.print_to_file(stringStream.str());
}

int main(int argc, char const * argv[])
{
  handleOpts(argc, argv);
  
  scalarField * potential = new scalarField(fieldFileName);
  boundary detector (geom);
  
  const int nPaths = 100;

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
