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
  
  const int nPaths = 50;

  std::thread * threads [nPaths];
  path * trajectories [nPaths];

  double yi = 0;

  std::vector <double> xi_space = linspace(0., 0.88, nPaths);

  for ( int i = 0; i < xi_space.size(); i++ ) {
    double xi = xi_space[i];
    threads[i] = new std::thread ( drift_path,
				   std::vector <double> {xi, yi, 0.5},
				   potential,
				   detector,
				   std::ref(trajectories[i]) );
  }

  for ( int i = 0; i < nPaths; i++ ) {
    double xi = xi_space[i];
    threads[i] -> join();
    path * thisTraj = trajectories[i];
    std::ostringstream stringStream;
    stringStream << "drift_from_" << xi << ".dat";
    thisTraj -> print_to_file(stringStream.str());
    delete thisTraj;
  }
  
  return 0;
}
