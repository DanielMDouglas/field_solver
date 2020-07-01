#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <math.h>
#include <sstream>

#include "drift.h"
#include "argParser.h"

std::string fieldFileName = "none";
std::string geom = "geometries/bulkPix.json";

void handleOpts(int argc, const char ** argv)
{
  argParser parser;

  parser.add_option("-f", [](arg_t ss) {*ss >> fieldFileName;});
  parser.add_option("-g", [](arg_t ss) {*ss >> geom;});

  parser.parse(argc, argv);
  
  if ( fieldFileName == "none" ) {
    std::cout << "Need a field file!" << std::endl;
    exit(1);
  }
}

void saySettings()
{
  std::cout << "Using arguments: \n"
	    << "potential field:  " << fieldFileName << '\n'
	    << "geom:             " << geom << std::endl;
}

// void drift_and_save(double xi, field * potential, boundary detector, path * trajectory)
// {
//   drift_path(std::vector <double> {xi, 0., 1.75}, potential, detector);
// }

int main(int argc, const char ** argv)
{
  handleOpts(argc, argv);
  
  field <double> * potential = new field <double> (fieldFileName);
  boundary detector (geom);
  
  const int nPaths = 25;

  std::thread * threads [nPaths];
  path * trajectories [nPaths];

  double yi = 0.;

  std::vector <double> xi_space = linspace(0.0, 0.19, nPaths);

  for ( int i = 0; i < xi_space.size(); i++ ) {
    double xi = xi_space[i];
    threads[i] = new std::thread ( drift_path,
				   std::vector <double> {xi, yi, 0.45},
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
